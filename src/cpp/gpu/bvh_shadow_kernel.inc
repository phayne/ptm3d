// bvh_shadow_kernel.inc -- the BVH shadow-ray kernel, SINGLE SOURCE.
//
// Included twice by bvh_opencl.cpp: once with PTM_BVH_KERNEL(...) stringizing
// its argument (the OpenCL program text), once expanding it as host C++ under
// small shims (the PTM_BVH_HOST_KERNEL=1 reference path). One text, so the
// host reference cannot drift from what the device compiles.
// Rules for this file: no preprocessor directives and no '//' comments inside
// the macro argument (comments are stripped before stringizing; directives
// would break it). Use block comments.
PTM_BVH_KERNEL(

/* Shadow ray kernel - returns 1 if blocked, 0 if visible */
/* Uses early exit on first hit for maximum performance */
/* Note: Inline functions removed to improve macOS AMD driver compatibility */
__kernel void shadow_ray_kernel(
    /* BVH nodes (SoA) */
    __global const float* node_min_x,
    __global const float* node_min_y,
    __global const float* node_min_z,
    __global const float* node_max_x,
    __global const float* node_max_y,
    __global const float* node_max_z,
    __global const uint* node_left,
    __global const uint* node_right,
    __global const uint* node_tri_start,
    __global const uint* node_tri_count,
    /* Triangle vertices (SoA) */
    __global const float* tri_v0_x,
    __global const float* tri_v0_y,
    __global const float* tri_v0_z,
    __global const float* tri_v1_x,
    __global const float* tri_v1_y,
    __global const float* tri_v1_z,
    __global const float* tri_v2_x,
    __global const float* tri_v2_y,
    __global const float* tri_v2_z,
    /* Triangle indices for leaf nodes */
    __global const uint* tri_indices,
    /* Rays (SoA) */
    __global const float* ray_org_x,
    __global const float* ray_org_y,
    __global const float* ray_org_z,
    __global const float* ray_dir_x,
    __global const float* ray_dir_y,
    __global const float* ray_dir_z,
    __global const float* ray_inv_dir_x,
    __global const float* ray_inv_dir_y,
    __global const float* ray_inv_dir_z,
    /* Per-ray max t along direction (1e10 = effectively infinite) */
    __global const float* ray_tmax_in,
    /* Per-ray min t (the CPU path's tnear rule; 0 for at-face rise-epsilon rays) */
    __global const float* ray_tmin_in,
    /* Per-ray source triangle to ignore (-1 = none): the CPU rise-epsilon */
    /* march skips the origin's own facet by primID. */
    __global const int* ray_skip_prim,
    /* Per-ray rise-epsilon filter: a hit whose out-of-plane rise */
    /* t * (dir . n_src) is below h_eps is self-relief and is ignored. */
    /* h_eps < 0 disables the filter (plain any-hit occlusion). */
    __global const float* ray_dir_dot_n,
    __global const float* ray_h_eps,
    /* Output */
    __global uchar* occluded,
    /* Ray count */
    int num_rays)
{
    int ray_id = get_global_id(0);
    if (ray_id >= num_rays) return;

    /* Load ray data */
    float org_x = ray_org_x[ray_id];
    float org_y = ray_org_y[ray_id];
    float org_z = ray_org_z[ray_id];
    float dir_x = ray_dir_x[ray_id];
    float dir_y = ray_dir_y[ray_id];
    float dir_z = ray_dir_z[ray_id];
    float inv_x = ray_inv_dir_x[ray_id];
    float inv_y = ray_inv_dir_y[ray_id];
    float inv_z = ray_inv_dir_z[ray_id];

    /* Watertight-test ray setup: kz = dominant direction axis; kx, ky the */
    /* others (swapped if dir[kz] < 0 to keep winding); shear constants. */
    float dvec[3];
    dvec[0] = dir_x; dvec[1] = dir_y; dvec[2] = dir_z;
    int kz = 0;
    if (fabs(dvec[1]) > fabs(dvec[kz])) kz = 1;
    if (fabs(dvec[2]) > fabs(dvec[kz])) kz = 2;
    int kx = kz + 1; if (kx == 3) kx = 0;
    int ky = kx + 1; if (ky == 3) ky = 0;
    if (dvec[kz] < 0.0f) { int tmp = kx; kx = ky; ky = tmp; }
    float Sx = dvec[kx] / dvec[kz];
    float Sy = dvec[ky] / dvec[kz];
    float Sz = 1.0f / dvec[kz];

    float ray_tmin = ray_tmin_in[ray_id];
    float ray_tmax = ray_tmax_in[ray_id];
    int   skip_prim = ray_skip_prim[ray_id];
    float dir_dot_n = ray_dir_dot_n[ray_id];
    float h_eps     = ray_h_eps[ray_id];

    /* Stack-based BVH traversal. Depth-first with both children pushed holds */
    /* at most (max depth + 1) entries; the builder caps depth at 32, so 32 */
    /* entries could overflow by one. 64 is provably enough. */
    uint stack[64];
    int stack_ptr = 0;
    stack[stack_ptr++] = 0;  /* Start at root node */

    while (stack_ptr > 0) {
        uint node_idx = stack[--stack_ptr];

        /* === Inline ray-box intersection (slab method) === */
        /* Each slab's exit distance is scaled by 1 + 2*gamma(3) = 1.0000004 */
        /* (Ize 2013, Robust BVH Ray Traversal, JCGT 2(2)), so float rounding */
        /* in the slab arithmetic can only ENLARGE a box, never cull one the */
        /* ray truly enters. Applied in all three slab tests below. */
        float t1, t2;
        float box_tmin = ray_tmin;
        float box_tmax = ray_tmax;

        t1 = (node_min_x[node_idx] - org_x) * inv_x;
        t2 = (node_max_x[node_idx] - org_x) * inv_x;
        box_tmin = max(box_tmin, min(t1, t2));
        box_tmax = min(box_tmax, max(t1, t2) * 1.0000004f);

        t1 = (node_min_y[node_idx] - org_y) * inv_y;
        t2 = (node_max_y[node_idx] - org_y) * inv_y;
        box_tmin = max(box_tmin, min(t1, t2));
        box_tmax = min(box_tmax, max(t1, t2) * 1.0000004f);

        t1 = (node_min_z[node_idx] - org_z) * inv_z;
        t2 = (node_max_z[node_idx] - org_z) * inv_z;
        box_tmin = max(box_tmin, min(t1, t2));
        box_tmax = min(box_tmax, max(t1, t2) * 1.0000004f);

        /* Skip if box not hit */
        if (box_tmin > box_tmax) {
            continue;
        }

        uint tri_count = node_tri_count[node_idx];

        if (tri_count > 0) {
            /* Leaf node - test triangles */
            uint tri_start = node_tri_start[node_idx];

            for (uint i = 0; i < tri_count; i++) {
                uint tri_idx = tri_indices[tri_start + i];

                /* === Watertight ray-triangle test (Woop, Benthin and Wald 2013, */
                /* JCGT 2(1); the algorithm Embree uses). Vertices are taken */
                /* relative to the ray origin and sheared so the ray is the +z */
                /* axis; the 2D edge functions U, V, W of a shared edge are then */
                /* computed from the SAME transformed vertices by both triangles */
                /* and are exact negations of each other (no FMA contraction: */
                /* see the program pragma and the host compile flags), so a ray */
                /* cannot pass between adjacent triangles. Two-sided: a hit needs */
                /* U, V, W of one sign (zeros allowed, so an edge or vertex hit */
                /* counts for both neighbors: conservative, never a gap). */
                float Av[3], Bv[3], Cv[3];
                Av[0] = tri_v0_x[tri_idx] - org_x; Av[1] = tri_v0_y[tri_idx] - org_y; Av[2] = tri_v0_z[tri_idx] - org_z;
                Bv[0] = tri_v1_x[tri_idx] - org_x; Bv[1] = tri_v1_y[tri_idx] - org_y; Bv[2] = tri_v1_z[tri_idx] - org_z;
                Cv[0] = tri_v2_x[tri_idx] - org_x; Cv[1] = tri_v2_y[tri_idx] - org_y; Cv[2] = tri_v2_z[tri_idx] - org_z;

                float Ax = Av[kx] - Sx * Av[kz];
                float Ay = Av[ky] - Sy * Av[kz];
                float Bx = Bv[kx] - Sx * Bv[kz];
                float By = Bv[ky] - Sy * Bv[kz];
                float Cx = Cv[kx] - Sx * Cv[kz];
                float Cy = Cv[ky] - Sy * Cv[kz];

                float pU1 = Cx * By; float pU2 = Cy * Bx; float U = pU1 - pU2;
                float pV1 = Ax * Cy; float pV2 = Ay * Cx; float V = pV1 - pV2;
                float pW1 = Bx * Ay; float pW2 = By * Ax; float W = pW1 - pW2;

                if ((U < 0.0f || V < 0.0f || W < 0.0f) &&
                    (U > 0.0f || V > 0.0f || W > 0.0f)) continue;
                float det = U + V + W;
                if (det == 0.0f) continue;

                float Az = Sz * Av[kz];
                float Bz = Sz * Bv[kz];
                float Cz = Sz * Cv[kz];
                float T = U * Az + V * Bz + W * Cz;
                float hit_t = T / det;

                if (hit_t >= ray_tmin && hit_t <= ray_tmax) {
                    /* Rise-epsilon filter (mirror of the CPU march): skip the */
                    /* source facet and hits rising less than h_eps out of the */
                    /* source plane; any other hit blocks the ray. */
                    if ((int)tri_idx == skip_prim) continue;
                    if (h_eps >= 0.0f && hit_t * dir_dot_n < h_eps) continue;
                    occluded[ray_id] = 1;
                    return;
                }
            }
        } else {
            /* Internal node - intersect both children's AABBs and push in */
            /* front-to-back order (farther first; nearer popped next). */
            /* Shadow rays benefit ~30-40% from this on incoherent / oblique */
            /* sun directions because the early-out fires sooner. */
            uint left_idx  = node_left[node_idx];
            uint right_idx = node_right[node_idx];

            /* Left child slab test */
            float lt1, lt2, l_tmin = ray_tmin, l_tmax = ray_tmax;
            lt1 = (node_min_x[left_idx] - org_x) * inv_x;
            lt2 = (node_max_x[left_idx] - org_x) * inv_x;
            l_tmin = max(l_tmin, min(lt1, lt2));
            l_tmax = min(l_tmax, max(lt1, lt2) * 1.0000004f);
            lt1 = (node_min_y[left_idx] - org_y) * inv_y;
            lt2 = (node_max_y[left_idx] - org_y) * inv_y;
            l_tmin = max(l_tmin, min(lt1, lt2));
            l_tmax = min(l_tmax, max(lt1, lt2) * 1.0000004f);
            lt1 = (node_min_z[left_idx] - org_z) * inv_z;
            lt2 = (node_max_z[left_idx] - org_z) * inv_z;
            l_tmin = max(l_tmin, min(lt1, lt2));
            l_tmax = min(l_tmax, max(lt1, lt2) * 1.0000004f);
            bool l_hit = (l_tmin <= l_tmax) && (left_idx != 0);

            /* Right child slab test */
            float rt1, rt2, r_tmin = ray_tmin, r_tmax = ray_tmax;
            rt1 = (node_min_x[right_idx] - org_x) * inv_x;
            rt2 = (node_max_x[right_idx] - org_x) * inv_x;
            r_tmin = max(r_tmin, min(rt1, rt2));
            r_tmax = min(r_tmax, max(rt1, rt2) * 1.0000004f);
            rt1 = (node_min_y[right_idx] - org_y) * inv_y;
            rt2 = (node_max_y[right_idx] - org_y) * inv_y;
            r_tmin = max(r_tmin, min(rt1, rt2));
            r_tmax = min(r_tmax, max(rt1, rt2) * 1.0000004f);
            rt1 = (node_min_z[right_idx] - org_z) * inv_z;
            rt2 = (node_max_z[right_idx] - org_z) * inv_z;
            r_tmin = max(r_tmin, min(rt1, rt2));
            r_tmax = min(r_tmax, max(rt1, rt2) * 1.0000004f);
            bool r_hit = (r_tmin <= r_tmax) && (right_idx != 0);

            /* Formerly `if (stack_ptr < 30)`: past 30 entries both children */
            /* were silently NOT pushed, dropping whole subtrees and leaking */
            /* light through real occluders (PTM-794). With a 64-entry stack */
            /* and builder depth <= 32 this cannot overflow; if it ever does, */
            /* report it (sentinel 2) instead of guessing. */
            if (stack_ptr > 62) {
                occluded[ray_id] = 2;
                return;
            }
            {
                if (l_hit && r_hit) {
                    /* Push farther first so nearer is popped next. */
                    if (l_tmin < r_tmin) {
                        stack[stack_ptr++] = right_idx;
                        stack[stack_ptr++] = left_idx;
                    } else {
                        stack[stack_ptr++] = left_idx;
                        stack[stack_ptr++] = right_idx;
                    }
                } else if (l_hit) {
                    stack[stack_ptr++] = left_idx;
                } else if (r_hit) {
                    stack[stack_ptr++] = right_idx;
                }
            }
        }
    }

    /* No hit - ray is visible */
    occluded[ray_id] = 0;
}

)
