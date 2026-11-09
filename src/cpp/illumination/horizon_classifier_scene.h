// horizon_classifier_scene.h -- the Embree scene and ray predicate the horizon classifier
// certifies against (history marcher M0). Port of the measurement tool's scene
// (the 2026-11-08 measurement tool hzcls.cpp and the horizon-cache probe hcommon.h; not in the repo),
// which replicates the engine's chord path under the production environment
// (PTM_ILLUM_TRUE_NORMALS=1 PTM_ILLUM_RISE_EPS=1 PTM_ILLUM_RISE_HABS=0):
//   * core [+ ring + seam curtain], float32 vertices recentered on the core+ring bbox center;
//     ROBUST | COMPACT | FILTER_FUNCTION_IN_ARGUMENTS, build quality HIGH;
//   * ray from the face CENTROID, tnear 0, any-hit filter rejecting only the source core
//     primitive (h_eps = 0);
//   * the R_occ sphere pre-check for open cores (in the classifier library).
// ⛔ M1 must re-verify this against the engine's chord predicate on the SAME build before the
// marcher uses it: the key's predicate descriptor names this predicate, so an engine-side
// change of predicate requires a new descriptor (and therefore new .hzc files).
#pragma once

#include "horizon_classifier.h"

#include <embree4/rtcore.h>

#include <array>
#include <string>
#include <vector>

namespace ptm {
namespace hzc {

using F3 = std::array<int, 3>;

class ClassifierScene {
public:
    // ring may be empty. orientation_sign has one entry per core face.
    ClassifierScene(const std::vector<V3>& Vc, const std::vector<F3>& Fc,
                    const std::vector<V3>& Vr, const std::vector<F3>& Fr,
                    std::vector<int8_t> orientation_sign, std::string orientation_rule);
    ~ClassifierScene();
    ClassifierScene(const ClassifierScene&) = delete;
    ClassifierScene& operator=(const ClassifierScene&) = delete;

    size_t numCoreFaces() const { return Fc_.size(); }
    size_t numRingFaces() const { return nring_; }
    size_t numCurtainFaces() const { return ncurtain_; }
    size_t numBoundaryEdges() const { return nbound_; }
    const SceneFacts& facts() const { return facts_; }
    FaceFrame frame(size_t core_face) const;
    // 0 clear, 1 blocked: Embree occlusion from the face centroid along d.
    int occluded(size_t core_face, const FaceFrame& f, const V3& d) const;

    uint64_t sceneHash() const { return scene_hash_; }
    static const char* predicateDescriptor();
    Key key(const Params& p) const;
    const std::string& orientationRule() const { return orientation_rule_; }

private:
    std::vector<V3> Vc_; std::vector<F3> Fc_;
    std::vector<int8_t> sign_; std::string orientation_rule_;
    size_t nring_ = 0, ncurtain_ = 0, nbound_ = 0;
    SceneFacts facts_;
    V3 off_{};
    uint64_t scene_hash_ = 0;
    RTCDevice dev_ = nullptr;
    RTCScene sc_ = nullptr;
};

}  // namespace hzc
}  // namespace ptm
