// gpu_dispatch_log.h -- PTM_LOG_GPU_DISPATCH=1: one line per GPU dispatch (test/ops only).
//
// The GPU wedge mitigations size pauses and batches around one budget: how
// long a single command buffer runs (the macOS/AMD watchdog kills a queue
// somewhere above ~2 s; memory gpu-dispatch-watchdog). The thermal path
// already times its sub-dispatches (PTM_LOG_SUBDISPATCH); this switch does the
// same for the illumination (BVH) and radiosity paths, and records how long
// each inter-batch pause actually slept, so a pause can be priced and a
// ramp-down judged from measurements. Off by default; costs nothing when off.
#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

namespace gpu {

inline bool gpuDispatchLogOn() {
    static const bool on = [] {
        const char* e = std::getenv("PTM_LOG_GPU_DISPATCH");
        return e && *e && std::strcmp(e, "0") != 0;
    }();
    return on;
}

// One whole line per call (a single fputs), so lines from concurrent devices
// do not interleave mid-line.
inline void gpuDispatchLog(const std::string& line) {
    const std::string s = "[gpu-dispatch] " + line + "\n";
    std::fputs(s.c_str(), stderr);
}

}  // namespace gpu
