// history_checkpoint.h -- pause / checkpoint / resume for the streamed history marcher (M1).
//
// Paul's requirement (2026-11-09): the per-site CPU history job runs on the MacBook alongside other
// work, so it must be cleanly pausable and resumable, and a resumed job must produce the SAME BYTES
// as an uninterrupted one. This module is the mechanism; test_history_checkpoint is the gate.
//
// What a checkpoint must hold (and why nothing else):
//   * the column temperatures T of every face (ThermalBatchSoA::advanceTimestep rebuilds k, cp and
//     kappa from T at the start of every step, and kc / rho are static per face from the grid and
//     the config, so the other per-batch arrays are DERIVED, not state -- the gate proves it by
//     restoring T alone);
//   * every vector the marcher carries between epochs (e.g. the pipelined first visible bounce);
//   * the next epoch index and any running counters (classifier census, audit);
//   * the job identity (below), checked field by field on resume.
//
// Files in the job directory: history.ckpt (latest), history.ckpt.prev (the one before), written
// atomically (temp file, fsync, rotate current -> .prev, rename temp -> current). A corrupt latest
// falls back to .prev LOUDLY; an identity mismatch never falls back (it is a different job).
//
// Pause: create the file PAUSE in the job directory, or send SIGUSR1 / SIGTERM / SIGINT. The marcher
// finishes the epoch in progress, writes a checkpoint, and exits with kExitPaused. Resume refuses to
// start while PAUSE exists (remove it to resume), so a pause cannot be undone by accident.
//
// Determinism the byte-identical gate relies on: the canonical CSR kernel reduces each row
// independently in a pinned order, and the thermal column of each face is independent of every
// other face, so neither result depends on the thread count -- a job may resume with a different
// thread cap (MacBook 18 -> 12 when Paul needs the machine) and still match.
#pragma once

#include <atomic>
#include <cstring>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace ptm {
namespace history {

constexpr int kExitPaused = 75;   // EX_TEMPFAIL: "paused, resume later"

// Ordered name -> value pairs naming every input the history depends on (CSR md5, .hzc key and
// checksum, config hash, epoch-grid hash, code version, ...). Compared field by field.
struct Identity {
    std::map<std::string, std::string> fields;
    void set(const std::string& k, const std::string& v) { fields[k] = v; }
    // "" if equal, else a human-readable list of the differing fields.
    std::string diff(const Identity& other) const;
};

struct Blob {
    std::string name;
    std::vector<unsigned char> bytes;
    template <class T> static Blob of(const std::string& n, const std::vector<T>& v) {
        Blob b; b.name = n; b.bytes.resize(v.size() * sizeof(T));
        if (!v.empty()) std::memcpy(b.bytes.data(), v.data(), b.bytes.size());
        return b;
    }
    template <class T> bool to(std::vector<T>& v) const {
        if (bytes.size() % sizeof(T)) return false;
        v.resize(bytes.size() / sizeof(T));
        if (!v.empty()) std::memcpy(v.data(), bytes.data(), bytes.size());
        return true;
    }
};

struct Checkpoint {
    Identity id;
    uint64_t next_epoch = 0;
    std::vector<Blob> blobs;
    const Blob* find(const std::string& name) const;
};

// Atomic write into dir (see above). "" or an error message.
std::string writeCheckpoint(const std::string& dir, const Checkpoint& c);

// Load the latest valid checkpoint whose identity equals `expect`. Returns "" and fills `out`,
// or an error. `found` is false (and the return "") when the directory holds no checkpoint at all
// (a fresh start). `note` receives any warning (e.g. fell back to .prev) -- the caller must print it.
std::string loadCheckpoint(const std::string& dir, const Identity& expect, Checkpoint& out, bool& found,
                           std::string& note);

// Pause requests: signals (SIGUSR1, SIGTERM, SIGINT) set a flag; the PAUSE file is polled.
class PauseControl {
public:
    explicit PauseControl(std::string dir);   // installs the signal handlers
    ~PauseControl();                           // restores the previous handlers
    bool requested() const;                    // signal received or dir/PAUSE exists
    bool pauseFilePresent() const;
    static void requestForTest();              // what a signal does (for the gate)
    static void clearForTest();
private:
    std::string dir_;
    void* old_[3] = {nullptr, nullptr, nullptr};
};

// Checkpoint cadence: every `every` epochs (e.g. N lunations x 120), and always on pause / at the end.
struct CheckpointPolicy {
    uint64_t every = 120;
    bool due(uint64_t epochs_done) const { return every && epochs_done && (epochs_done % every == 0); }
};

// 64-bit FNV-1a over little-endian u64 words (tail zero-padded): content digests for Identity.
uint64_t digest(const void* p, size_t n);
std::string hex(uint64_t v);

}  // namespace history
}  // namespace ptm
