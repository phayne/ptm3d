// history_checkpoint.cpp -- see history_checkpoint.h.
#include "history_checkpoint.h"

#include <csignal>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

namespace ptm {
namespace history {

namespace {
constexpr char kMagic[8] = {'P', 'T', 'M', 'H', 'C', 'K', '0', '1'};
constexpr uint32_t kVersion = 1;
constexpr uint64_t kFnvOff = 1469598103934665603ULL, kFnvPrime = 1099511628211ULL;

std::atomic<bool> g_pause{false};
extern "C" void onPauseSignal(int) { g_pause.store(true); }

void putU64(std::vector<unsigned char>& b, uint64_t v) { unsigned char t[8]; std::memcpy(t, &v, 8); b.insert(b.end(), t, t + 8); }
void putStr(std::vector<unsigned char>& b, const std::string& s) { putU64(b, s.size()); b.insert(b.end(), s.begin(), s.end()); }
struct Rd {
    const unsigned char* p; size_t n, o = 0;
    bool u64(uint64_t& v) { if (o + 8 > n) return false; std::memcpy(&v, p + o, 8); o += 8; return true; }
    bool str(std::string& s) { uint64_t L; if (!u64(L) || L > n - o) return false; s.assign(reinterpret_cast<const char*>(p + o), L); o += L; return true; }
};

std::string ckptPath(const std::string& d) { return d + "/history.ckpt"; }
std::string prevPath(const std::string& d) { return d + "/history.ckpt.prev"; }

bool exists(const std::string& p) { struct stat st; return ::stat(p.c_str(), &st) == 0; }

// Parse one file. "" or an error; identity is NOT checked here.
std::string parse(const std::string& path, Checkpoint& c) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return "cannot open " + path;
    std::vector<unsigned char> b;
    unsigned char buf[1 << 16]; size_t r;
    while ((r = std::fread(buf, 1, sizeof buf, f)) > 0) b.insert(b.end(), buf, buf + r);
    std::fclose(f);
    if (b.size() < 32 || (b.size() & 7)) return path + ": truncated";
    if (std::memcmp(b.data(), kMagic, 8) != 0) return path + ": not a history checkpoint";
    uint64_t cs; std::memcpy(&cs, b.data() + b.size() - 8, 8);
    if (digest(b.data(), b.size() - 8) != cs) return path + ": checksum mismatch (corrupt or truncated)";
    Rd R{b.data(), b.size() - 8, 8};
    uint64_t ver, nid, nb;
    if (!R.u64(ver) || ver != kVersion) return path + ": unsupported checkpoint version";
    if (!R.u64(c.next_epoch) || !R.u64(nid)) return path + ": truncated";
    c.id.fields.clear();
    for (uint64_t i = 0; i < nid; ++i) { std::string k, v; if (!R.str(k) || !R.str(v)) return path + ": truncated identity"; c.id.fields[k] = v; }
    if (!R.u64(nb)) return path + ": truncated";
    c.blobs.assign(nb, Blob());
    for (uint64_t i = 0; i < nb; ++i) {
        std::string data;
        if (!R.str(c.blobs[i].name) || !R.str(data)) return path + ": truncated blob";
        c.blobs[i].bytes.assign(data.begin(), data.end());
    }
    return "";
}
}  // namespace

uint64_t digest(const void* p, size_t n) {
    const unsigned char* c = static_cast<const unsigned char*>(p);
    uint64_t h = kFnvOff; size_t i = 0;
    for (; i + 8 <= n; i += 8) { uint64_t w; std::memcpy(&w, c + i, 8); h ^= w; h *= kFnvPrime; }
    if (i < n) { uint64_t w = 0; std::memcpy(&w, c + i, n - i); h ^= w; h *= kFnvPrime; }
    return h;
}
std::string hex(uint64_t v) { char b[17]; std::snprintf(b, sizeof b, "%016llx", (unsigned long long)v); return b; }

std::string Identity::diff(const Identity& o) const {
    std::string d;
    for (const auto& kv : fields) {
        auto it = o.fields.find(kv.first);
        if (it == o.fields.end()) d += "  " + kv.first + ": '" + kv.second + "' vs (absent)\n";
        else if (it->second != kv.second) d += "  " + kv.first + ": '" + kv.second + "' vs '" + it->second + "'\n";
    }
    for (const auto& kv : o.fields) if (!fields.count(kv.first)) d += "  " + kv.first + ": (absent) vs '" + kv.second + "'\n";
    return d;
}

const Blob* Checkpoint::find(const std::string& name) const {
    for (const auto& b : blobs) if (b.name == name) return &b;
    return nullptr;
}

std::string writeCheckpoint(const std::string& dir, const Checkpoint& c) {
    std::vector<unsigned char> b(kMagic, kMagic + 8);
    putU64(b, kVersion); putU64(b, c.next_epoch); putU64(b, c.id.fields.size());
    for (const auto& kv : c.id.fields) { putStr(b, kv.first); putStr(b, kv.second); }
    putU64(b, c.blobs.size());
    for (const auto& bl : c.blobs) { putStr(b, bl.name); putU64(b, bl.bytes.size()); b.insert(b.end(), bl.bytes.begin(), bl.bytes.end()); }
    while (b.size() & 7) b.push_back(0);
    putU64(b, digest(b.data(), b.size()));
    const std::string tmp = ckptPath(dir) + ".tmp";
    const int fd = ::open(tmp.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) return "cannot create " + tmp;
    size_t off = 0;
    while (off < b.size()) {
        const ssize_t w = ::write(fd, b.data() + off, b.size() - off);
        if (w <= 0) { ::close(fd); return "write failed: " + tmp; }
        off += (size_t)w;
    }
    if (::fsync(fd) != 0) { ::close(fd); return "fsync failed: " + tmp; }
    ::close(fd);
    // Rotate: a crash between the two renames leaves .prev (valid) + .tmp, and load() finds .prev.
    if (exists(ckptPath(dir)) && std::rename(ckptPath(dir).c_str(), prevPath(dir).c_str()) != 0)
        return "cannot rotate " + ckptPath(dir);
    if (std::rename(tmp.c_str(), ckptPath(dir).c_str()) != 0) return "cannot rename " + tmp;
    const int dfd = ::open(dir.c_str(), O_RDONLY);   // persist the directory entries
    if (dfd >= 0) { ::fsync(dfd); ::close(dfd); }
    return "";
}

std::string loadCheckpoint(const std::string& dir, const Identity& expect, Checkpoint& out, bool& found,
                           std::string& note) {
    found = false; note.clear();
    const bool have_cur = exists(ckptPath(dir)), have_prev = exists(prevPath(dir));
    if (!have_cur && !have_prev) return "";
    std::string e1 = "absent";
    if (have_cur) {
        e1 = parse(ckptPath(dir), out);
        if (e1.empty()) {
            const std::string d = expect.diff(out.id);
            if (!d.empty()) return ckptPath(dir) + " belongs to a DIFFERENT job (identity differs: run vs file):\n" + d;
            found = true; return "";
        }
    }
    if (have_prev) {
        const std::string e2 = parse(prevPath(dir), out);
        if (e2.empty()) {
            const std::string d = expect.diff(out.id);
            if (!d.empty()) return prevPath(dir) + " belongs to a DIFFERENT job (identity differs: run vs file):\n" + d;
            note = "WARNING: latest checkpoint unusable (" + e1 + "); resuming from the PREVIOUS one, epoch " +
                   std::to_string(out.next_epoch) + " -- the history redoes the epochs since then.";
            found = true; return "";
        }
        return "no usable checkpoint: latest: " + e1 + "; previous: " + e2;
    }
    return "no usable checkpoint: latest: " + e1 + "; no previous checkpoint";
}

PauseControl::PauseControl(std::string dir) : dir_(std::move(dir)) {
    const int sigs[3] = {SIGUSR1, SIGTERM, SIGINT};
    for (int i = 0; i < 3; ++i) {
        struct sigaction sa{}, old{};
        sa.sa_handler = onPauseSignal; sigemptyset(&sa.sa_mask); sa.sa_flags = 0;
        sigaction(sigs[i], &sa, &old);
        auto* keep = new struct sigaction(old); old_[i] = keep;
    }
}
PauseControl::~PauseControl() {
    const int sigs[3] = {SIGUSR1, SIGTERM, SIGINT};
    for (int i = 0; i < 3; ++i) {
        auto* keep = static_cast<struct sigaction*>(old_[i]);
        if (keep) { sigaction(sigs[i], keep, nullptr); delete keep; }
    }
}
bool PauseControl::pauseFilePresent() const { return exists(dir_ + "/PAUSE"); }
bool PauseControl::requested() const { return g_pause.load() || pauseFilePresent(); }
void PauseControl::requestForTest() { g_pause.store(true); }
void PauseControl::clearForTest() { g_pause.store(false); }

}  // namespace history
}  // namespace ptm
