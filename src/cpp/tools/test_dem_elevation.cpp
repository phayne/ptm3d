// Checks src/cpp/tools/dem_elevation.h against scripts/dem_elevation_cases.tsv,
// the same table src/python/tests/test_dem_elevation_rule.py runs through
// scripts/dem_elevation.py. No GDAL needed.
//
//   ./build/bin/test_dem_elevation [path/to/dem_elevation_cases.tsv]
//
// Exit 0 when every case matches, 1 otherwise.
#include "dem_elevation.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#ifndef PTM_DEM_ELEVATION_CASES
#define PTM_DEM_ELEVATION_CASES "scripts/dem_elevation_cases.tsv"
#endif

namespace {
std::optional<double> opt(const std::string& s) {
    if (s == "-") return std::nullopt;
    return std::stod(s);
}
bool close(double a, double b) { return std::abs(a - b) <= 1e-9 * std::max(1.0, std::abs(b)); }
}  // namespace

int main(int argc, char** argv) {
    const std::string path = argc > 1 ? argv[1] : PTM_DEM_ELEVATION_CASES;
    std::ifstream in(path);
    if (!in) {
        std::cerr << "test_dem_elevation: cannot open " << path << "\n";
        return 1;
    }
    int n = 0, fails = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> c;
        std::stringstream ss(line);
        std::string f;
        while (std::getline(ss, f, '\t')) c.push_back(f);
        if (c.size() != 11) {
            std::cerr << "FAIL: malformed case line: " << line << "\n";
            ++fails;
            continue;
        }
        ++n;
        const std::string unit = c[3] == "-" ? "" : c[3];
        std::string got;
        bool ok = false;
        try {
            auto t = ptm::dem::resolve(opt(c[1]), opt(c[2]), unit, c[4] == "1", opt(c[5]),
                                       opt(c[6]), std::stod(c[7]));
            got = "OK scale " + ptm::dem::fmt(t.scale) + " offset " + ptm::dem::fmt(t.offset);
            ok = c[8] == "OK" && close(t.scale, std::stod(c[9])) && close(t.offset, std::stod(c[10]));
        } catch (const ptm::dem::DemElevationError& e) {
            got = std::string("ERROR ") + e.what();
            ok = c[8] == "ERROR";
        }
        std::cout << (ok ? "ok:   " : "FAIL: ") << c[0] << ": " << got << "\n";
        if (!ok) ++fails;
    }
    if (n == 0) {
        std::cerr << "FAIL: no cases read from " << path << "\n";
        return 1;
    }
    std::cout << (fails ? std::to_string(fails) + " FAILURE(S)" : "PASS") << " (" << n << " cases)\n";
    return fails ? 1 : 0;
}
