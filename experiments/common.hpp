// Shared helpers for the experiment programs.
#pragma once

#include <cstdio>
#include <filesystem>
#include <string>

inline std::FILE* open_csv(const std::string& name, const char* header) {
    std::filesystem::create_directories("docs/data");
    const std::string path = "docs/data/" + name;
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) {
        std::perror(path.c_str());
        std::exit(1);
    }
    std::fprintf(f, "%s\n", header);
    std::printf("writing %s\n", path.c_str());
    return f;
}
