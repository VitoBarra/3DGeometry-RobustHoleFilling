#ifndef FILE_H
#define FILE_H
#include <filesystem>
#include <vector>
#include <string>

namespace fs = std::filesystem;

// Function to get all .off file paths from a folder
inline std::vector<std::string> GetFilesInFolderExtension(const std::string& folder_path,const std::string& extension = ".ply" ) {
    std::vector<std::string> off_files;
    for (const auto& entry : fs::directory_iterator(folder_path)) {
        if (entry.is_regular_file() && entry.path().extension() == extension) {
            off_files.push_back(entry.path().string());
        }
    }
    return off_files;
}
#endif
