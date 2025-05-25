#ifndef FILE_H
#define FILE_H
#include <filesystem>
#include <vector>
#include <string>
#include <sys/stat.h>

namespace fs = std::filesystem;

// Function to get all .off file paths from a folder
inline std::vector<std::string> GetFilesInFolderExtension(const fs::path& folder_path,const std::string& extension = ".ply" ) {
    std::vector<std::string> off_files;
    //working directory + path
    try {
        if (!fs::exists(folder_path)) {
            std::cerr << "Directory does not exist: " << folder_path << std::endl;
            return off_files;
        }

        for (const auto& entry : fs::directory_iterator(folder_path)) {
            if (entry.is_regular_file() && entry.path().extension() == extension) {
                off_files.push_back(entry.path().string());
            }
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
    }

    return off_files;
}



template<class MeshType>
void ExportMesh(const MeshType &mesh,fs::path path)
{
    std::cout << "exporting:"<< path.filename() << std::endl;

    fs::create_directories(path.parent_path());
    // Save mesh with flags (mask) enabled
    int result = vcg::tri::io::ExporterPLY<MeshType>::Save(mesh, path.string().c_str(), vcg::tri::io::Mask::IOM_ALL);
    if (result != 0)
        std::cerr << "Error saving mesh: " << vcg::tri::io::ExporterPLY<MeshType>::ErrorMsg(result) << std::endl;
    else
        std::cout << "Mesh successfully saved to " << path << std::endl;
}
#endif
