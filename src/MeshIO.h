//
// Created by vbvit on 27/05/2025.
//

#ifndef MESHOP_H
#define MESHOP_H

#include <filesystem>
#include <string>

using namespace  vcg;
namespace fs = std::filesystem;

template<class MeshType>
void ExportMesh(const MeshType &mesh, const fs::path& path)
{
    std::cout << "exporting:"<< path.filename() << std::endl;

    fs::create_directories(path.parent_path());
    // Save mesh with flags (mask) enabled
    int result = tri::io::ExporterPLY<MeshType>::Save(mesh, path.string().c_str(), vcg::tri::io::Mask::IOM_ALL);
    if (result != 0)
        std::cerr << "Error saving mesh: " << tri::io::ExporterPLY<MeshType>::ErrorMsg(result) << std::endl;
    else
        std::cout << "Mesh successfully saved to " << path << std::endl;
}

template <class MeshType>
void ExportMeshInFolder(MeshType &mesh, const fs::path& folder ,const std::string& subFolderName , const std::string& meshName)
{
    fs::path path = folder;
    if (!subFolderName.empty())
        path /= subFolderName;
    path /= meshName;

    ExportMesh<MeshType>(mesh, path);
}

template <class MeshType>
static void CloneMesh(MeshType& source, MeshType& target)
{
    tri::Append<MeshType, MeshType>::MeshCopy(target, source);
}


#endif //MESHOP_H