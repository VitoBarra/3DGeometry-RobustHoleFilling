//
// Created by vbvit on 27/05/2025.
//

#ifndef DEBUG_H
#include <string>

//Global variable and function needed for debug
fs::path meshesFolder = fs::current_path().parent_path() / "Mesh";
std::string meshName = "";

template <class MeshType>
void ExportMeshInFolder(MeshType &mesh, std::string foldername)
{
    fs::path path = meshesFolder;
    if (!foldername.empty())
        path /= foldername;
    path /= meshName;

    ExportMesh<MeshType>(mesh, path);
}

#define DEBUG_H

#endif //DEBUG_H