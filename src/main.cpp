#include <chrono>
#include <iostream>

#include <cmath>
#include <vector>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/isotropic_remeshing.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/simplex/face/topology.h>
#include <wrap/io_trimesh/import_ply.h>
#include "vcg/complex/algorithms/hole.h"

#include "HoleFilling.h"
#include "file.h"
#include "MeshIO.h"


using namespace  vcg;
class MyVertex; class MyEdge; class MyFace;




// Base types definition for the mesh components
struct MyUsedTypes : public UsedTypes<Use<MyVertex>::AsVertexType,
                                      Use<MyEdge>::AsEdgeType,
                                      Use<MyFace>::AsFaceType>
{
};

// Vertex class with 3D coordinates, normals, vertex-face adjacency and flags
class MyVertex : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::VFAdj,vertex::VEAdj, vertex::BitFlags, vertex::Mark, vertex::Qualityf>
{
};

// Face class with face-face adjacency, vertex-face adjacency, vertex references, normals and flags
class MyFace : public Face<MyUsedTypes, face::FFAdj, face::VFAdj, face::VertexRef, face::Normal3f, face::BitFlags>
{
};

// Basic edge class
class MyEdge : public Edge<MyUsedTypes>
{
};

// Main mesh class composed of vertices and faces stored in vectors
class MyMesh    : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};

//Global variable and function needed for debug
fs::path g_MeshesFolder = fs::current_path().parent_path() / "Mesh";
std::string g_MeshName;

template <class MeshType>
void ExportMeshInFolder(MeshType &mesh, const std::string& folderName)
{
    ExportMeshInFolder(mesh, g_MeshesFolder,folderName,g_MeshName);
}



std::string getWeightTypeSuffix(TypeOfWeight tow)
{
    switch (tow)
    {
    case uniform: return "_uniform";
    case distance: return "_distance";
    case cotangent: return "_cotangent";
    default: return "_uniform";
    }
}

int main( int argc, char **argv )
{
    for (auto &v : GetFilesInFolderExtension(g_MeshesFolder ,".ply"))
    {
        MyMesh mesh;
        g_MeshName =fs::path(v).filename().string();

        if(tri::io::ImporterPLY<MyMesh>::Open(mesh,v.c_str())!=0)
        {
            std::cout << "Error reading file " << g_MeshName << std::endl;
            continue;
        }
        size_t originalFaceNumber = mesh.FN();
        std::cout << "mesh " << g_MeshName << " has vertexes:" << mesh.VN() << " faces:" << originalFaceNumber <<
            std::endl;
        // Counting the number of edges using FF adjacency

        //Clean up unreferenced vertices
        tri::Clean<MyMesh>::RemoveUnreferencedVertex(mesh);
        tri::UpdateFlags<MyMesh>::FaceBorderFromNone(mesh); // Initialize border flags
        tri::UpdateTopology<MyMesh>::FaceFace(mesh); // Compute information for face-to-face adjacency
        tri::UpdateTopology<MyMesh>::VertexFace(mesh); // Compute information for vertex-to-face adjacency

        // Fill holes in the mesh using ear cutting algorithm with minimum weight criterion
        tri::Hole<MyMesh>::EarCuttingFill<tri::MinimumWeightEar< MyMesh> >(mesh,500,false,nullptr);
        assert(tri::Clean<MyMesh>::IsFFAdjacencyConsistent(mesh));

        // update mesh topology information
        tri::UpdateTopology<MyMesh>::FaceFace(mesh);
        tri::UpdateTopology<MyMesh>::VertexFace(mesh);

        // Clear selection
        tri::UpdateSelection<MyMesh>::FaceClear(mesh);
        tri::UpdateSelection<MyMesh>::VertexClear(mesh);

        ExportMeshInFolder(mesh, "HoleFilled");


        size_t vertexSelected = 0;
        size_t faceIndex = originalFaceNumber;
        // First, select all new faces
        for (; faceIndex < mesh.FN(); ++faceIndex)
        {
            if (mesh.face[faceIndex].IsD())
                continue;
            mesh.face[faceIndex].SetS();
            // Select all vertices of this face
            for (int v = 0; v < mesh.face[faceIndex].VN(); ++v)
            {
                mesh.face[faceIndex].V(v)->SetS();
                vertexSelected++;
            }
        }

        std::cout << "number of face selected " << faceIndex - originalFaceNumber << std::endl;
        std::cout << "number of vertex selected " << vertexSelected << std::endl;
        ExportMeshInFolder(mesh, "HoleFilledAndSelected");

        HolePatchRefinement<MyMesh>(mesh);
        // Compact mesh to remove deleted elements
        tri::Allocator<MyMesh>::CompactFaceVector(mesh);
        tri::Allocator<MyMesh>::CompactVertexVector(mesh);
        tri::UpdateSelection<MyMesh>::VertexClear(mesh);
        tri::UpdateTopology<MyMesh>::FaceFace(mesh);
        tri::UpdateTopology<MyMesh>::VertexFace(mesh);

        // Then select vertices that have all adjacent faces selected
        faceIndex = originalFaceNumber;
        for (; faceIndex < mesh.FN(); ++faceIndex)
        {
            if (mesh.face[faceIndex].IsD())
                continue;

            for (int v=0 ; v<mesh.face[faceIndex].VN(); ++v)
            {
                auto vertex = mesh.face[faceIndex].V(v);
                if (vertex->IsD()) continue;


                vcg::face::VFIterator<MyFace> vfi(vertex);

                int faceCount = 0;
                bool allFacesSelected = true;
                for (; !vfi.End(); ++vfi)
                {
                    if (vfi.F()->IsD() || vfi.F() == nullptr)
                        continue;

                    if (!vfi.F()->IsS())
                    {
                        faceCount++;
                        allFacesSelected = false;
                    }
                }

                if (allFacesSelected)
                {
                    vertex->SetS();
                    vertexSelected++;
                }
            }
        }

        std::cout << "number of vertex selected " << vertexSelected << std::endl;

        ExportMeshInFolder(mesh, "HoleFilledAndRefined");


        // Apply cotangent smoothing
        for (int i =0 ; i<3; i++)
        {
            MyMesh meshCopy;
            TypeOfWeight tow = static_cast<TypeOfWeight>(i);
            std::string weightTypeSuffix = getWeightTypeSuffix(tow);

            CloneMesh(mesh, meshCopy);

            auto start = std::chrono::high_resolution_clock::now();
            HoleFairing<MyMesh>(meshCopy, 100,1, true, tow);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end - start;

            std::cout << "HoleFairing Duration (" << weightTypeSuffix << "): " << duration.count() << " seconds" <<
                std::endl;
            ExportMeshInFolder(meshCopy,"HoleFilledAndRefined"+ weightTypeSuffix);
        }
    }
    return 0;
}

