#include <direct.h>
#include <iostream>

#include <vcg/complex/complex.h>
#include <wrap/io_trimesh/import_ply.h>
#include <vcg/complex/algorithms/clean.h>
#include <wrap/io_trimesh/export_ply.h>
#include <vcg/complex/algorithms/update/topology.h>
#include "file.h"
#include "vcg/complex/algorithms/hole.h"
#include <vcg/complex/algorithms/isotropic_remeshing.h>

using namespace  vcg;
class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public UsedTypes<Use<MyVertex>   ::AsVertexType,
                                           Use<MyEdge>     ::AsEdgeType,
                                           Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::VFAdj,vertex::BitFlags>{};
class MyFace    : public Face<   MyUsedTypes, face::FFAdj, face::VFAdj ,face::VertexRef,face::Normal3f ,face::BitFlags> {};
class MyEdge    : public Edge<   MyUsedTypes> {};

class MyMesh    : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};

bool DO_CLEAN_UP = false;

void execute_mesh_refinement(MyMesh& mesh, float targetEdgeLen)
{
  tri::IsotropicRemeshing<MyMesh>::Params params;
  params.SetFeatureAngleDeg(181.0f);
  params.adapt        = false;
  params.selectedOnly = true;
  params.splitFlag    = true;
  params.collapseFlag = true;
  params.swapFlag     = true;
  params.smoothFlag   = true;
  params.projectFlag  = false;
  params.surfDistCheck= false;

  // Refinement and smoothing can be tricky. Usually it is good to
  // 1) start with large tris to get fast convergence to the min surf
  // 2) switch a bit to small tri to unfold bad things at the boundary
  // 3) go for the desired edge len
  // Rinse and repeat.
  for (int k = 0; k < 3; k++)
  {
    params.SetTargetLen(targetEdgeLen * 3.0);
    params.iter = 5;
    tri::IsotropicRemeshing<MyMesh>::Do(mesh, params);

    params.SetTargetLen(targetEdgeLen / 3.0);
    params.iter = 3;
    tri::IsotropicRemeshing<MyMesh>::Do(mesh, params);

    params.SetTargetLen(targetEdgeLen);
    params.iter = 2;
    tri::IsotropicRemeshing<MyMesh>::Do(mesh, params);
  }
}

int main( int argc, char **argv )
{
  fs::path meshesFolder = fs::current_path().parent_path() / "TestMesh";
  for (auto &v : GetFilesInFolderExtension(meshesFolder ,".ply"))
  {
    MyMesh mesh;
    std::string meshName =fs::path(v).filename().string();

    if(tri::io::ImporterPLY<MyMesh>::Open(mesh,v.c_str())!=0)
    {
      printf("Error reading file  %s\n",meshName);
      continue;
    }
    size_t originalFaceNumber= mesh.FN();
    printf("mesh %s has vertexes:%i faces:%i\n",meshName.c_str(),mesh.VN(),originalFaceNumber);
    // Counting the number of edges using FF adjacency

    //Clean up unreferenced vertices if cleanup flag is set
    if (DO_CLEAN_UP)
      tri::Clean<MyMesh>::RemoveUnreferencedVertex(mesh);
    tri::UpdateTopology<MyMesh>::FaceFace(mesh); // Compute information for face-to-face adjacency
    tri::UpdateTopology<MyMesh>::VertexFace(mesh); // Compute information for vertex-to-face adjacency
    tri::UpdateFlags<MyMesh>::FaceBorderFromNone(mesh);
    
    // Fill holes in the mesh using ear cutting algorithm with minimum weight criterion
    tri::Hole<MyMesh>::EarCuttingFill<vcg::tri::MinimumWeightEar< MyMesh> >(mesh,500,false,nullptr);
    assert(tri::Clean<MyMesh>::IsFFAdjacencyConsistent(mesh));
    tri::UpdateNormal<MyMesh>::NormalizePerFaceByArea(mesh);

    tri::UpdateSelection<MyMesh>::FaceClear(mesh);
    size_t faceIndex = originalFaceNumber;
    for (; faceIndex < mesh.FN(); ++faceIndex)
      if (!mesh.face[faceIndex].IsD()) mesh.face[faceIndex].SetS();
    std::cout << "number of face selected " <<   faceIndex - originalFaceNumber<<std::endl;

    execute_mesh_refinement(mesh,  0.1555f);

    // Apply cotangent smoothing
    //vcg::tri::Smooth<MyMesh>::VertexCoordTaubin(mesh, 20, 0.5,0.5,false);


    int mask = 0; // Use this to specify what attributes to save (colors, normals, etc.)
    fs::path outputFilename = meshesFolder/"exportWithHoleRefinement"/meshName;
    mkdir(outputFilename.parent_path().string().c_str());
    int result = tri::io::ExporterPLY<MyMesh>::Save(mesh, outputFilename.string().c_str(), true);
    if (result != 0)
      std::cerr << "Error saving mesh: " << tri::io::ExporterPLY<MyMesh>::ErrorMsg(result) << std::endl;
    else
      std::cout << "Mesh successfully saved to " << outputFilename << std::endl;

  }
  return 0;
}
