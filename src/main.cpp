#include <direct.h>
#include <iostream>
#include <chrono>

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


// LaplacianInfo: Helper class for Laplacian smoothing operations
// Stores the sum of neighboring vertices' positions and the count of neighbors
// Used for calculating the average position in smoothing operations
template <typename MeshType>
class LaplacianInfo
{
public:
    LaplacianInfo(const MeshType::CoordType& _p, const int _n) : sum(_p), cnt(_n)
    {
    }

    LaplacianInfo()
    {
    }

    typename MeshType::CoordType sum; // Sum of neighboring vertices' positions
    typename MeshType::ScalarType cnt; // Count of neighboring vertices
};

// Base types definition for the mesh components
struct MyUsedTypes : public UsedTypes<Use<MyVertex>::AsVertexType,
                                      Use<MyEdge>::AsEdgeType,
                                      Use<MyFace>::AsFaceType>
{
};

// Vertex class with 3D coordinates, normals, vertex-face adjacency and flags
class MyVertex : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::VFAdj, vertex::BitFlags, vertex::Mark>
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

bool DO_CLEAN_UP = false;






template <class MeshType>
void execute_mesh_refinement(MeshType& mesh, float targetEdgeLen)
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



// Weight types for mesh operations:
// uniform - all edges have equal weight
// distance - weight based on edge length
// cotangent - weight based on angles (for preserving geometric features)
enum TypeOfWeight
{
    uniform = 0, // Equal weights
    distance = 1, // Distance-based weights
    cotangent = 2, // Cotangent weights (angle-based)
};


// calculate the cotangent function, it is stable around 0 and pi but can have some instabilities around pi/2
float cot(float angle)
{
    return tan(M_PI *0.5 - angle);
}


template <class MeshType>
static float calculateWeight(TypeOfWeight tow,
                             const typename MeshType::FaceType* f,
                             const int edge)
{
    switch (tow)
    {
    default:
    case uniform:
        return 1.0f;
    case distance:
        return Distance(f->V0(edge)->P(), f->V1(edge)->P());
    case cotangent:
        {

            float cotWeight = 0.0f;
            // Calculate cotangent for current face
            float alpha = Angle(f->P1(edge) - f->P2(edge), f->P0(edge) - f->P2(edge));
            cotWeight += cot(alpha);

            typename MeshType::FaceType* ff = f->FFp(edge);
            // Add cotangent from the adjacent face if it exists
            if (ff != nullptr && ff != f) // if f it is connected to itself, then the FF face does not exist (it is a border)
            {
                int ffIndex = f->FFi(edge);
                float beta = Angle(ff->P1(ffIndex) - ff->P2(ffIndex),
                                   ff->P0(ffIndex) - ff->P2(ffIndex));
                cotWeight += cot(beta);
            }
            return cotWeight;
        }
    }
}

template <class MeshType>
static void AccumulateLaplacianInfo(MeshType &m, SimpleTempData<typename MeshType::VertContainer, LaplacianInfo<MeshType>> &TD, TypeOfWeight weightType = uniform)
{
    if (weightType == cotangent)
        tri::RequireFFAdjacency(m);

    float weight = 1.0f;

    typename MeshType::FaceIterator fi;
    for (fi = m.face.begin(); fi != m.face.end(); ++fi)
    {
        if ((*fi).IsD())
            continue;
        for (int edge = 0; edge < 3; ++edge)
        {
            if ((*fi).IsB(edge))
                continue;

            weight = calculateWeight<MyMesh>(weightType, &(*fi), edge);

            TD[(*fi).V0(edge)].sum += (*fi).P1(edge) * weight;
            TD[(*fi).V1(edge)].sum += (*fi).P0(edge) * weight;
            TD[(*fi).V0(edge)].cnt += weight;
            TD[(*fi).V1(edge)].cnt += weight;
        }
    }
    // si azzaera i dati per i vertici di bordo
    for (fi = m.face.begin(); fi != m.face.end(); ++fi)
    {
        if ((*fi).IsD())
            continue;
        for (int j = 0; j < 3; ++j)
            if ((*fi).IsB(j))
            {
                TD[(*fi).V0(j)].sum = (*fi).P0(j);
                TD[(*fi).V1(j)].sum = (*fi).P1(j);
                TD[(*fi).V0(j)].cnt = 1;
                TD[(*fi).V1(j)].cnt = 1;
            }
    }

    // se l'edge j e' di bordo si deve mediare solo con gli adiacenti
    for (fi = m.face.begin(); fi != m.face.end(); ++fi)
    {
        if (!(*fi).IsD())
            for (int j = 0; j < 3; ++j)
                if ((*fi).IsB(j))
                {
                    TD[(*fi).V(j)].sum += (*fi).V1(j)->P();
                    TD[(*fi).V1(j)].sum += (*fi).V(j)->P();
                    ++TD[(*fi).V(j)].cnt;
                    ++TD[(*fi).V1(j)].cnt;
                }
    }
}

// Parameters:
//   m - Input mesh to be smoothed
//   step - Number of smoothing iterations
//   alpha - Smoothing factor (0-1), controls the intensity of smoothing
//   SmoothSelected - If true, only smooth selected vertices
template <class MeshType>
static void CloneMesh(MeshType& source, MeshType& target)
{
    tri::Append<MeshType, MeshType>::MeshCopy(target, source);
}

template <class MeshType>
static void HoleFairing(MeshType &m, int step, float alpha, bool SmoothSelected = false , TypeOfWeight weightType = uniform)
{
    typename MeshType::VertexIterator vi;
    LaplacianInfo<MyMesh> lpz(typename MeshType::CoordType(0, 0, 0), 0);
    assert(alpha <= 1.0f && alpha > 0.f);
    SimpleTempData<typename MeshType::VertContainer, LaplacianInfo<MeshType>> TD(m.vert);

    for (int i = 0; i < step; ++i)
    {
        TD.Init(lpz);
        AccumulateLaplacianInfo(m, TD,weightType);
        for (vi = m.vert.begin(); vi != m.vert.end(); ++vi)
            if (!(*vi).IsD() && TD[*vi].cnt > 0)
            {
                if (!SmoothSelected || (*vi).IsS())
                {
                    typename MeshType::CoordType Delta = TD[*vi].sum / TD[*vi].cnt - (*vi).P();
                    (*vi).P() = (*vi).P() + Delta * alpha;
                }
            }
    }
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
    fs::path meshesFolder = fs::current_path().parent_path() / "TestMesh";
    fs::path outputFilename ="";
    for (auto &v : GetFilesInFolderExtension(meshesFolder ,".ply"))
    {
        MyMesh Mesh;
        std::string meshName =fs::path(v).filename().string();

        if(tri::io::ImporterPLY<MyMesh>::Open(Mesh,v.c_str())!=0)
        {
            printf("Error reading file  %s\n",meshName);
            continue;
        }
        size_t originalFaceNumber= Mesh.FN();
        printf("mesh %s has vertexes:%i faces:%i\n",meshName.c_str(),Mesh.VN(),originalFaceNumber);
        // Counting the number of edges using FF adjacency

        //Clean up unreferenced vertices if cleanup flag is set
        if (DO_CLEAN_UP)
            tri::Clean<MyMesh>::RemoveUnreferencedVertex(Mesh);
        tri::UpdateFlags<MyMesh>::FaceBorderFromNone(Mesh); // Initialize border flags
        tri::UpdateTopology<MyMesh>::FaceFace(Mesh); // Compute information for face-to-face adjacency
        tri::UpdateTopology<MyMesh>::VertexFace(Mesh); // Compute information for vertex-to-face adjacency

        // Fill holes in the mesh using ear cutting algorithm with minimum weight criterion
        tri::Hole<MyMesh>::EarCuttingFill<tri::MinimumWeightEar< MyMesh> >(Mesh,500,false,nullptr);
        assert(tri::Clean<MyMesh>::IsFFAdjacencyConsistent(Mesh));

        // update mesh topology information
        tri::UpdateNormal<MyMesh>::NormalizePerFaceByArea(Mesh);
        tri::UpdateTopology<MyMesh>::FaceFace(Mesh);
        tri::UpdateTopology<MyMesh>::VertexFace(Mesh);

        // Clear selection
        tri::UpdateSelection<MyMesh>::FaceClear(Mesh);
        tri::UpdateSelection<MyMesh>::VertexClear(Mesh);

        outputFilename = meshesFolder / "HoleFilled"/ meshName;
        ExportMesh<MyMesh>(Mesh, outputFilename);


        size_t vertexSelected = 0;
        size_t faceIndex = originalFaceNumber;
        // First, select all new faces
        for (; faceIndex < Mesh.FN(); ++faceIndex)
        {
            if (Mesh.face[faceIndex].IsD())
                continue;
            Mesh.face[faceIndex].SetS();

        }

        std::cout << "number of face selected " << faceIndex - originalFaceNumber << std::endl;
        std::cout << "number of vertex selected " << vertexSelected << std::endl;
        execute_mesh_refinement(Mesh,0.01555);

        tri::UpdateSelection<MyMesh>::VertexClear(Mesh);
        tri::UpdateNormal<MyMesh>::NormalizePerFaceByArea(Mesh);
        tri::UpdateTopology<MyMesh>::FaceFace(Mesh);
        tri::UpdateTopology<MyMesh>::VertexFace(Mesh);

        // Then select vertices that have all adjacent faces selected
        faceIndex = originalFaceNumber;
        // First, select all new faces
        for (; faceIndex < Mesh.FN(); ++faceIndex)
        {
            if (Mesh.face[faceIndex].IsD())
                continue;

            for (int v=0 ; v<Mesh.face[faceIndex].VN(); ++v)
            {
                auto vertex = Mesh.face[faceIndex].V(v);
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

                std::cout << "vertex: " << vertex<< " | "<< faceCount<< " faces not selected " << std::endl;
                if (allFacesSelected)
                {
                    vertex->SetS();
                    vertexSelected++;
                }
            }
        }

        std::cout << "number of vertex selected " << vertexSelected << std::endl;

        outputFilename = meshesFolder / "HoleFilledAndRefined"/ meshName;
        ExportMesh<MyMesh>(Mesh, outputFilename);

        // Apply cotangent smoothing
        for (int i =0 ; i<3; i++)
        {
            MyMesh meshCopy;
            TypeOfWeight tow = static_cast<TypeOfWeight>(i);
            std::string weightTypeSuffix = getWeightTypeSuffix(tow);

            CloneMesh(Mesh, meshCopy);

            auto start = std::chrono::high_resolution_clock::now();
            HoleFairing<MyMesh>(meshCopy, 100,1, true, tow);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end - start;

            std::cout << "HoleFairing Duration (" << weightTypeSuffix << "): " << duration.count() << " seconds" <<
                std::endl;
            outputFilename = meshesFolder / ("HoleFilledAndRefined"+ weightTypeSuffix) / meshName;
            ExportMesh<MyMesh>(meshCopy, outputFilename);
        }
    }
    return 0;
}
