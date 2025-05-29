#include <chrono>
#include <iostream>

#include <cmath>
#include <vector>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/isotropic_remeshing.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/simplex/face/topology.h>
#include <vcg/space/point3.h>
#include <vcg/space/triangle3.h>
#include <wrap/io_trimesh/import_ply.h>
#include "vcg/complex/algorithms/hole.h"


#include "file.h"
#include "utility.h"
#include "debug.h"


using namespace  vcg;
class MyVertex; class MyEdge; class MyFace;




// Base types definition for the mesh components
struct MyUsedTypes : public UsedTypes<Use<MyVertex>::AsVertexType,
                                      Use<MyEdge>::AsEdgeType,
                                      Use<MyFace>::AsFaceType>
{
};

// Vertex class with 3D coordinates, normals, vertex-face adjacency and flags
class MyVertex : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::VFAdj,vertex::VEAdj, vertex::BitFlags, vertex::Mark>
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



struct portableFace
{
    MyVertex v1;
    MyVertex v2;
    MyVertex v3;

};



template <class MeshType>
void HolePatchRefinement(MeshType &mesh, float densityFactor = 1.414213562373f) //square root 2
{
    using namespace vcg::tri;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::FaceType FaceType;

    UpdateTopology<MeshType>::FaceFace(mesh);
    UpdateTopology<MeshType>::VertexFace(mesh);

    std::map<VertexType*, float> meanLenghtMap;

    // Step 1: Compute the mean edge length laundry for each boundary vertex
    for (VertexType &v : mesh.vert)
    {
        if (v.IsD() || !v.IsS()) continue;
        float sumLen = 0.0f;
        float count = 0;
        std::vector<VertexType *> verts;

        // compute the VV adjacency 
        vcg::face::VVStarVF<FaceType>(&v,verts) ;

        // Calculate average edge length using VV adjacency
        for (VertexType* adj : verts)
        {        vcg::face::VVStarVF<FaceType>(&v,verts) ;

            sumLen += Distance(v.P(), adj->P());
            count++;
        }

        // Store average edge length in scale attribute
        meanLenghtMap[&v] = (count > 0) ? sumLen / count : 0.0f;
    }

    bool updated = true;
    while (true)
    {
        updated = false;

        // Step 2: For each triangle
        std::vector<std::tuple<VertexType*,VertexType*,VertexType*,VertexType*>> faceToAdd;
        for (int i = 0  ; i < mesh.face.size(); i++)
        {
            FaceType* f = &mesh.face[i];
            if (f->IsD() || !f->IsS() ) continue;

            Point3f vertexCentroid = (f->V(0)->P() + f->V(1)->P() + f->V(2)->P()) / 3.0f;
            // Compute centroid
            float c_MeanLength = (meanLenghtMap[f->V(0)] + meanLenghtMap[f->V(1)] + meanLenghtMap[f->V(2)]) / 3.0f;

            bool needToAddTriangle = true;
            for (int i = 0; i < 3; ++i)
            {
                float weightedDistance = densityFactor*Distance(vertexCentroid, f->V(i)->P());
                if (weightedDistance <= c_MeanLength || weightedDistance <= meanLenghtMap[f->V(i)]) {
                    needToAddTriangle = false;
                    break;
                }
            }

            if (needToAddTriangle)
            {
                updated = true;
                // Insert centroid vertex
                auto vC = Allocator<MeshType>::AddVertex(mesh,vertexCentroid);
                VertexType* v = &*vC;

                // Replace triangle with 3 new triangles
                VertexType* v0f = f->V(0);
                VertexType* v1f = f->V(1);
                VertexType* v2f = f->V(2);
                faceToAdd.push_back(std::make_tuple(v0f, v1f, v2f, v));

                // Delete the original face after storing vertices
                tri::Allocator<MeshType>::DeleteFace(mesh, *f);
            }
        }
        tri::Allocator<MeshType>::CompactFaceVector(mesh);

        for (const auto &face : faceToAdd)
        {
            VertexType *v0f = std::get<0>(face);
            VertexType *v1f = std::get<1>(face);
            VertexType *v2f = std::get<2>(face);
            VertexType *v   = std::get<3>(face);
            // first face:  (v, v1f, v2f)
            auto f1 = Allocator<MeshType>::AddFace(mesh, v, v1f, v2f);
            f1->SetS();

            // Second face: (v0f, v, v2f)
            auto f2 = Allocator<MeshType>::AddFace(mesh, v0f, v, v2f);
            f2->SetS();

            // Third face: (v0f, v1f, v)
            auto f3 = Allocator<MeshType>::AddFace(mesh, v0f, v1f, v);
            f3->SetS();
        }

        if (!updated)
            return;

        //recompute adjacency
        tri::UpdateTopology<MeshType>::FaceFace(mesh);
        tri::UpdateTopology<MeshType>::VertexFace(mesh);


        bool swappedEdge = true;
        while (swappedEdge)
        {
            tri::UpdateTopology<MeshType>::FaceFace(mesh);
            tri::UpdateTopology<MeshType>::VertexFace(mesh);
            // Step 4: Relax all interior edges
            swappedEdge= false;
            for (FaceType &f : mesh.face)
            {
                if (f.IsD() || !f.IsS()) continue;

                for (int edge = 0; edge < 3; ++edge)
                {
                    // Get the adjacent face across edge i
                    FaceType *adjF = f.FFp(edge);
                    int adjEdgeIdx = f.FFi(edge);

                    if (adjF ==nullptr || adjF == &f || adjF->IsD()) continue;

                    // Get shared edge vertices
                    VertexType *v0 = f.V0(edge);
                    VertexType *v1 = f.V1(edge);

                    // Get opposing vertices
                    VertexType *vOppF = f.V2(edge);
                    VertexType *vOppAdj = adjF->V2(adjEdgeIdx);

                    // Compute circumcircle of triangle (v0, v1, vOppF)
                    vcg::Triangle3 tri(v0->P(), v1->P(), vOppF->P());
                    vcg::Point3f circumcenter = Circumcenter(tri);
                    float circumcenter_radius= Distance(circumcenter, v0->P());

                    // Print face references

                    // std::cout << "------------------------------------------------------------------"<< std::endl;
                    // std::cout << "Face   1: " << &f << " | Face   2: " << adjF << std::endl;
                    // std::cout << "vertex 0: " << f.cV0(edge) << " | vertex 0: " << adjF->cV0(edge) << std::endl;
                    // std::cout << "vertex 1: " << f.cV1(edge) << " | vertex 1: " << adjF->cV1(edge) << std::endl;
                    // std::cout << "vertex 2: " << f.cV2(edge) << " | vertex 2: " << adjF->cV2(edge) << std::endl;
                    // std::cout << std::flush;



                    if (Distance(circumcenter, vOppAdj->P()) < circumcenter_radius) {
                        try
                        {
                            face::FlipEdge<FaceType>(f, edge);
                            swappedEdge= true;
                        }
                        catch (std::exception &e)
                        {
                            tri::UpdateSelection<MeshType>::FaceClear(mesh);
                            f.SetS();
                            adjF->SetS();
                            ExportMeshInFolder(mesh, "failed");
                            //throw e;
                        }
                    }
                }


            }


        }

    }
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



template <class MeshType>
static void HoleFairing(MeshType &m, int step, float alpha, bool SmoothSelected = false , TypeOfWeight weightType = uniform)
{
    typedef typename MeshType::VertexIterator VertexType;
    typedef typename MeshType::CoordType CoordType;

    VertexType vi;
    LaplacianInfo<MyMesh> lpz(CoordType(0, 0, 0), 0);
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
                    CoordType Delta = TD[*vi].sum / TD[*vi].cnt - (*vi).P();
                    (*vi).P() = (*vi).P() + Delta * alpha;
                }
            }
    }
}




int main( int argc, char **argv )
{
    for (auto &v : GetFilesInFolderExtension(meshesFolder ,".ply"))
    {
        MyMesh mesh;
        meshName =fs::path(v).filename().string();

        if(tri::io::ImporterPLY<MyMesh>::Open(mesh,v.c_str())!=0)
        {
            printf("Error reading file  %s\n",meshName);
            continue;
        }
        size_t originalFaceNumber= mesh.FN();
        printf("mesh %s has vertexes:%i faces:%i\n",meshName.c_str(),mesh.VN(),originalFaceNumber);
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
        tri::UpdateNormal<MyMesh>::NormalizePerFaceByArea(mesh);
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
        HolePatchRefinement(mesh);
        // Compact mesh to remove deleted elements
        tri::Allocator<MyMesh>::CompactFaceVector(mesh);
        tri::Allocator<MyMesh>::CompactVertexVector(mesh);
        tri::UpdateSelection<MyMesh>::VertexClear(mesh);
        tri::UpdateNormal<MyMesh>::NormalizePerFaceByArea(mesh);
        tri::UpdateTopology<MyMesh>::FaceFace(mesh);
        tri::UpdateTopology<MyMesh>::VertexFace(mesh);

        // Then select vertices that have all adjacent faces selected
        faceIndex = originalFaceNumber;
        // First, select all new faces
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

