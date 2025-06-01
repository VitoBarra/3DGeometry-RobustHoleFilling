//
// Created by vbvit on 01/06/2025.
//

#ifndef HOLEFILLING_H
#define HOLEFILLING_H
#include "MeshIO.h"
#include "vcg/space/point.h"
#include <cmath>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/complex.h>
#include <vcg/simplex/face/topology.h>
#include <vector>

using namespace vcg;
template<class VertexType>
struct portableFace
{
    VertexType v1;
    VertexType v2;
    VertexType v3;

};

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

// LaplacianInfo: Helper class for Laplacian smoothing operations
// Stores the sum of neighboring vertices' positions and the count of neighbors
// Used for calculating the average position in smoothing operations
template <typename MeshType>
class LaplacianInfo
{
public:
    LaplacianInfo(const MeshType::CoordType& _p, const int _n) : sum(_p), cnt(_n){}

    LaplacianInfo() = default;

    typename MeshType::CoordType sum; // Sum of neighboring vertices' positions
    typename MeshType::ScalarType cnt; // Count of neighboring vertices
};



template <class MeshType>
void HolePatchRefinement(MeshType &mesh, float densityFactor = 1.414213562373f) //square root 2
{
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::FaceType FaceType;
    int step_counter = 0;

    tri::UpdateTopology<MeshType>::FaceFace(mesh);
    tri::UpdateTopology<MeshType>::VertexFace(mesh);

    std::map<VertexType*, float> meanLengthMap;

    // Step 1: Compute the average edges length laundry for each boundary vertex
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
        {
            sumLen += Distance(v.P(), adj->P());
            count++;
        }

        // Store average edge length in scale attribute
        meanLengthMap[&v] = (count > 0) ? sumLen / count : 0.0f;
    }

    bool updated = true;
    while (true)
    {
        step_counter++;
        updated = false;

        // Step 2: For each triangle
        std::vector<std::tuple<VertexType*,VertexType*,VertexType*,VertexType*>> faceToAdd;
        for (int i = 0  ; i < mesh.face.size(); i++)
        {
            FaceType* f = &mesh.face[i];
            if (f->IsD() || !f->IsS() ) continue;

            vcg::Point3f vertexCentroid = (f->V(0)->P() + f->V(1)->P() + f->V(2)->P()) / 3.0f;
            // Compute centroid
            float c_MeanLength = (meanLengthMap[f->V(0)] + meanLengthMap[f->V(1)] + meanLengthMap[f->V(2)]) / 3.0f;

            bool needToAddTriangle = true;
            for (int i = 0; i < 3; ++i)
            {
                float weightedDistance = densityFactor*Distance(vertexCentroid, f->V(i)->P());
                if (weightedDistance <= c_MeanLength || weightedDistance <= meanLengthMap[f->V(i)]) {
                    needToAddTriangle = false;
                    break;
                }
            }

            if (needToAddTriangle)
            {
                updated = true;
                // Insert centroid vertex
                auto vC = tri::Allocator<MeshType>::AddVertex(mesh,vertexCentroid);
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
            auto f1 = tri::Allocator<MeshType>::AddFace(mesh, v, v1f, v2f);
            f1->SetS();

            // Second face: (v0f, v, v2f)
            auto f2 = tri::Allocator<MeshType>::AddFace(mesh, v0f, v, v2f);
            f2->SetS();

            // Third face: (v0f, v1f, v)
            auto f3 = tri::Allocator<MeshType>::AddFace(mesh, v0f, v1f, v);
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
                    Triangle3 tri(v0->P(), v1->P(), vOppF->P());
                    Point3f circumcenter = Circumcenter(tri);
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

            weight = calculateWeight<MeshType>(weightType, &(*fi), edge);

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
    typedef typename MeshType::VertContainer VertContainer;

    VertexType vi;
    LaplacianInfo<MeshType> lpz(CoordType(0, 0, 0), 0);
    assert(alpha <= 1.0f && alpha > 0.f);
    SimpleTempData<VertContainer, LaplacianInfo<MeshType>> TD(m.vert);

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


#endif //HOLEFILLING_H
