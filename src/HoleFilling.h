//
// Created by vbvit on 01/06/2025.
//

#ifndef HOLEFILLING_H
#define HOLEFILLING_H
#include "vcg/space/point.h"
#include <cmath>
#include <vector>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/simplex/face/topology.h>


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

    std::map<int, float> meanEdgeLength;

    // Compute the average edges length for each boundary vertex ignoring internal edges
    for (VertexType &v : mesh.vert)
    {
        if (v.IsD() || !v.IsS()) continue;
        float sumLen = 0.0f;
        float count = 0;
        std::vector<FaceType *> faces;
        std::vector<int> verIND;

        // compute the VV adjacency
        vcg::face::VFStarVF<FaceType>(&v,faces,verIND) ;

        // Calculate average edge length using VV adjacency
        for (FaceType* f : faces)
        {
            if (f->IsD() || f->IsS()) continue;
            for (int j = 0; j < 3; ++j)
            {
                if (f->V(j) == &v) continue;
                sumLen += Distance(v.P(), f->V(j)->P());
                count++;
            }
        }

        // Store average edge length in scale attribute
        float meanLength =(count > 0) ? sumLen / count : std::numeric_limits<float>::infinity();
        meanEdgeLength[tri::Index(mesh,v)] = meanLength;
        v.Q()= meanLength;
    }

    bool updated = true;
    while (true)
    {
        step_counter++;
        updated = false;

        std::vector<std::tuple<int,int, int,int >> faceToAdd;
        for (FaceType &f : mesh.face)
        {
            if (f.IsD() || !f.IsS() ) continue;

            // Compute centroid
            vcg::Point3f centroidPos = (f.V(0)->P() + f.V(1)->P() + f.V(2)->P()) / 3.0f;
            const float centroidMeanLength = ( meanEdgeLength[tri::Index(mesh,f.V(0))] +meanEdgeLength[tri::Index(mesh,f.V(1))] + meanEdgeLength[tri::Index(mesh,f.V(2))]) / 3.0f;



            bool needToAddTriangle = true;
            for (int j = 0; j < 3; ++j)
            {
                float weightedDistance = densityFactor*Distance(centroidPos, f.V(j)->P());
                if (!(weightedDistance > centroidMeanLength && weightedDistance > meanEdgeLength[tri::Index(mesh,f.V(j))])) {
                    needToAddTriangle = false;
                    break;
                }
            }

            if (needToAddTriangle)
            {
                updated = true;
                VertexType* centroidVert = &*tri::Allocator<MeshType>::AddVertex(mesh,centroidPos);
                meanEdgeLength[tri::Index(mesh,centroidVert)] = centroidMeanLength;
                centroidVert->Q()= centroidMeanLength;

                // save vertexes indices
                int v0fi = tri::Index(mesh,f.V(0)); //vertex 0 face index
                int v1fi = tri::Index(mesh,f.V(1));
                int v2fi = tri::Index(mesh,f.V(2));
                int vci = tri::Index(mesh,centroidVert) ; // vertex centroid index
                faceToAdd.push_back(std::make_tuple(v0fi, v1fi, v2fi, vci));

                // Delete the original face after storing vertices
                tri::Allocator<MeshType>::DeleteFace(mesh, f);
            }
        }

        tri::Allocator<MeshType>::CompactFaceVector(mesh);

        //Add faces
        for (const auto &[v0fi,v1fi,v2fi,vci] : faceToAdd)
        {
            VertexType* v0f = &mesh.vert[v0fi];
            VertexType* v1f = &mesh.vert[v1fi];
            VertexType* v2f = &mesh.vert[v2fi];
            VertexType* v   = &mesh.vert[vci];


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
        faceToAdd.clear();

        //Exit if not face is added
        if (!updated)
            return;

        ExportMeshInFolder<MeshType>(mesh, "HoleFilledPreSwap_step"+std::to_string(step_counter));

        //Relax all interior edges
        for (int retry = 0; retry<3;)
        {
            int edgeSwapped = 0;
            //recompute adjacency
            tri::UpdateTopology<MeshType>::FaceFace(mesh);
            tri::UpdateTopology<MeshType>::VertexFace(mesh);
            for (FaceType &f : mesh.face)
            {
                if (f.IsD() || !f.IsS()) continue;

                for (int edge = 0; edge < 3; ++edge)
                {
                    // Get the adjacent face across the edge
                    FaceType *adjF = f.FFp(edge);
                    int adjEdgeIdx = f.FFi(edge);

                    if (adjF->IsD() || !adjF->IsS() || adjF ==nullptr || adjF == &f   ) continue;

                    // Get shared-edge vertices
                    VertexType *v0 = f.V0(edge);
                    VertexType *v1 = f.V1(edge);

                    // Get opposing vertices
                    VertexType *vOppF   = f.V2(edge);
                    VertexType *vOppAdj = adjF->V2(adjEdgeIdx);

                    // Compute circumcircle of the opposite triangle (v0, v1, vOppF)
                    Triangle3 tri(v0->P(), v1->P(), vOppF->P());
                    Point3f circumcenter = Circumcenter(tri);
                    const float circumSphere_radius = Distance(circumcenter,vOppF ->P());
                    const float oppositeAdjVertexDistance = Distance(circumcenter, vOppAdj->P());


                    if (oppositeAdjVertexDistance < circumSphere_radius) {
                        face::FlipEdge<FaceType>(f, edge);
                        edgeSwapped++;
                    }
                }
            }
            if (edgeSwapped==0) // if no edge is swapped, there is nothing left to do
                break;
            if (edgeSwapped<5)   // if some few edges were swapped, maybe almost finished or the swap is in a loop
                retry++;
        }

        ExportMeshInFolder<MeshType>(mesh, "HoleFilledPostSwap_step"+std::to_string(step_counter));
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
