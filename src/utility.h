//
// Created by vbvit on 25/05/2025.
//
#pragma once

#include <vcg/complex/complex.h>
#include <cmath>

#ifndef STRUCTURE_H
#define STRUCTURE_H


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


// Parameters:
//   m - Input mesh to be smoothed
//   step - Number of smoothing iterations
//   alpha - Smoothing factor (0-1), controls the intensity of smoothing
//   SmoothSelected - If true, only smooth selected vertices
template <class MeshType>
static void CloneMesh(MeshType& source, MeshType& target)
{
    vcg::tri::Append<MeshType, MeshType>::MeshCopy(target, source);
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

// calculate the cotangent function, it is stable around 0 and pi but can have some instabilities around pi/2
float cot(float angle)
{
    return tan(M_PI *0.5 - angle);
}



#endif //STRUCTURE_H
