/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/*! \file trimesh_base.cpp
\ingroup code_sample

\brief the minimal example of using the lib

This file contain a minimal example of the library, showing how to load a mesh and how to compute per vertex normals on it.

*/
#include<unordered_map>
#include<vcg/complex/complex.h>
#include<wrap/io_trimesh/import_obj.h>
#include<wrap/io_trimesh/import_off.h>
#include<wrap/io_trimesh/export_off.h>
#include<vcg/space/polygon3.h>
using namespace vcg;

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public UsedTypes<Use<MyVertex>   ::AsVertexType,
                                           Use<MyEdge>     ::AsEdgeType,
                                           Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face<   MyUsedTypes, face::PolyInfo, face::PFVAdj, face::BitFlags > {};
class MyEdge    : public Edge<   MyUsedTypes> {};

class MyMesh    : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge>  > {};

int main( int argc, char **argv )
{
  if(argc<2)
  {
    printf("Usage trimesh_base <meshfilename.off>\n");
    return -1;
  }
  /*!
    */
  MyMesh m;
  MyMesh out;
  
  // for the output mesh we keep a counter for each vertex as an additional attribute
  auto v_cnH = tri::Allocator<MyMesh>::AddPerVertexAttribute<int>(out,"counter");



  int mask;
  if(tri::io::ImporterOBJ<MyMesh>::Open(m,argv[1],mask)!=tri::io::ImporterOFF<MyMesh>::NoError)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }

  tri::RequirePerVertexNormal(m);
  tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
  printf("Input mesh  vn:%i fn:%i\n",m.VN(),m.FN());

// step 0 copy all the vertexes in the output mesh
  for(auto &v : m.vert)
    tri::Allocator<MyMesh>::AddVertex(out,v.P());
printf( "Mesh has %i vert and %i faces\n", out.VN(), out.FN() );

    std::map<int, int> faceMap;
// for each face we add a new vertex in the center of the face
  for(auto &f : m.face)
  {
    Point3f center = PolyBarycenter(f);    
    auto vp = tri::Allocator<MyMesh>::AddVertex(out,center);
    faceMap[tri::Index(m,f)] = tri::Index(out,*vp);
  }
  
  printf( "Mesh has %i vert and %i faces\n", out.VN(), out.FN() );

  // second step. Create a vertex for each edge.
  // we use and unordered_map to store the vertex already created to avoid duplication

  std::map<std::pair<int,int>, int> edgeMap;
  
  // loop over all faces
  for(auto &f : m.face)
  {
    // loop over all edges
    for(int i=0;i<f.VN();i++)
    {
      int v0 = tri::Index(m,f.V0(i));
      int v1 = tri::Index(m,f.V1(i));
      if(v0>v1) std::swap(v0,v1);
      Point3f center = (f.V(i)->P() + f.V((i+1)%f.VN())->P()) / 2.0f;
       
      // check if the edge is already in the map
      auto it = edgeMap.find(std::make_pair(v0,v1));
      int edgeVertIndex;
      if(it == edgeMap.end())
      {
        // if not, create a new vertex in the middle of the edge
        auto vp = tri::Allocator<MyMesh>::AddVertex(out,center*6.0f);        
        edgeMap[std::make_pair(v0,v1)] = tri::Index(out,*vp);
        edgeVertIndex = tri::Index(out,*vp);
      }
      else
      { 
        edgeVertIndex = it->second;
        out.vert[edgeVertIndex].P() += center*6.0f;
      }
      
       out.vert[edgeVertIndex].P() += f.V2(i)->P();
       out.vert[edgeVertIndex].P() += f.V( (i+f.VN())%f.VN() )->P();

      v_cnH[edgeVertIndex]+=8; // increment the counter of the vertex for the edge  
      
    }
  }

printf( "Mesh has %i vert and %i faces\n", out.VN(), out.FN() );

  // third step. Create a face for each wedge of each face of the original mesh
  for(auto &f : m.face)
  {
    // loop over all wedge of the face 
    for(int i=0;i<f.VN();i++)
    {
      int e00 = tri::Index(m,f.V0(i));
      int e01 = tri::Index(m,f.V1(i));
      int e10 = tri::Index(m,f.V0(i));
      int e11 = tri::Index(m,f.V((i+f.VN()-1)%f.VN()));
      if(e00>e01) std::swap(e00,e01);
      if(e10>e11) std::swap(e10,e11);

      // retrieve the edge vertices 
      auto e1it = edgeMap.find(std::make_pair(e00,e01));
      auto e2it = edgeMap.find(std::make_pair(e10,e11));

      // retrieve the face vertex
      auto fit = faceMap.find(tri::Index(m,f));
      
      auto v0p = &out.vert[tri::Index(m,f.V(i))];
      auto v1p = &out.vert[e1it->second];
      auto v2p = &out.vert[fit->second];
      auto v3p = &out.vert[e2it->second];
      

      tri::Allocator<MyMesh>::AddQuadFace(out, v0p, v1p, v2p, v3p);
    }
  }
  printf( "Mesh has %i vert and %i faces\n", out.VN(), out.FN() );

  // final step; Normalize all vertices position using the handle counter
  for(auto &v : out.vert)
    if(v_cnH[&v]>0) v.P() /= v_cnH[&v];

  tri::io::ExporterOFF<MyMesh>::Save(out,"output.off",tri::io::Mask::IOM_VERTNORMAL);
  return 0;
}
