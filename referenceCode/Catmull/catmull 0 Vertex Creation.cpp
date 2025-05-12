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
#include<wrap/io_trimesh/import_off.h>
#include<wrap/io_trimesh/export_off.h>
#include<vcg/space/polygon3.h>
using namespace vcg;

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public UsedTypes<Use<MyVertex>   ::AsVertexType,
                                           Use<MyEdge>     ::AsEdgeType,
                                           Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face<   MyUsedTypes, face::FFAdj,  face::VertexRef, face::BitFlags > {};
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

  if(tri::io::ImporterOFF<MyMesh>::Open(m,argv[1])!=tri::io::ImporterOFF<MyMesh>::NoError)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }

  tri::RequirePerVertexNormal(m);
  tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
  printf("Input mesh  vn:%i fn:%i\n",m.VN(),m.FN());

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
      // check if the edge is already in the map
      auto it = edgeMap.find(std::make_pair(v0,v1));
      if(it == edgeMap.end())
      {
        // if not, create a new vertex in the middle of the edge
        Point3f center = (f.V(i)->P() + f.V((i+1)%f.VN())->P()) / 2.0f;
        auto vp = tri::Allocator<MyMesh>::AddVertex(out,center);        
        edgeMap[std::make_pair(v0,v1)] = tri::Index(out,*vp);
      }
    }
  }

printf( "Mesh has %i vert and %i faces\n", out.VN(), out.FN() );


  tri::io::ExporterOFF<MyMesh>::Save(out,"output.off",tri::io::Mask::IOM_VERTNORMAL);
  return 0;
}
