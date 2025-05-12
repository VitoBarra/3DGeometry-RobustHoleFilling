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

#include<vcg/complex/complex.h>
#include<wrap/io_trimesh/import_stl.h>
#include<vcg/complex/algorithms/clean.h>

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public vcg::UsedTypes<vcg::Use<MyVertex>   ::AsVertexType,
                                           vcg::Use<MyEdge>     ::AsEdgeType,
                                           vcg::Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags  >{};
class MyFace    : public vcg::Face<   MyUsedTypes, vcg::face::FFAdj,  vcg::face::VertexRef, vcg::face::BitFlags > {};
class MyEdge    : public vcg::Edge<   MyUsedTypes> {};

class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge>  > {};

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

  // The Mask is a bitfield that tells the importer which kind of data has found in the file
  int LoadMask=0;

  if(vcg::tri::io::ImporterSTL<MyMesh>::Open(m,argv[1],LoadMask)!=0)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }

  vcg::tri::RequirePerVertexNormal(m);
  vcg::tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
  printf("Input mesh  vn:%i fn:%i\n",m.VN(),m.FN());

  vcg::tri::Clean<MyMesh>::RemoveDuplicateVertex(m);
  printf("After duplicate removal      vert vector size :%i face vector:%i\n",m.vert.size(),m.face.size());
 
  vcg::tri::Allocator<MyMesh>::CompactEveryVector(m);
  printf("After compacting the vectors vert vector size :%i face vector:%i\n",m.vert.size(),m.face.size());

  printf( "Mesh has %i vert and %i faces\n", m.VN(), m.FN() );

  return 0;
}
