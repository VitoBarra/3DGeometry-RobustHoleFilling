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
#include<wrap/io_trimesh/import_off.h>
#include<vcg/complex/algorithms/clean.h>

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public vcg::UsedTypes<vcg::Use<MyVertex>   ::AsVertexType,
                                           vcg::Use<MyEdge>     ::AsEdgeType,
                                           vcg::Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags  >{};
class MyFace    : public vcg::Face<   MyUsedTypes, vcg::face::FFAdj,  vcg::face::VertexRef, vcg::face::BitFlags > {};
class MyEdge    : public vcg::Edge<   MyUsedTypes> {};

class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};

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

  if(vcg::tri::io::ImporterOFF<MyMesh>::Open(m,argv[1])!=0)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }

  printf("Input mesh  vn:%i fn:%i\n",m.VN(),m.FN());
  // Counting the number of edges using FF adjacency
  vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
  int boundaryEdgesNum=0;
  int boundaryLoopsNum=0;

  vcg::tri::UpdateFlags<MyMesh>::FaceClearV(m);

  for(auto &f : m.face)
  {
    for(int i=0;i<3;++i)
    {
      if(f.FFp(i)==&f)
      {      
        boundaryEdgesNum++;
        
        if(!f.IsV())
        { // Use the face::Pos  to navigate around the boundary loop
          boundaryLoopsNum++;
          int boundaryEdgesNum=0;
          vcg::face::Pos<MyFace> pos(&f,i,f.V(i));
          vcg::face::Pos<MyFace> start=pos;
          // Loop around the boundary
          do 
          { 
            assert(pos.IsBorder());
            assert(!pos.F()->IsV());

            pos.F()->SetV();
            // Loop around the vertex
            do {
              pos.FlipE();
              pos.FlipF();
            } while (!pos.IsBorder());
            boundaryEdgesNum++;

            pos.FlipV();
          } while(pos !=start);
          printf("Boundary size %i\n",boundaryEdgesNum);
  
        }
      }

    }        
  }

  printf("Boundary loops %i\n",boundaryLoopsNum);

  int EN = (m.fn *3 + boundaryEdgesNum)/2;
  printf("Euler characteristic %i\n", m.vn - EN + m.fn);
  
  return 0;
}
