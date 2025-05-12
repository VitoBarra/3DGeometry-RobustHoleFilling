#include <iostream>

#include<vcg/complex/complex.h>
#include<wrap/io_trimesh/import_off.h>
#include<vcg/complex/algorithms/clean.h>



using namespace  vcg;
class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public UsedTypes<Use<MyVertex>   ::AsVertexType,
                                           Use<MyEdge>     ::AsEdgeType,
                                           Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face<   MyUsedTypes, face::FFAdj,  face::VertexRef, face::BitFlags > {};
class MyEdge    : public Edge<   MyUsedTypes> {};

class MyMesh    : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};

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

  if(tri::io::ImporterOFF<MyMesh>::Open(m,argv[1])!=0)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }

  printf("Input mesh  vn:%i fn:%i\n",m.VN(),m.FN());
  // Counting the number of edges using FF adjacency
  tri::UpdateTopology<MyMesh>::FaceFace(m);
  int boundaryEdgesNum=0;
  int boundaryLoopsNum=0;

  tri::UpdateFlags<MyMesh>::FaceClearV(m);

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
          face::Pos<MyFace> pos(&f,i,f.V(i));
          face::Pos<MyFace> start=pos;
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