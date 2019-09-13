#include "PlyReader.h"
#include "UnstructuredGrid3D.h"
#include "ply.h"
#include <string.h>
#include <iostream>


typedef struct Vertex 
{
  float x,y,z;             /* the usual 3-space position of a vertex */
} Vertex;

typedef struct Face 
{
  unsigned char nverts;    /* number of vertex indices in list */
  int *verts;              /* vertex index list */
} Face;


PlyProperty vert_props[] = { /* list of property information for a vertex */
  {"x", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,x), 0, 0, 0, 0},
  {"y", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,y), 0, 0, 0, 0},
  {"z", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,z), 0, 0, 0, 0},
};

PlyProperty face_props[] = { /* list of property information for a vertex */
  {"vertex_indices", PLY_INT, PLY_INT, offsetof(Face,verts),
   1, PLY_UCHAR, PLY_UCHAR, offsetof(Face,nverts)},
};




UnstructuredGrid3D* PlyReader::read(const char* filename)
{
  int i,j,k;
  int nelems;
  char **elist;
  int file_type;
  float version;
  int nprops;
  int num_elems;
  PlyProperty **plist;
  char *elem_name;

  PlyFile* ply = ply_open_for_reading((char*)filename, &nelems, &elist, &file_type, &version);		//Open file for reading

  int num_vertices, num_faces;
 
  ply_get_element_description(ply, "vertex", &num_vertices, &nprops);
  plist = ply_get_element_description (ply, "face", &num_faces, &nprops);
  
  UnstructuredGrid3D* grid = new UnstructuredGrid3D(num_vertices,num_faces);

  for (i = 0; i < nelems; ++i) 
  {
    elem_name = elist[i];
    plist = ply_get_element_description(ply, elem_name, &num_elems, &nprops);

    if (equal_strings("vertex", elem_name)) 
	{
      ply_get_property(ply, elem_name, &vert_props[0]);
      ply_get_property(ply, elem_name, &vert_props[1]);
      ply_get_property(ply, elem_name, &vert_props[2]);

     for (j = 0; j < num_vertices; j++) 
	  {
		Vertex v; float V[3]; 
        ply_get_element(ply,&v);
		V[0]=v.x; V[1]=v.y; V[2]=v.z;
		grid->setPoint(j,V);
      }
    }

    if (equal_strings("face", elem_name)) 
	{
      ply_get_property (ply, elem_name, &face_props[0]);

      for (j = 0; j < num_faces; j++) 
	  {
	    Face face;
        ply_get_element(ply,&face);		
		grid->setCell(j,face.verts);
      }
    }    
  }
  
  ply_close (ply);
  return grid;
}
