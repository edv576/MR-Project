#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <list>
#include <math.h>
#include <algorithm>
#include "OFF_PLYConverter.h"




OFF_PLYConverter::OFF_PLYConverter() {


}

OFF_PLYConverter::~OFF_PLYConverter() {


}

void OFF_PLYConverter::Convert_OFF_PLY(FILE *fo, FILE *fd){

	char buffer[101];
	std::list<Point> points;
	std::list<Face> faces;

	//Scan first line - type of file
	fscanf(fo, "%100s\n", buffer);
	//Temporal value for the number of vertices
	int nv;
	//Temporal value for the number of faces
	int nf;
	//Temporal value for the number of edges
	int ne;
	//Temporal values for an edge
	float x, y, z;
	//Temporal values for type of faces a points of a face
	int tf, p1, p2, p3;
	//Point to store the centroid (barycenter)
	Point centroid;
	//Initializing the centroid in (0, 0, 0)
	centroid.x = 0;
	centroid.y = 0;
	centroid.z = 0;
	//Here the values for the translation will be stored
	Point correction;
	//Variables for the minimum values of x, y and z
	float minX, minY, minZ;
	//Variables for the maximum values of x, y and z
	float maxX, maxY, maxZ;

	if (strcmp(buffer, "OFF")) {

		return;
	}
	else
	{
	
		//Scan number of vertex, faces and edges
		fscanf(fo, "%d %d %d", &nv, &nf, &ne);


		for (int i = 0; i < nv; i++)
		{
			Point point;
			fscanf(fo, "%f %f %f", &x, &y, &z);
			point.x = x;
			point.y = y;
			point.z = z;
			//Adding all vertices to a list
			points.push_back(point);
			//Adding the values of a point to the Centroid
			centroid.x += point.x;
			centroid.y += point.y;
			centroid.z += point.z;			
		}

		std::list<Point>::iterator itPoints = points.begin();

		for (int i = 0; i < nf; i++)
		{
			Face face;
			fscanf(fo, "%d %d %d %d", &tf, &p1, &p2, &p3);
			face.type_face = tf;
			face.point1 = p1;
			face.point2 = p2;
			face.point3 = p3;
			//Adding the values of the faces to a list
			faces.push_back(face);			
		}

		std::list<Face>::iterator itFaces = faces.begin();

		//Calculating the actual Centroid
		centroid.x /= nv;
		centroid.y /= nv;
		centroid.z /= nv;

		//Initializing the values for the correction of the Centroid. Needed to put it in (0,0,0)
		correction.x = -1 * centroid.x;
		correction.y = -1 * centroid.y;
		correction.z = -1 * centroid.z;


		//Print header in the ply file
		fprintf(fd, "ply\n");
		fprintf(fd, "format ascii 1.0\n");
		fprintf(fd, "comment %s\n", " ");
		fprintf(fd, "element vertex %d\n", nv);
		fprintf(fd, "property float x\n");
		fprintf(fd, "property float y\n");
		fprintf(fd, "property float z\n");
		fprintf(fd, "element face %d\n", nf);
		fprintf(fd, "property list uchar int vertex_indices\n");
		fprintf(fd, "end_header\n");


		itPoints = points.begin();

		float ex = 1.0e6;

		//Getting the first values for the minimum and maximum values
		//minX = itPoints->x;
		//minY = itPoints->y;
		//minZ = itPoints->z;
		//maxX = itPoints->x;
		//maxY = itPoints->y;
		//maxZ = itPoints->z;

		minX = ex;
		minY = ex;
		minZ = ex;
		maxX = -ex;
		maxY = -ex;
		maxZ = -ex;


		//Correcting the Centroid and moving every vertex
		//It also calculates the minimum and maximum values
		for (int i = 0; i < nv; i++)
		{

			itPoints = points.begin();
			std::advance(itPoints, i);
			itPoints->x = itPoints->x + correction.x;
			itPoints->y = itPoints->y + correction.y;
			itPoints->z = itPoints->z + correction.z;
			
			minX = min(itPoints->x, minX);
			minY = min(itPoints->y, minY);
			minZ = min(itPoints->z, minZ);
			maxX = max(itPoints->x, maxX);
			maxY = max(itPoints->y, maxY);
			maxZ = max(itPoints->z, maxZ);

			//if (itPoints->x < minX) {
			//	minX = itPoints->x;				
			//}
			//if (itPoints->y < minY) {
			//	minY = itPoints->y;
			//}
			//if (itPoints->z < minZ) {
			//	minZ = itPoints->z;
			//}
			//if (itPoints->x > maxX) {
			//	maxX = itPoints->x;
			//}
			//if (itPoints->y > maxY) {
			//	maxY = itPoints->y;
			//}
			//if (itPoints->z > maxZ) {
			//	maxZ = itPoints->z;
			//}

			
		}

		float sizeX = maxX - minX;								//2. Compute a single scaling factor that best fits the grid
		sizeX = (sizeX) ? 1 / sizeX : 1;							//   in the [-1,1] cube. Using a single factor for x,y, and z
		float sizeY = maxY - minY;								//   ensures that the object is scaled while keeping its
		sizeY = (sizeY) ? 1 / sizeY : 1;							//   aspect ratio.
		float sizeZ = maxZ - minZ;
		sizeZ = (sizeZ) ? 1 / sizeZ : 1;

		float scale = min(sizeX, min(sizeY, sizeZ));

		//Initializing the extreme points for normalization. Bounding box of length 1.
		//Point extreme1;
		//Point extreme2;

		//extreme1.x = -0.5;
		//extreme1.y = -0.5;
		//extreme1.z = -0.5;

		//extreme2.x = 0.5;
		//extreme2.y = 0.5;
		//extreme2.z = 0.5;

		

		//extreme1.x = -ex;
		//extreme1.y = -ex;
		//extreme1.z = -ex;

		//extreme2.x = ex;
		//extreme2.y = ex;
		//extreme2.z = ex;

		//float minT;
		//float scale;

		////Getting the scale
		//minT = std::min(1 / (maxX - minX), 1 / (maxY - minY));
		//scale = std::min(minT, 1 / (maxZ - minZ));

		for (int i = 0; i < nv; i++)
		{

			itPoints = points.begin();
			float xt, yt, zt;
			std::advance(itPoints, i);

			//itPoints->x = 2.0*(itPoints->x - minX) / (maxX - minX) - 1.0;
			//itPoints->y = 2.0*(itPoints->y - minY) / (maxY - minY) - 1.0;
			//itPoints->z = 2.0*(itPoints->z - minZ) / (maxZ - minZ) - 1.0;

			//Doing the scaling of the mesh. Its done keeping the aspect ratio
			itPoints->x = 2 * ((itPoints->x - minX)*scale - 0.5);
			itPoints->y = 2 * ((itPoints->y - minY)*scale - 0.5);
			itPoints->z = 2 * ((itPoints->z - minZ)*scale - 0.5);
			//itPoints->x = (itPoints->x - 0.5*(minX + maxX))*scale;
			//itPoints->y = (itPoints->y - 0.5*(minY + maxY))*scale;
			//itPoints->z = (itPoints->z - 0.5*(minZ + maxZ))*scale;

			xt = itPoints->x;
			yt = itPoints->y;
			zt = itPoints->z;

			//Writing the vertices values in the .ply file
			fprintf(fd, "%f %f %f\n", xt, yt, zt);

		}

		for (int i = 0; i < nf; i++)
		{
			itFaces = faces.begin();
			int type_face_t, point1_t, point2_t, point3_t;
			std::advance(itFaces, i);
			type_face_t = itFaces->type_face;
			point1_t = itFaces->point1;
			point2_t = itFaces->point2;
			point3_t = itFaces->point3;
			//Writing the faces values in the .ply file
			fprintf(fd, "%d %d %d %d\n", type_face_t, point1_t, point2_t, point3_t);
		}

		

	}

}
