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

	fscanf(fo, "%100s\n", buffer);
	int nv;
	int nf;
	int ne;
	float x, y, z;
	int tf, p1, p2, p3;
	Point centroid;
	centroid.x = 0;
	centroid.y = 0;
	centroid.z = 0;
	Point correction;
	float minX, minY, minZ;
	float maxX, maxY, maxZ;

	if (strcmp(buffer, "OFF")) {

		return;
	}
	else
	{
	
		fscanf(fo, "%d %d %d", &nv, &nf, &ne);


		for (int i = 0; i < nv; i++)
		{
			Point point;
			fscanf(fo, "%f %f %f", &x, &y, &z);
			point.x = x;
			point.y = y;
			point.z = z;
			points.push_back(point);
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
			faces.push_back(face);			
		}

		std::list<Face>::iterator itFaces = faces.begin();

		centroid.x /= nv;
		centroid.y /= nv;
		centroid.z /= nv;

		correction.x = -1 * centroid.x;
		correction.y = -1 * centroid.y;
		correction.z = -1 * centroid.z;


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
		minX = itPoints->x;
		minY = itPoints->y;
		minZ = itPoints->z;
		maxX = itPoints->x;
		maxY = itPoints->y;
		maxZ = itPoints->z;

		for (int i = 0; i < nv; i++)
		{

			itPoints = points.begin();
			std::advance(itPoints, i);
			itPoints->x = itPoints->x + correction.x;
			itPoints->y = itPoints->y + correction.y;
			itPoints->z = itPoints->z + correction.z;
			if (itPoints->x < minX) {
				minX = itPoints->x;
			}
			if (itPoints->y < minY) {
				minY = itPoints->y;
			}
			if (itPoints->z < minZ) {
				minZ = itPoints->z;
			}
			if (itPoints->x > maxX) {
				maxX = itPoints->x;
			}
			if (itPoints->y > maxY) {
				maxY = itPoints->y;
			}
			if (itPoints->z > maxZ) {
				maxZ = itPoints->z;
			}

			
		}

		//int j = 0;
		Point extreme1;
		Point extreme2;

		extreme1.x = -0.5;
		extreme1.y = -0.5;
		extreme1.z = -0.5;

		extreme2.x = 0.5;
		extreme2.y = 0.5;
		extreme2.z = 0.5;

		float minT;
		float scale;

		minT = std::min(1 / (maxX - minX), 1 / (maxY - minY));
		scale = std::min(minT, 1 / (maxZ - minZ));

		for (int i = 0; i < nv; i++)
		{

			itPoints = points.begin();
			float xt, yt, zt;
			std::advance(itPoints, i);

			//itPoints->x = 2.0*(itPoints->x - minX) / (maxX - minX) - 1.0;
			//itPoints->y = 2.0*(itPoints->y - minY) / (maxY - minY) - 1.0;
			//itPoints->z = 2.0*(itPoints->z - minZ) / (maxZ - minZ) - 1.0;

			itPoints->x = (itPoints->x - 0.5*(minX + maxX))*scale;
			itPoints->y = (itPoints->y - 0.5*(minY + maxY))*scale;
			itPoints->z = (itPoints->z - 0.5*(minZ + maxZ))*scale;

			xt = itPoints->x;
			yt = itPoints->y;
			zt = itPoints->z;

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
			fprintf(fd, "%d %d %d %d\n", type_face_t, point1_t, point2_t, point3_t);
		}

		

	}

}
