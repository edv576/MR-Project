#include <string.h>
#include <iostream>
#include <stdlib.h>
#include "OFF_PLYConverter.h"

OFF_PLYConverter::OFF_PLYConverter() {


}

OFF_PLYConverter::~OFF_PLYConverter() {


}

void OFF_PLYConverter::Convert_OFF_PLY(FILE *fo, FILE *fd){

	char buffer[101];

	fscanf(fo, "%100s\n", buffer);
	int nv;
	int nf;
	int ne;
	float x, y, z;
	int tf, p1, p2, p3;

	if (strcmp(buffer, "OFF")) {

		return;
	}
	else
	{
	
		fscanf(fo, "%d %d %d", &nv, &nf, &ne);

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

		for (int i = 0; i < nv; i++)
		{
			fscanf(fo, "%f %f %f", &x, &y, &z);
			fprintf(fd, "%f %f %f\n", x, y, z);
		}

		for (int i = 0; i < nf; i++)
		{
			fscanf(fo, "%d %d %d %d", &tf, &p1, &p2, &p3);
			fprintf(fd, "%d %d %d %d\n", tf, p1, p2, p3);
		}

	}

}
