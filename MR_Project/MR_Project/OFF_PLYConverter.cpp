#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <list>
#include <math.h>
#include <algorithm>
#include "OFF_PLYConverter.h"
#include <Eigen/Dense>


using namespace Eigen;

OFF_PLYConverter::OFF_PLYConverter() {


}

OFF_PLYConverter::~OFF_PLYConverter() {


}


void OFF_PLYConverter::Convert_OFF_PLY(FILE *fo, FILE *fd){

	char buffer[101];
	std::list<Point> points;
	std::list<Face> faces;
	MatrixXf allPoints(2, 2);



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

		allPoints.resize(nv, 3);

		for (int i = 0; i < nv; i++)
		{
			Point point;
			fscanf(fo, "%f %f %f", &x, &y, &z);

			allPoints(i, 0) = x;
			allPoints(i, 1) = y;
			allPoints(i, 2) = z;

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

		//Applying PCA

		Matrix3f covarianceMatrixXYZ;

		float varX, varY, varZ;
		float covXY, covXZ, covYZ;
		float meanX, meanY, meanZ;
		float dummySum1 = 0;
		float dummySum2 = 0;
		float dummySum3 = 0;

		meanX = 0;
		meanY = 0;
		meanZ = 0;

		for (int i = 0; i < nv; i++)
		{
			itPoints = points.begin();
			std::advance(itPoints, i);

			meanX += itPoints->x;
			meanY += itPoints->y;
			meanZ += itPoints->z;

		}

		meanX /= nv;
		meanY /= nv;
		meanZ /= nv;

		for (int i = 0; i < nv; i++)
		{
			itPoints = points.begin();
			std::advance(itPoints, i);
			dummySum1 += (itPoints->x - meanX)*(itPoints->x - meanX);
			dummySum2 += (itPoints->y - meanY)*(itPoints->y - meanY);
			dummySum3 += (itPoints->z - meanZ)*(itPoints->z - meanZ);

		}

		varX = dummySum1 / (nv - 1);
		varY = dummySum2 / (nv - 1);
		varZ = dummySum3 / (nv - 1);

		dummySum1 = 0;
		dummySum2 = 0;
		dummySum3 = 0;

		for (int i = 0; i < nv; i++)
		{
			itPoints = points.begin();
			std::advance(itPoints, i);
			dummySum1 += (itPoints->x - meanX)*(itPoints->y - meanY);
			dummySum2 += (itPoints->x - meanX)*(itPoints->z - meanZ);
			dummySum3 += (itPoints->y - meanY)*(itPoints->z - meanZ);
		}

		covXY = dummySum1 / (nv - 1);
		covXZ = dummySum2 / (nv - 1);
		covYZ = dummySum3 / (nv - 1);

		covarianceMatrixXYZ(0, 0) = varX;
		covarianceMatrixXYZ(1, 1) = varY;
		covarianceMatrixXYZ(2, 2) = varZ;
		covarianceMatrixXYZ(0, 1) = covXY;
		covarianceMatrixXYZ(1, 0) = covXY;
		covarianceMatrixXYZ(0, 2) = covXZ;
		covarianceMatrixXYZ(2, 0) = covXZ;
		covarianceMatrixXYZ(1, 2) = covYZ;
		covarianceMatrixXYZ(2, 1) = covYZ;

		EigenSolver<Matrix3f> solver(covarianceMatrixXYZ);

		Vector3cf eigenValues = solver.eigenvalues();
		Matrix3cf eigenVectors = solver.eigenvectors();

		Vector3f eigen1, eigen2, eigen3;



		eigen1(0) = eigenVectors(0, 0).real();
		eigen1(1) = eigenVectors(0, 1).real();
		eigen1(2) = eigenVectors(0, 2).real();
		eigen2(0) = eigenVectors(1, 0).real();
		eigen2(1) = eigenVectors(1, 1).real();
		eigen2(2) = eigenVectors(1, 2).real();
		eigen3(0) = eigenVectors(2, 0).real();
		eigen3(1) = eigenVectors(2, 1).real();
		eigen3(2) = eigenVectors(2, 2).real();

		Vector3i indexes;
		Vector3f eValues;
		int indexTemp;
		float eValueTemp;

		indexes(0) = 0;
		indexes(1) = 1;
		indexes(2) = 2;
		eValues(0) = eigenValues(0).real();
		eValues(1) = eigenValues(1).real();
		eValues(2) = eigenValues(2).real();

		for (int i = 0; i < 2; i++)
		{
			for (int j = i+1; j < 3; j++)
			{
				if (eValues(j) > eValues(i)) {
					eValueTemp = eValues(i);
					eValues(i) = eValues(j);
					eValues(j) = eValueTemp;
					indexTemp = indexes(i);
					indexes(i) = indexes(j);
					indexes(j) = indexTemp;
				}
			}
		}

		Vector3f xAxis, yAxis, zAxis;

		xAxis(0) = eigenVectors(indexes(0), 0).real();
		xAxis(1) = eigenVectors(indexes(0), 1).real();
		xAxis(2) = eigenVectors(indexes(0), 2).real();

		yAxis(0) = eigenVectors(indexes(1), 0).real();
		yAxis(1) = eigenVectors(indexes(1), 1).real();
		yAxis(2) = eigenVectors(indexes(1), 2).real();

		zAxis(0) = eigenVectors(indexes(2), 0).real();
		zAxis(1) = eigenVectors(indexes(2), 1).real();
		zAxis(2) = eigenVectors(indexes(2), 2).real();

		for (int i = 0; i < nv; i++)
		{
			itPoints = points.begin();
			std::advance(itPoints, i);

			Vector3f tempVec;
			tempVec(0) = itPoints->x;
			tempVec(1) = itPoints->y;
			tempVec(2) = itPoints->z;

			itPoints->x = xAxis.dot(tempVec);
			itPoints->y = yAxis.dot(tempVec);
			itPoints->z = zAxis.dot(tempVec);
		}


		//Normalizing to size 1

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
