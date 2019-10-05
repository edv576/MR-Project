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

float OFF_PLYConverter::DistanceBetweenPoints(Point p1, Point p2) {

	float distance = sqrtf(pow(p2.x - p1.x, 2) +
		pow(p2.y - p1.y, 2) +
		pow(p2.z - p1.z, 2));

	return distance;
}

float OFF_PLYConverter::CalculateDiameter()
{
	Point dummyPoint1;
	Point dummyPoint2;

	dummyPoint1.x = allPoints(0, 0);
	dummyPoint1.y = allPoints(0, 1);
	dummyPoint1.z = allPoints(0, 2);

	dummyPoint2.x = allPoints(1, 0);
	dummyPoint2.y = allPoints(1, 1);
	dummyPoint2.z = allPoints(1, 2);

	float maxDistance = DistanceBetweenPoints(dummyPoint1, dummyPoint2);
	float tempDistance;

	for (int i = 0; i < allPoints.rows() - 1; i++)
	{
		for (int j = i+1; j < allPoints.rows(); j++)
		{
			dummyPoint1.x = allPoints(i, 0);
			dummyPoint1.y = allPoints(i, 1);
			dummyPoint1.z = allPoints(i, 2);

			dummyPoint2.x = allPoints(j, 0);
			dummyPoint2.y = allPoints(j, 1);
			dummyPoint2.z = allPoints(j, 2);

			tempDistance = DistanceBetweenPoints(dummyPoint1, dummyPoint2);

			if (tempDistance > maxDistance) {
				maxDistance = tempDistance;
			}
		}
	}
	return maxDistance;
}





void OFF_PLYConverter::Convert_OFF_PLY(FILE *fo, FILE *fd){

	char buffer[101];
	/*std::list<Point> points;
	std::list<Face> faces;*/
	//MatrixXf allPoints(2, 2);
	//MatrixXi allFaces(2, 2);



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
		allFaces.resize(nf, 4);

		for (int i = 0; i < nv; i++)
		{
			Point point;
			fscanf(fo, "%f %f %f", &x, &y, &z);

			allPoints(i, 0) = x;
			allPoints(i, 1) = y;
			allPoints(i, 2) = z;

			centroid.x += x;
			centroid.y += y;
			centroid.z += z;			
		}

		//std::list<Point>::iterator itPoints = points.begin();

		for (int i = 0; i < nf; i++)
		{
			Face face;
			fscanf(fo, "%d %d %d %d", &tf, &p1, &p2, &p3);
			allFaces(i, 0) = tf;
			allFaces(i, 1) = p1;
			allFaces(i, 2) = p2;
			allFaces(i, 3) = p3;
	
		}

		//std::list<Face>::iterator itFaces = faces.begin();

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


		//itPoints = points.begin();

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

			allPoints(i, 0) += correction.x;
			allPoints(i, 1) += correction.y;
			allPoints(i, 2) += correction.z;
			
			minX = min(allPoints(i, 0), minX);
			minY = min(allPoints(i, 1), minY);
			minZ = min(allPoints(i, 2), minZ);
			maxX = max(allPoints(i, 0), maxX);
			maxY = max(allPoints(i, 1), maxY);
			maxZ = max(allPoints(i, 2), maxZ);

			
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
			meanX += allPoints(i, 0);
			meanY += allPoints(i, 1);
			meanZ += allPoints(i, 2);


		}

		meanX /= nv;
		meanY /= nv;
		meanZ /= nv;

		for (int i = 0; i < nv; i++)
		{
			dummySum1 += (allPoints(i, 0) - meanX)*(allPoints(i, 0) - meanX);
			dummySum2 += (allPoints(i, 1) - meanY)*(allPoints(i, 1) - meanY);
			dummySum3 += (allPoints(i, 2) - meanZ)*(allPoints(i, 2) - meanZ);

		}

		varX = dummySum1 / (nv - 1);
		varY = dummySum2 / (nv - 1);
		varZ = dummySum3 / (nv - 1);

		dummySum1 = 0;
		dummySum2 = 0;
		dummySum3 = 0;

		for (int i = 0; i < nv; i++)
		{
			dummySum1 += (allPoints(i, 0) - meanX)*(allPoints(i, 1) - meanY);
			dummySum2 += (allPoints(i, 0) - meanX)*(allPoints(i, 2) - meanZ);
			dummySum3 += (allPoints(i, 1) - meanY)*(allPoints(i, 2) - meanZ);
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
			Vector3f tempVec;

			tempVec(0) = allPoints(i, 0);
			tempVec(1) = allPoints(i, 1);
			tempVec(2) = allPoints(i, 2);

			allPoints(i, 0) = xAxis.dot(tempVec);
			allPoints(i, 1) = yAxis.dot(tempVec);
			allPoints(i, 2) = zAxis.dot(tempVec);
		}




		float fX, fY, fZ;
		fX = 0;
		fY = 0;
		fZ = 0;

		for (int i = 0; i < nf; i++)
		{
			Vector3f triangleCenter;
			triangleCenter(0) = (allPoints(allFaces(i, 1), 0) + allPoints(allFaces(i, 2), 0) + allPoints(allFaces(i, 3), 0))/3;
			triangleCenter(1) = (allPoints(allFaces(i, 1), 1) + allPoints(allFaces(i, 2), 1) + allPoints(allFaces(i, 3), 1))/3;
			triangleCenter(2) = (allPoints(allFaces(i, 1), 2) + allPoints(allFaces(i, 2), 2) + allPoints(allFaces(i, 3), 2))/3;

			fX += (triangleCenter(0) / abs(triangleCenter(0)))*pow(triangleCenter(0), 2);
			fY += (triangleCenter(1) / abs(triangleCenter(1)))*pow(triangleCenter(1), 2);
			fZ += (triangleCenter(2) / abs(triangleCenter(2)))*pow(triangleCenter(2), 2);

		}

		float sign_fX, sign_fY, sign_FZ;
		sign_fX = fX / abs(fX);
		sign_fY = fY / abs(fY);
		sign_FZ = fZ / abs(fZ);

		for (int i = 0; i < nv; i++)
		{
			allPoints(i, 0) *= sign_fX;
			allPoints(i, 1) *= sign_fY;
			allPoints(i, 2) *= sign_FZ;
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

			//itPoints = points.begin();
			float xt, yt, zt;
			//std::advance(itPoints, i);

			////itPoints->x = 2.0*(itPoints->x - minX) / (maxX - minX) - 1.0;
			////itPoints->y = 2.0*(itPoints->y - minY) / (maxY - minY) - 1.0;
			////itPoints->z = 2.0*(itPoints->z - minZ) / (maxZ - minZ) - 1.0;

			////Doing the scaling of the mesh. Its done keeping the aspect ratio
			//itPoints->x = 2 * ((itPoints->x - minX)*scale - 0.5);
			//itPoints->y = 2 * ((itPoints->y - minY)*scale - 0.5);
			//itPoints->z = 2 * ((itPoints->z - minZ)*scale - 0.5);
			////itPoints->x = (itPoints->x - 0.5*(minX + maxX))*scale;
			////itPoints->y = (itPoints->y - 0.5*(minY + maxY))*scale;
			////itPoints->z = (itPoints->z - 0.5*(minZ + maxZ))*scale;

			xt = allPoints(i, 0);
			yt = allPoints(i, 1);
			zt = allPoints(i, 2);

			allPoints(i, 0) = 2 * ((allPoints(i, 0) - minX)*scale - 0.5);
			allPoints(i, 1) = 2 * ((allPoints(i, 1) - minY)*scale - 0.5);
			allPoints(i, 2) = 2 * ((allPoints(i, 2) - minZ)*scale - 0.5);
			//itPoints->x = (itPoints->x - 0.5*(minX + maxX))*scale;
			//itPoints->y = (itPoints->y - 0.5*(minY + maxY))*scale;
			//itPoints->z = (itPoints->z - 0.5*(minZ + maxZ))*scale;

			xt = allPoints(i, 0);
			yt = allPoints(i, 1);
			zt = allPoints(i, 2);

			//Writing the vertices values in the .ply file
			fprintf(fd, "%f %f %f\n", xt, yt, zt);

		}

		for (int i = 0; i < nf; i++)
		{
			int type_face_t, point1_t, point2_t, point3_t;

			type_face_t = allFaces(i, 0);
			point1_t = allFaces(i, 1);
			point2_t = allFaces(i, 2);
			point3_t = allFaces(i, 3);
			//Writing the faces values in the .ply file
			fprintf(fd, "%d %d %d %d\n", type_face_t, point1_t, point2_t, point3_t);
		}

		float maxDistance = CalculateDiameter();

		int t = 0;

	}

	


	
	

}


