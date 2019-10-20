#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <list>
#include <math.h>
#include <algorithm>
#include "OFF_PLYConverter.h"
#include <Eigen/Dense>
#include <vector>
#include <random>
#include <chrono>
#include <cassert>

#define _USE_MATH_DEFINES
#include <math.h>


using namespace Eigen;



OFF_PLYConverter::OFF_PLYConverter() {
	numberRegionsShape = 64;
	minMaxPoints.resize(6);
	pointsPerRegion.resize(numberRegionsShape, 2);

	for (int i = 0; i < numberRegionsShape; i++)
	{
		pointsPerRegion(i, 0) = i;
		pointsPerRegion(i, 1) = 0;
	}
}

OFF_PLYConverter::~OFF_PLYConverter() {


}

void OFF_PLYConverter::SetNumberRegionsShape(int sp) {
	numberRegionsShape = sp * sp*sp;
	singlePass = sp;
}

int OFF_PLYConverter::getNumberRegionsShape() {
	return numberRegionsShape;
}

VectorXf OFF_PLYConverter::GetMinMaxPoints() {

	float ex = 1.0e6;


	float minX = ex;
	float minY = ex;
	float minZ = ex;
	float maxX = -ex;
	float maxY = -ex;
	float maxZ = -ex;

	for (int i = 0; i < allPoints.rows(); i++)
	{
		minX = min(allPoints(i, 0), minX);
		minY = min(allPoints(i, 1), minY);
		minZ = min(allPoints(i, 2), minZ);
		maxX = max(allPoints(i, 0), maxX);
		maxY = max(allPoints(i, 1), maxY);
		maxZ = max(allPoints(i, 2), maxZ);
	}

	minMaxPoints(0) = minX;
	minMaxPoints(1) = minY;
	minMaxPoints(2) = minZ;
	minMaxPoints(3) = maxX;
	minMaxPoints(4) = maxY;
	minMaxPoints(5) = maxZ;

	return minMaxPoints;
}

MatrixXi OFF_PLYConverter::CalculatePointsPerRegion() {

	float lengthX;
	float lengthY;
	float lengthZ;

	lengthX = minMaxPoints(3) - minMaxPoints(0);
	lengthY = minMaxPoints(4) - minMaxPoints(1);
	lengthZ = minMaxPoints(5) - minMaxPoints(2);
	bool found = false;

	for (int i = 0; i < allPoints.rows(); i++) {
		found = false;
		for (int j = 0; j < singlePass && !found; j++) {
			for (int k = 0; k < singlePass && !found; k++) {
				for (int l = 0; l < singlePass && !found; l++) {
					if ((allPoints(i, 0) >= minMaxPoints(0) + (lengthX / singlePass)*j) && (allPoints(i, 0) < minMaxPoints(0) + (lengthX / singlePass)*(j + 1))
						&& (allPoints(i, 1) >= minMaxPoints(1) + (lengthY / singlePass)*k) && (allPoints(i, 1) < minMaxPoints(1) + (lengthY / singlePass)*(k + 1))
						&& (allPoints(i, 2) >= minMaxPoints(2) + (lengthZ / singlePass)*l) && (allPoints(i, 2) < minMaxPoints(2) + (lengthZ / singlePass)*(l + 1))){
						pointsPerRegion(l + k + j)++;
						pointXRegion(i, 0) = i;
						pointXRegion(i, 1) = l + k + j;
						found = true;
					}
				}
			}
		}

	}

	return pointsPerRegion;

}


float OFF_PLYConverter::DistanceBetweenPoints(Point p1, Point p2) {

	float distance = sqrtf(pow(p2.x - p1.x, 2) +
		pow(p2.y - p1.y, 2) +
		pow(p2.z - p1.z, 2));

	return distance;

	
}

Point Substract(Point p1, Point p2) {

	Point substraction;

	substraction.x = p1.x - p2.x;
	substraction.y = p1.y - p2.y;
	substraction.z = p1.z - p2.z;

	return substraction;

}

float Determinant_3x3(Matrix3f m) {

	return (m(0, 0) * (m(1, 1)*m(2, 2) - m(1, 2)*m(2, 1)) -
		m(1, 0) * (m(0, 1)*m(2, 2) - m(0, 2)*m(2, 1)) +
		m(2, 0) * (m(0, 1)*m(1, 2) - m(0, 2)*m(1, 1)));

}

void combinationUtil(int arr[], int n, int r,
	int index, int data[], int i)
{
	// Current cobination is ready, print it  
	if (index == r)
	{
		for (int j = 0; j < r; j++)
			cout << data[j] << " ";
		cout << endl;
		return;
	}

	// When no more elements are there to put in data[]  
	if (i >= n)
		return;

	// current is included, put next at next location  
	data[index] = arr[i];
	combinationUtil(arr, n, r, index + 1, data, i + 1);

	// current is excluded, replace it with next (Note that  
	// i+1 is passed, but index is not changed)  
	combinationUtil(arr, n, r, index, data, i + 1);
}

float VolumeOfTetrahedron(Point p1, Point p2, Point p3, Point p4) {

	Matrix3f matForDeterminant;

	Point subs1 = Substract(p1, p2);
	Point subs2 = Substract(p2, p3);
	Point subs3 = Substract(p3, p4);

	matForDeterminant(0, 0) = subs1.x;
	matForDeterminant(0, 1) = subs1.y;
	matForDeterminant(0, 2) = subs1.z;

	matForDeterminant(1, 0) = subs2.x;
	matForDeterminant(1, 1) = subs2.y;
	matForDeterminant(1, 2) = subs2.z;

	matForDeterminant(2, 0) = subs3.x;
	matForDeterminant(2, 1) = subs3.y;
	matForDeterminant(2, 2) = subs3.z;
	
	return(abs(Determinant_3x3(matForDeterminant))/6);

}

float SignedVolumeOfTriangle(Point p1, Point p2, Point p3) {
	float v321 = p3.x*p2.y*p1.z;
	float v231 = p2.x*p3.y*p1.z;
	float v312 = p3.x*p1.y*p2.z;
	float v132 = p1.x*p3.y*p2.z;
	float v213 = p2.x*p1.y*p3.z;
	float v123 = p1.x*p2.y*p3.z;
	return (1.0f / 6.0f)*(-v321 + v231 + v312 - v132 - v213 + v123);
}

float FullVolumeOfMesh(MatrixXi faces, MatrixXf vertices) {

	Point dummyPoint1;
	Point dummyPoint2;
	Point dummyPoint3;

	float fullVolume = 0;

	

	for (int i = 0; i < faces.rows(); i++)
	{
		dummyPoint1.x = vertices(faces(i, 1), 0);
		dummyPoint1.y = vertices(faces(i, 1), 1);
		dummyPoint1.z = vertices(faces(i, 1), 2);

		dummyPoint2.x = vertices(faces(i, 2), 0);
		dummyPoint2.y = vertices(faces(i, 2), 1);
		dummyPoint2.z = vertices(faces(i, 2), 2);

		dummyPoint3.x = vertices(faces(i, 3), 0);
		dummyPoint3.y = vertices(faces(i, 3), 1);
		dummyPoint3.z = vertices(faces(i, 3), 2);

		fullVolume += SignedVolumeOfTriangle(dummyPoint1, dummyPoint2, dummyPoint3);

	}

	return abs(fullVolume);


}

template <typename T>
MatrixXi Combination(const std::vector<T>& v, std::size_t count)
{
	assert(count <= v.size());
	std::vector<bool> bitset(v.size() - count, 0);
	bitset.resize(v.size(), 1);

	do {
		int vertices = 0;
		for (std::size_t i = 0; i != v.size(); ++i) {
			if (bitset[i]) {
				
				
				std::cout << v[i] << " ";
			}
		}
		std::cout << std::endl;
	} while (std::next_permutation(bitset.begin(), bitset.end()));

	MatrixXi response;

	return response;
}




VectorXi OFF_PLYConverter::GetRandomIndexes(int first, int sizeSample, int sizePopulation) {

	VectorXi randomIndexes(sizeSample);
	std::vector<int> numbers;

	for (int i = first; i < sizePopulation; i++)      
		numbers.push_back(i);

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle(numbers.begin(), numbers.end(), std::default_random_engine(seed));

	for (int i = 0; i < sizeSample; i++) {
		randomIndexes(i) = numbers[i];
	}

	return randomIndexes;


}

VectorXf OFF_PLYConverter::GetFeatureVector(VectorXf samples, int numberBins, float minValue, float maxValue) {

	float binInterval = (maxValue - minValue)/numberBins;

	VectorXf featureVector(numberBins);
	featureVector.setZero();

	for (int i = 0; i < samples.size(); i++) {
		for (int j = 0; j < numberBins; j++) {
			if ((samples(i) >= j*binInterval) && (samples(i) < (j + 1)*binInterval)) {
				featureVector(j)++;
			}
		}
	}

	//Normalizing feature vector

	float totalFeatureVector = 0;

	for (int i = 0; i < numberBins; i++) {
		totalFeatureVector += featureVector(i);
	}

	for (int i = 0; i < numberBins; i++) {
		featureVector(i) /= totalFeatureVector;
	}

	return featureVector;


}

VectorXf OFF_PLYConverter::CalculateHistogram_Bary_RandVert(int sampleSize, int numberBins) {

	float maxDistance;
	Point dummyPoint;
	VectorXi randomIndexes;
	VectorXf sampleDistances(2);

	randomIndexes = GetRandomIndexes(0, sampleSize, allPoints.rows());

	dummyPoint.x = allPoints(randomIndexes(0), 0);
	dummyPoint.y = allPoints(randomIndexes(0), 1);
	dummyPoint.z = allPoints(randomIndexes(0), 2);

	maxDistance = DistanceBetweenPoints(centroid, dummyPoint);


	sampleDistances.resize(sampleSize);
	VectorXf featureVector;

	

	for (int i = 0; i < randomIndexes.size(); i++)
	{
		dummyPoint.x = allPoints(randomIndexes(i), 0);
		dummyPoint.y = allPoints(randomIndexes(i), 1);
		dummyPoint.z = allPoints(randomIndexes(i), 2);

		sampleDistances(i) = DistanceBetweenPoints(centroid, dummyPoint);

		if (sampleDistances(i) > maxDistance) {
			maxDistance = sampleDistances(i);
		}
	
	}

	featureVector = GetFeatureVector(sampleDistances, 10, 0, 0.5 * sqrt(3));

	

	return featureVector;
}

VectorXf OFF_PLYConverter::CalculateHistogram_2_RandVert(int sampleSize, int numberBins)
{
	VectorXi randomIndexes;
	Point dummyPoint1;
	Point dummyPoint2;
	int verticesXsample = 5;
	//VectorXf sampleDistances(sampleSize*verticesXsample);
	VectorXf sampleDistances((allPoints.rows() - 1) * sampleSize);
	VectorXf featureVector;
	

	randomIndexes = GetRandomIndexes(0, sampleSize, allPoints.rows());
	int actualSizeSamples = 0;
	
	int range = sampleSize - 1;

	for (int i = 0; i < randomIndexes.rows(); i++) {
		dummyPoint1.x = allPoints(randomIndexes(i), 0);
		dummyPoint1.y = allPoints(randomIndexes(i), 1);
		dummyPoint1.z = allPoints(randomIndexes(i), 2);

		//int j = 0;
		//int randomIndex2;


		//while (j < verticesXsample) {
		//	randomIndex2 = rand() % (allPoints.rows());			
		//	if (randomIndexes(i) != randomIndex2) {
		//		dummyPoint2.x = allPoints(randomIndex2, 0);
		//		dummyPoint2.y = allPoints(randomIndex2, 1);
		//		dummyPoint2.z = allPoints(randomIndex2, 2);
		//		sampleDistances(actualSizeSamples) = DistanceBetweenPoints(dummyPoint1, dummyPoint2);
		//		actualSizeSamples++;
		//		j++;
		//	}

		//}

		for (int j = 0; j < allPoints.rows(); j++)
		{
			if (randomIndexes(i) != j) {
				dummyPoint2.x = allPoints(j, 0);
				dummyPoint2.y = allPoints(j, 1);
				dummyPoint2.z = allPoints(j, 2);
				sampleDistances(actualSizeSamples) = DistanceBetweenPoints(dummyPoint1, dummyPoint2);
				actualSizeSamples++;
			}
		}
	}

	sampleDistances.resize(actualSizeSamples);

	featureVector = GetFeatureVector(sampleDistances, 10, 0, 1*sqrt(3));

	return featureVector;

}

VectorXf OFF_PLYConverter::CalculateHistogram_Tetra_4_RandVert(int sampleSize, int numberBins) {

	VectorXi randomIndexes;
	
	VectorXf featureVector;

	randomIndexes = GetRandomIndexes(0, sampleSize*4, allPoints.rows());
	//int verticesXsample = 50;
	int verticesXsample = allPoints.rows();
	VectorXf sampleVolumes(sampleSize*verticesXsample*verticesXsample*verticesXsample);
	Point dummyPoint1, dummyPoint2, dummyPoint3, dummyPoint4;
	int actualVolumeSamples = 0;
	int dummyIndex1, dummyIndex2, dummyIndex3;

	for (int i = 0; i < sampleSize; i++) {
		dummyPoint1.x = allPoints(randomIndexes(i), 0);
		dummyPoint1.y = allPoints(randomIndexes(i), 1);
		dummyPoint1.z = allPoints(randomIndexes(i), 2);
		
		/*for (int j = 0; j < allPoints.rows(); j++)
		{
			if (randomIndexes(i) != j) {
				dummyPoint2.x = allPoints(j, 0);
				dummyPoint2.y = allPoints(j, 1);
				dummyPoint2.z = allPoints(j, 2);

				for (int k = 0; k < allPoints.rows(); k++)
				{
					if (randomIndexes(i) != k && j != k) {
						dummyPoint3.x = allPoints(k, 0);
						dummyPoint3.y = allPoints(k, 1);
						dummyPoint3.z = allPoints(k, 2);

						for (int k = 0; k < allPoints.rows(); k++)
						{
							if (randomIndexes(i) != k && j != k) {
								dummyPoint3.x = allPoints(k, 0);
								dummyPoint3.y = allPoints(k, 1);
								dummyPoint3.z = allPoints(k, 2);



							}
						}

					}
				}


				
			}
		}*/





		//sampleVolumes(actualVolumeSamples) = DistanceBetweenPoints(dummyPoint1, dummyPoint2);
		//actualVolumeSamples++;

		VectorXi randomIndexes2 = GetRandomIndexes(sampleSize, verticesXsample, sampleSize * 2);

		for (int j = 0; j < verticesXsample; j++) {


			dummyPoint2.x = allPoints(randomIndexes2(j), 0);
			dummyPoint2.y = allPoints(randomIndexes2(j), 1);
			dummyPoint2.z = allPoints(randomIndexes2(j), 2);

			VectorXi randomIndexes3 = GetRandomIndexes(sampleSize*2, verticesXsample, sampleSize * 3);

			for (int k = 0; k < verticesXsample; k++) {
				dummyPoint3.x = allPoints(randomIndexes3(k), 0);
				dummyPoint3.y = allPoints(randomIndexes3(k), 1);
				dummyPoint3.z = allPoints(randomIndexes3(k), 2);

				VectorXi randomIndexes4 = GetRandomIndexes(sampleSize * 3, verticesXsample, sampleSize * 4);

				for (int l = 0; l < verticesXsample; l++) {
					dummyPoint4.x = allPoints(randomIndexes4(l), 0);
					dummyPoint4.y = allPoints(randomIndexes4(l), 1);
					dummyPoint4.z = allPoints(randomIndexes4(l), 2);

					sampleVolumes(actualVolumeSamples) = VolumeOfTetrahedron(dummyPoint1, dummyPoint2, dummyPoint3, dummyPoint4);
					actualVolumeSamples++;
				}
			}
		}

	}

	sampleVolumes.resize(actualVolumeSamples);

	Point p1, p2, p3, p4;

	GetMinMaxPoints();

	p1.x = minMaxPoints(0);
	p1.y = minMaxPoints(1);
	p1.z = minMaxPoints(2);

	p2.x = minMaxPoints(3);
	p2.y = minMaxPoints(4);
	p2.z = minMaxPoints(5);

	p3.x = p2.x;
	p3.y = p1.y;
	p3.z = p1.z;

	p4.x = p1.x;
	p4.y = p2.y;
	p4.z = p1.z;


	featureVector = GetFeatureVector(sampleVolumes, 10, 0, VolumeOfTetrahedron(p1, p2, p3, p4));

	return featureVector;

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


float OFF_PLYConverter::CalculateCompactness()
{
	float volume = FullVolumeOfMesh(allFaces, allPoints);
	float surfaceArea = SurfaceArea(&allFaces, &allPoints);


	return powf(surfaceArea, 3) / (36 * M_PI * powf(volume, 2));

}

float OFF_PLYConverter::SurfaceArea(MatrixXi* faces, MatrixXf* vertices) {

	//total area
	float area = 0.0f;

	for (int i = 0; i < faces->rows(); i++)
	{
		Vector4i face = faces->row(i);

		//get the vertices of the face (triangle)
		Vector3f a = vertices->row(face[1]);
		Vector3f b = vertices->row(face[2]);
		Vector3f c = vertices->row(face[3]);

		//https://math.stackexchange.com/questions/128991/how-to-calculate-the-area-of-a-3d-triangle
		//area of a triangle = cross(ab,ac).norm() / 2
		area += (a - b).cross(a - c).norm() / 2;
	}

	return area;
}

void OFF_PLYConverter::Read_PLY_File(FILE* f) {

	char buffer[101];
	char buffer1[101], buffer2[101];

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
	//Point centroid;
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

	fscanf(f, "%100s\n", buffer);

	if (strcmp(buffer, "ply")) {

		return;
	}
	else
	{
		fgets(buffer, 100, f);
		fgets(buffer, 100, f);
		fscanf(f, "%100s %100s %d\n", buffer1, buffer2, &nv);
		fgets(buffer, 100, f);
		fgets(buffer, 100, f);
		fgets(buffer, 100, f);
		fscanf(f, "%100s %100s %d\n", buffer1, buffer2, &nf);
		fgets(buffer, 100, f);
		fgets(buffer, 100, f);

		allPoints.resize(nv, 3);
		allFaces.resize(nf, 4);

		for (int i = 0; i < nv; i++)
		{
			Point point;
			fscanf(f, "%f %f %f", &x, &y, &z);

			allPoints(i, 0) = x;
			allPoints(i, 1) = y;
			allPoints(i, 2) = z;

			centroid.x += x;
			centroid.y += y;
			centroid.z += z;
		}

		for (int i = 0; i < nf; i++)
		{
			Face face;
			fscanf(f, "%d %d %d %d", &tf, &p1, &p2, &p3);
			allFaces(i, 0) = tf;
			allFaces(i, 1) = p1;
			allFaces(i, 2) = p2;
			allFaces(i, 3) = p3;

		}

		//Calculating the actual Centroid
		centroid.x /= nv;
		centroid.y /= nv;
		centroid.z /= nv;
	}
}

void OFF_PLYConverter::Process_Post_Norm_PLY(FILE* fo, FILE* fd)
{
	char buffer[101];
	char buffer1[101], buffer2[101];

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
	//Point centroid;
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

	fscanf(fo, "%100s\n", buffer);

	if (strcmp(buffer, "ply")) {

		return;
	}
	else
	{
		fgets(buffer, 100, fo);
		fgets(buffer, 100, fo);
		fscanf(fo, "%100s %100s %d\n", buffer1, buffer2, &nv);
		fgets(buffer, 100, fo);
		fgets(buffer, 100, fo);
		fgets(buffer, 100, fo);
		fscanf(fo, "%100s %100s %d\n", buffer1, buffer2, &nf);
		fgets(buffer, 100, fo);
		fgets(buffer, 100, fo);

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

		for (int i = 0; i < nf; i++)
		{
			Face face;
			fscanf(fo, "%d %d %d %d", &tf, &p1, &p2, &p3);
			allFaces(i, 0) = tf;
			allFaces(i, 1) = p1;
			allFaces(i, 2) = p2;
			allFaces(i, 3) = p3;

		}

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

		float ex = 1.0e6;

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

		centroid.x = 0;
		centroid.y = 0;
		centroid.z = 0;

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
			dummySum1 += (allPoints(i, 0) - meanX) * (allPoints(i, 0) - meanX);
			dummySum2 += (allPoints(i, 1) - meanY) * (allPoints(i, 1) - meanY);
			dummySum3 += (allPoints(i, 2) - meanZ) * (allPoints(i, 2) - meanZ);

		}

		varX = dummySum1 / (nv - 1);
		varY = dummySum2 / (nv - 1);
		varZ = dummySum3 / (nv - 1);

		dummySum1 = 0;
		dummySum2 = 0;
		dummySum3 = 0;

		for (int i = 0; i < nv; i++)
		{
			dummySum1 += (allPoints(i, 0) - meanX) * (allPoints(i, 1) - meanY);
			dummySum2 += (allPoints(i, 0) - meanX) * (allPoints(i, 2) - meanZ);
			dummySum3 += (allPoints(i, 1) - meanY) * (allPoints(i, 2) - meanZ);
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
			for (int j = i + 1; j < 3; j++)
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
			triangleCenter(0) = (allPoints(allFaces(i, 1), 0) + allPoints(allFaces(i, 2), 0) + allPoints(allFaces(i, 3), 0)) / 3;
			triangleCenter(1) = (allPoints(allFaces(i, 1), 1) + allPoints(allFaces(i, 2), 1) + allPoints(allFaces(i, 3), 1)) / 3;
			triangleCenter(2) = (allPoints(allFaces(i, 1), 2) + allPoints(allFaces(i, 2), 2) + allPoints(allFaces(i, 3), 2)) / 3;

			fX += (triangleCenter(0) / abs(triangleCenter(0))) * pow(triangleCenter(0), 2);
			fY += (triangleCenter(1) / abs(triangleCenter(1))) * pow(triangleCenter(1), 2);
			fZ += (triangleCenter(2) / abs(triangleCenter(2))) * pow(triangleCenter(2), 2);

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



		Vector3f A, B, A0, B0;

		A(0) = -0.5;
		A(1) = -0.5;
		A(2) = -0.5;

		B(0) = 0.5;
		B(1) = 0.5;
		B(2) = 0.5;

		A0(0) = minX;
		A0(1) = minY;
		A0(2) = minZ;

		B0(0) = maxX;
		B0(1) = maxY;
		B0(2) = maxZ;

		float scale = min((B(0) - A(0)) / (B0(0) - A0(0)), min((B(1) - A(1)) / (B0(1) - A0(1)), (B(2) - A(2)) / (B0(2) - A0(2))));



		for (int i = 0; i < nv; i++)
		{

			//itPoints = points.begin();
			float xt, yt, zt;


			xt = allPoints(i, 0);
			yt = allPoints(i, 1);
			zt = allPoints(i, 2);

			allPoints(i, 0) = (allPoints(i, 0) - 0.5 * (A0(0) + B0(0))) * scale + 0.5 * (A(0) + B(0));
			allPoints(i, 1) = (allPoints(i, 1) - 0.5 * (A0(1) + B0(1))) * scale + 0.5 * (A(1) + B(1));
			allPoints(i, 2) = (allPoints(i, 2) - 0.5 * (A0(2) + B0(2))) * scale + 0.5 * (A(2) + B(2));

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

	}

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
	//Point centroid;
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

		centroid.x = 0;
		centroid.y = 0;
		centroid.z = 0;

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

		////Normalizing to size 1

		//float sizeX = maxX - minX;								//2. Compute a single scaling factor that best fits the grid
		//sizeX = (sizeX) ? 1 / sizeX : 0.5;							//   in the [-1,1] cube. Using a single factor for x,y, and z
		//float sizeY = maxY - minY;								//   ensures that the object is scaled while keeping its
		//sizeY = (sizeY) ? 1 / sizeY : 0.5;							//   aspect ratio.
		//float sizeZ = maxZ - minZ;
		//sizeZ = (sizeZ) ? 1 / sizeZ : 0.5;

		//float scale = min(sizeX, min(sizeY, sizeZ));

		////Initializing the extreme points for normalization. Bounding box of length 1.


		//for (int i = 0; i < nv; i++)
		//{

		//	//itPoints = points.begin();
		//	float xt, yt, zt;


		//	xt = allPoints(i, 0);
		//	yt = allPoints(i, 1);
		//	zt = allPoints(i, 2);

		//	allPoints(i, 0) = 2 * ((allPoints(i, 0) - minX)*scale - 0.5);
		//	allPoints(i, 1) = 2 * ((allPoints(i, 1) - minY)*scale - 0.5);
		//	allPoints(i, 2) = 2 * ((allPoints(i, 2) - minZ)*scale - 0.5);


		//	xt = allPoints(i, 0);
		//	yt = allPoints(i, 1);
		//	zt = allPoints(i, 2);

		//	//Writing the vertices values in the .ply file
		//	fprintf(fd, "%f %f %f\n", xt, yt, zt);

		//}

		Vector3f A, B, A0, B0;

		A(0) = -0.5;
		A(1) = -0.5;
		A(2) = -0.5;

		B(0) = 0.5;
		B(1) = 0.5;
		B(2) = 0.5;

		A0(0) = minX;
		A0(1) = minY;
		A0(2) = minZ;

		B0(0) = maxX;
		B0(1) = maxY;
		B0(2) = maxZ;

		float scale = min((B(0) - A(0)) / (B0(0) - A0(0)), min((B(1) - A(1)) / (B0(1) - A0(1)), (B(2) - A(2)) / (B0(2) - A0(2))));



		for (int i = 0; i < nv; i++)
		{

			//itPoints = points.begin();
			float xt, yt, zt;


			xt = allPoints(i, 0);
			yt = allPoints(i, 1);
			zt = allPoints(i, 2);

			allPoints(i, 0) = (allPoints(i, 0) - 0.5*(A0(0) + B0(0)))*scale + 0.5*(A(0) + B(0));
			allPoints(i, 1) = (allPoints(i, 1) - 0.5*(A0(1) + B0(1)))*scale + 0.5*(A(1) + B(1));
			allPoints(i, 2) = (allPoints(i, 2) - 0.5*(A0(2) + B0(2)))*scale + 0.5*(A(2) + B(2));

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

		pointXRegion.resize(allPoints.rows(), 2);

		//VectorXf hist_Bary_RandVert;
		//hist_Bary_RandVert = CalculateHistogram_Bary_RandVert(200, 10);

		//VectorXf hist_2_RandVert = CalculateHistogram_2_RandVert(200, 10);

		// 

		////CalculateHistogram_Bary_RandVert(200, 10);

		//VectorXf hist_Tetra_4_RandVert = CalculateHistogram_Tetra_4_RandVert(50, 10);

		float maxDistance = CalculateDiameter();
		float compactness = CalculateCompactness();

		int t = 0;

		Point p1, p2, p3, p4;

		p1.x = -0.5;
		p1.y = -0.5;
		p1.z = 0.5;

		p2.x = -0.5;
		p2.y = 0.5;
		p2.z = 0.5;

		p3.x = 0.5;
		p3.y = 0.5;
		p3.z = 0.5;

		p4.x = -0.5;
		p4.y = 0.5;
		p4.z = -0.5;

		float volumeTest = VolumeOfTetrahedron(p1, p2, p3, p4);
	

	}

	


	
	

}


