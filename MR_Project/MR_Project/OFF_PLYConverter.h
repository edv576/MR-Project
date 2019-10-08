#pragma once
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//Struct representing a point in 3D
struct Point
{
	float x;
	float y;
	float z;
};

//Struct containing the values for a face
struct Face
{
	int type_face;
	int point1;
	int point2;
	int point3;
};

class OFF_PLYConverter {

private:
	MatrixXf allPoints;
	MatrixXi allFaces;
	Point centroid;

public:
	OFF_PLYConverter();
	~OFF_PLYConverter();

	void Convert_OFF_PLY(FILE *fo, FILE *fd);
	VectorXi GetRandomIndexes(int size);
	VectorXi GetFeatureVector(VectorXf samples, int numberBins, float minValue, float maxValue);
	float CalculateDiameter();
	float CalculateCompactness(MatrixXi faces, MatrixXf vertices);
	float SurfaceArea(MatrixXi* faces, MatrixXf* vertices);
	float DistanceBetweenPoints(Point p1, Point p2);
	VectorXi CalculateHistogram_Bary_RandVert(int sampleSize, int numberBins);
	//VectorXi CalculateHistogram_2_RandVert(int sampleSize, int numberBins);
	



};