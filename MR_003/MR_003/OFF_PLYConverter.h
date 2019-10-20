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
	int numberRegionsShape;
	int singlePass;
	VectorXf minMaxPoints;
	//Number of points per region
	MatrixXi pointsPerRegion;
	//Point and its corresponding region
	MatrixXi pointXRegion;

public:
	OFF_PLYConverter();
	~OFF_PLYConverter();

	void SetNumberRegionsShape(int sp);
	int getNumberRegionsShape();
	void Convert_OFF_PLY(FILE *fo, FILE *fd);
	VectorXi GetRandomIndexes(int first, int sizeSample, int sizePopulation);
	VectorXf GetFeatureVector(VectorXf samples, int numberBins, float minValue, float maxValue);
	float CalculateDiameter();
	float CalculateCompactness();
	float SurfaceArea(MatrixXi* faces, MatrixXf* vertices);
	float DistanceBetweenPoints(Point p1, Point p2);
	VectorXf CalculateHistogram_Bary_RandVert(int sampleSize, int numberBins);
	VectorXf CalculateHistogram_2_RandVert(int sampleSize, int numberBins);
	VectorXf CalculateHistogram_Tetra_4_RandVert(int sampleSize, int numberBins);
	VectorXf GetMinMaxPoints();
	MatrixXi CalculatePointsPerRegion();
	void Process_Post_Norm_PLY(FILE* fo, FILE* fd);
	void Read_PLY_File(FILE* f);

	//VectorXi CalculateHistogram_2_RandVert(int sampleSize, int numberBins);
	



};