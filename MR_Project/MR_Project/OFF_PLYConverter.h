#pragma once

using namespace std;

struct Point
{
	float x;
	float y;
	float z;
};

struct Face
{
	int type_face;
	int point1;
	int point2;
	int point3;
};

class OFF_PLYConverter {

public:
	OFF_PLYConverter();
	~OFF_PLYConverter();

	void Convert_OFF_PLY(FILE *fo, FILE *fd);

};