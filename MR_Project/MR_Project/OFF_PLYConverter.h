#pragma once

using namespace std;

class OFF_PLYConverter {

public:
	OFF_PLYConverter();
	~OFF_PLYConverter();

	void Convert_OFF_PLY(FILE *fo, FILE *fd);

};