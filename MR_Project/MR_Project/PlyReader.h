#pragma once



class UnstructuredGrid3D;


class PlyReader
{
public:

					PlyReader() { }
UnstructuredGrid3D*	read(const char* filename);
};


			

