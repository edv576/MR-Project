#pragma once

#include "UnstructuredGrid.h"
#include "VectorAttributes.h"

using namespace std;




class UnstructuredGrid3D : public UnstructuredGrid
{
public: 
					UnstructuredGrid3D(int P,int C): UnstructuredGrid(P,C),pointNormals(P),faceNormals(C), pointsZ(P)
					{ }	 
					 	
void				getPoint(int i,float* p);

void				setPoint(int i,float* p);		

void				normalize();

void				computeFaceNormals();

void				computeVertexNormals();

VectorAttributes&	getFaceNormals()
					{ return faceNormals; }

VectorAttributes&	getPointNormals()
					{ return pointNormals; }

private:

vector<float>		pointsZ;
VectorAttributes    pointNormals;
VectorAttributes    faceNormals;
};



