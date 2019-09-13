#include "UnstructuredGrid.h"




void UnstructuredGrid::getPoint(int i,float* p)
{
	p[0] = pointsX[i];
	p[1] = pointsY[i];
}


void UnstructuredGrid::setPoint(int i,float* p)
{
	pointsX[i] = p[0];
	pointsY[i] = p[1];
}

int	UnstructuredGrid::getCell(int cell,int* vertices)
{
	vertices[0] = cells[3*cell];
	vertices[1] = cells[3*cell+1];
	vertices[2] = cells[3*cell+2];
	
	return 3;
}

void UnstructuredGrid::setCell(int cell,int* vertices)
{
	cells[3*cell]   = vertices[0];
	cells[3*cell+1] = vertices[1];
	cells[3*cell+2] = vertices[2];	
}


int	UnstructuredGrid::findCell(float* p)
{
    return -1;
}
