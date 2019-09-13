#include "StructuredGrid.h"


void StructuredGrid::getPoint(int i,float* p)
{
	p[0] = pointsX[i];
	p[1] = pointsY[i];
}


void StructuredGrid::setPoint(int i,float* p)
{
	pointsX[i] = p[0];
	pointsY[i] = p[1];
}

			
int	StructuredGrid::findCell(float* p)
{
    return -1;
}







