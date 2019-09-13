#include "RectilinearGrid.h"
#include <math.h>




void RectilinearGrid::getPoint(int i, float* p)
{
	//Compute vertex grid-coordinates row, col

	int row = i/N1;
	int col = i%N1;	

	//Compute actual vertex coordinates using the sampling-steps stored in dX[],dY[] for the X and Y axes
	p[0] = dX[col];  
	p[1] = dY[row];
}

