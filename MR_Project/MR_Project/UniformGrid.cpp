#include "UniformGrid.h"
#include <math.h>


UniformGrid::~UniformGrid()
{  }


int	UniformGrid::getCell(int i,int* v)
{  
	int cell_row = i / (N1-1);
	int cell_col = i % (N1-1);
	
    v[0] = i + cell_row;
	v[1] = v[0]+1;
	v[2] = v[1]+N1;
	v[3] = v[0]+N1;
	
	return 4;
}





int UniformGrid::findCell(float* p)
{
	int C[2];

	//compute	c e l l	coordinates	C[ 0 ] ,C[ 1 ]
	C[0] = floor((p[0]-m1)*N1/d1);
	C[1] = floor((p[1]-m2)*N2/d2);


	//test if p is inside the dataset
	if (C[0]<0 || C[0]>=N1-1 || C[1]<0 || C[1]>=N2-1)
	  return -1;

	//go from cell coordinates to cell index
	return C[0] + C[1]*N1;
}


	