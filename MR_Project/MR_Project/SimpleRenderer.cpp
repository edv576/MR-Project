#include "SimpleRenderer.h"
#include "Grid.h"
#include <GL/glut.h>									//GLUT library



void SimpleRenderer::draw(Grid& g)
{
	glClearColor(1,1,1,1);								//1. Clear the frame and depth buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glColor3f(0,0,0);									//2. Use black to render cells
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);			//3. Render cells as outlines
	

	for(int i=0;i<g.numCells();++i)						//4. Draw all cells in the grid as black outlines
	{
		int   cell[10];									//Should have enough storage space for our largest cell
		float x[2];
		
		int nvertices = g.getCell(i,cell);				//4.1. Get the 'nvertices' vertex-IDs of i-th cell
		if (nvertices!=3 && nvertices!=4)				//     We only handle here drawing of triangle and quad cells.
		   continue;									//     It is quite simple to extend this to other cell types.

		if (nvertices==3)								//     Triangle cells:
			glBegin(GL_TRIANGLES);						
		else //nvertices==4	
			glBegin(GL_QUADS);							//     Quad cells:
		
		for(int j=0;j<nvertices;++j)					//4.2. Render current cell
		{
			g.getPoint(cell[j],x);						//     Get vertex coordinates of j-th vertex of i-th cell		
			glVertex3f(x[0],x[1],0);					//     Pass this coordinate to OpenGL
		}	
		glEnd();                                        //4.3. Finish rendering current cell
	}
	
	glPointSize(3);										//5.   Draw all vertices in the grid as large red points
	glColor3f(1,0,0);   
	glBegin(GL_POINTS);									//5.1. Draw vertices as an OpenGL point-set	
	for(int i=0;i<g.numPoints();++i)						
	{
		float v[2];
		g.getPoint(i,v);								//5.2. Get coordinates of i-th vertex	
		glVertex3f(v[0],v[1],0);						//     Pass this coordinate to OpenGL
	}
	glEnd();											//5.3. Finish rendering the entire point-set
	
	glutSwapBuffers();
}




