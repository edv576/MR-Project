#include <GLUT/glut.h>									//GLUT library
#include "StructuredGrid.h"
#include "SimpleRenderer.h"
#include <stdlib.h>										//For random functions



float Xmin	= -5, Xmax	= 5;							//Range of variable x that we sample
float Ymin	= -5, Ymax	= 5;							//Range of variable y that we sample
int   Nx	= 20, Ny	= 15;							//Number of samples taken along the x and y axes


float fov	 = 60;										//Perspective projection parameters
float aspect = 1; 
float z_near = 0.1; 
float z_far  = 30;

float eye_x = 10, eye_y = 10, eye_z = 10;				//Modelview (camera extrinsic) parameters
float c_x   = 0,  c_y   = 0,  c_z   = 0;
float up_x  = 0,  up_y  = 0,  up_z  = 1;


SimpleRenderer	renderer;
StructuredGrid  grid(Nx,Ny);



void viewing(int W, int H)								//Window resize function, sets up viewing parameters (GLUT callback function)
{
	glMatrixMode(GL_MODELVIEW);							//1. Set the modelview matrix (including the camera position and view direction)
	glLoadIdentity ();										
	gluLookAt(eye_x,eye_y,eye_z,c_x,c_y,c_z,up_x,up_y,up_z);

	glMatrixMode (GL_PROJECTION);						//2. Set the projection matrix
	glLoadIdentity ();
	gluPerspective(fov,float(W)/H,z_near,z_far);

	glViewport(0,0,W,H);								//3. Set the viewport to the entire size of the rendering window
}



void draw()												//Render the height plot (GLUT callback function)
{
	renderer.draw(grid);
}


void createGrid()										//Create a structured grid. For illustration purposes, we create
{														//a grid whose vertices are those of an uniform grid covering the same domain,
	srand(time(0));										//but slightly jittered in the X and Y directions

	int counter=0;										//Vertex index in the storage-order
	for(int j=0;j<Ny;++j)
	  for(int i=0;i<Nx;++i)
	  {
		float p[2];

		float factor_x = 0.3*(Xmax-Xmin)/Nx;			//Maximal jitter factors are a fraction of the uniform-grid cell size
		float factor_y = 0.3*(Ymax-Ymin)/Ny; 
		float rx = factor_x*float(rand())/RAND_MAX;		//Compute jitter amounts along X and Y
		float ry = factor_y*float(rand())/RAND_MAX;

		p[0] = Xmin + i*(Xmax-Xmin)/Nx + rx;				
		p[1] = Ymin + j*(Ymax-Ymin)/Ny + ry;
		grid.setPoint(counter++,p);						//Set vertex coordinates of next vertex in the grid
	  }
}	  

int main(int argc,char* argv[])							//Main program
{
   glutInit(&argc, argv);								//1. Initialize the GLUT toolkit
   glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);			//2. Ask GLUT to create next windows with a RGB framebuffer and a Z-buffer too
   glutInitWindowSize(500,500);							//3. Tell GLUT how large are the windows we want to create next
   glutCreateWindow("3. Structured grid");				//4. Create our window
   
   createGrid();										//5. Create the structured grid
   
   glutDisplayFunc(draw);								//6. Add a drawing callback to the window	
   glutReshapeFunc(viewing);							//7. Add a resize callback to the window
   glutMainLoop();										//8. Start the event loop that displays the graph and handles window-resize events
   
   return 0;
}




