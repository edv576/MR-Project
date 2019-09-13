#include <math.h>										//For abs(), pow()
#include <GLUT/glut.h>									//GLUT library
#include "RectilinearGrid.h"
#include "SimpleRenderer.h"



float Xmin	= -5, Xmax	= 5;							//Range of variable x that we sample
float Ymin	= -5, Ymax	= 5;							//Range of variable y that we sample
int   Nx	= 40, Ny	= 30;							//Number of samples taken along the x and y axes


float fov	 = 60;										//Perspective projection parameters
float aspect = 1; 
float z_near = 0.1; 
float z_far  = 30;

float eye_x = 10, eye_y = 10, eye_z = 10;				//Modelview (camera extrinsic) parameters
float c_x   = 0,  c_y   = 0,  c_z   = 0;
float up_x  = 0,  up_y  = 0,  up_z  = 1;



RectilinearGrid*	grid = 0;
SimpleRenderer		renderer;




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
	renderer.draw(*grid);
}



RectilinearGrid* createGrid()							//Construct a rectilinear grid
{
	vector<float> samplesX(Nx);
	vector<float> samplesY(Ny);
	
	float xc = (Xmax+Xmin)/2;							//x coordinate of grid center
	float xh = xc-Xmin;									//half-size of grid x dimension
	float yc = (Ymax+Ymin)/2;							//x coordinate of grid center
	float yh = yc-Ymin;									//half-size of grid x dimension
	float K  = 2;										//parameter determining the increase-of-density of grid points towards its center
	
	for(int i=0;i<Nx;++i)
	{
	   float x  = Xmin + i*(Xmax-Xmin)/Nx;				//in [Xmin,Xmax]
	   float nx = pow(fabs(xc-x)/xh,K);					//normalized distance (in [0,1]) to xc
	   samplesX[i] = (x<xc)? xc-xh*nx : xc+xh*nx;
	}   

	for(int i=0;i<Ny;++i)
	{
	   float y  = Ymin + i*(Ymax-Ymin)/Ny;				//in [Ymin,Ymax]
	   float ny = pow(fabs(yc-y)/yh,K);					//normalized distance (in [0,1]) to yc
	   samplesY[i] = (y<yc)? yc-yh*ny : yc+yh*ny;
	}   
	   

	RectilinearGrid* grid = new RectilinearGrid(samplesX,samplesY);

	return grid;
}



int main(int argc,char* argv[])							//Main program
{
   glutInit(&argc, argv);								//1. Initialize the GLUT toolkit
   glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);			//2. Ask GLUT to create next windows with a RGB framebuffer and a Z-buffer too
   glutInitWindowSize(500,500);							//3. Tell GLUT how large are the windows we want to create next
   glutCreateWindow("2. Rectilinear grid");				//4. Create our window
   
   grid = createGrid();									//5. Create a rectilinear grid
   
   glutDisplayFunc(draw);								//6. Add a drawing callback to the window	
   glutReshapeFunc(viewing);							//7. Add a resize callback to the window
   glutMainLoop();										//8. Start the event loop that displays the graph and handles window-resize events
   
   return 0;
}




