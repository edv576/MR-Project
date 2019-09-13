#include <GL/glut.h>									//GLUT library
#include "UniformGrid.h"
#include "SimpleRenderer.h"



float Xmin	= -5, Xmax	= 5;							//Range of variable x that we sample
float Ymin	= -5, Ymax	= 5;							//Range of variable y that we sample
int   Nx	= 30, Ny	= 20;							//Number of samples taken along the x and y axes


float fov	 = 60;										//Perspective projection parameters
float aspect = 1; 
float z_near = 0.1; 
float z_far  = 30;

float eye_x = 10, eye_y = 10, eye_z = 10;				//Modelview (camera extrinsic) parameters
float c_x   = 0,  c_y   = 0,  c_z   = 0;
float up_x  = 0,  up_y  = 0,  up_z  = 1;


UniformGrid		grid(Nx,Ny,Xmin,Ymin,Xmax,Ymax);
SimpleRenderer	renderer;




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



int main(int argc,char* argv[])							//Main program
{
   glutInit(&argc, argv);								//1. Initialize the GLUT toolkit
   glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);			//2. Ask GLUT to create next windows with a RGB framebuffer and a Z-buffer too
   glutInitWindowSize(500,500);							//3. Tell GLUT how large are the windows we want to create next
   glutCreateWindow("1. Uniform grid");					//4. Create our window
   
   glutDisplayFunc(draw);								//5. Add a drawing callback to the window	
   glutReshapeFunc(viewing);							//6. Add a resize callback to the window
   glutMainLoop();										//7. Start the event loop that displays the graph and handles window-resize events
   
   return 0;
}




