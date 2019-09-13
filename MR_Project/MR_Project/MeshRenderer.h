#pragma once

#include "Renderer.h"

class UnstructuredGrid3D;



class MeshRenderer : public Renderer
{
public:

enum DRAW_STYLE
{
    DRAW_GRID     =0,					//Draw the grid only. See the SimpleRenderer class.
	DRAW_POINTS,						//Draw points only
	DRAW_C0_CELLS,						//Draw cells, using flat shading
	DRAW_C1_CELLS						//Draw cells, using smooth shading
};	

				MeshRenderer(): draw_style(DRAW_GRID) {}

void			draw(Grid&);

void			setDrawingStyle(DRAW_STYLE s)
				{ draw_style = s; }

protected:

void			drawPoints(UnstructuredGrid3D&);
void			drawGrid(UnstructuredGrid3D&);
void			drawC0Cells(UnstructuredGrid3D&);
void			drawC1Cells(UnstructuredGrid3D&);

DRAW_STYLE		draw_style;
};



						

