#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

void SwapCoord(float* v1, float* v2)
{
	float tempx = v1[X];
	float tempy = v1[Y];
	float tempz = v1[Z];
	v1[X] = v2[X];
	v1[Y] = v2[Y];
	v1[Z] = v2[Z];
	v2[X] = tempx;
	v2[Y] = tempy;
	v2[Z] = tempz;
}

void SetupTri(GzCoord* vertices);
void LEE(GzRender* render, GzCoord* vertices);
float EdgeSide(const float* start, const float* end, int x, int y, bool right);
void GetZPlane(const GzCoord* vertice, float* A, float* B, float* C, float* D);
float InterpolateZ(const GzCoord* vertices, int x, int y);

int GzNewRender(GzRender **render, GzDisplay *display)
{
/* 
- malloc a renderer struct
- span interpolator needs pointer to display for pixel writes
*/
	*render = (GzRender*)malloc(sizeof(GzRender));
	(*render)->display = display;
	return GZ_SUCCESS;
}


int GzFreeRender(GzRender *render)
{
/* 
-free all renderer resources
*/
	free(render);
	return GZ_SUCCESS;
}


int GzBeginRender(GzRender	*render)
{
/* 
- set up for start of each frame - init frame buffer
*/
	for (int i = 0; i < render->display->xres * render->display->yres; i++)
	{
		render->display->fbuf[i].alpha = 1;
		render->display->fbuf[i].blue = 100;
		render->display->fbuf[i].green = 100;
		render->display->fbuf[i].red = 100;
		render->display->fbuf[i].z = MAXINT;
	}
	return GZ_SUCCESS;
}


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer *valueList) /* void** valuelist */
{
/*
- set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
- later set shaders, interpolaters, texture maps, and lights
*/
	for (int i = 0; i < numAttributes; i++) {
		if (*nameList == GZ_RGB_COLOR) {
			GzColor* a = (GzColor*)valueList[0];
			render->flatcolor[0] = (*a)[0];
			render->flatcolor[1] = (*a)[1];
			render->flatcolor[2] = (*a)[2];
		}
	}
	return GZ_SUCCESS;
}


int GzPutTriangle(GzRender *render, int	numParts, GzToken *nameList,
	GzPointer *valueList) 
/* numParts - how many names and values */
{
/* 
- pass in a triangle description with tokens and values corresponding to
      GZ_NULL_TOKEN:		do nothing - no values
      GZ_POSITION:		3 vert positions 
- Invoke the scan converter and return an error code
*/
	GzCoord* vertices;
	for (int i = 0; i < numParts; i++){
		if (nameList[i] == GZ_NULL_TOKEN){
			continue;
		}
		if (nameList[i] == GZ_POSITION){
			vertices = (GzCoord*)valueList[i];
		}
	}
	SetupTri(vertices);
	LEE(render, vertices);
	return GZ_SUCCESS;
}

/* NOT part of API - just for general assistance */
short	ctoi(float color)		/* convert float color to GzIntensity short */
{
  return(short)((int)(color * ((1 << 12) - 1)));
}

//Make CW edges always 0-1,1-2,2-0
void SetupTri(GzCoord* vertices) {
	//Determine Top/Bot relationship
	for (int i = 2; i > 0; i--) {
		if (vertices[i][Y] < vertices[i - 1][Y]) {
			SwapCoord(vertices[i], vertices[i - 1]);
		}
	}
	if (vertices[2][Y] < vertices[1][Y]) {
		SwapCoord(vertices[2], vertices[1]);
	}
	//Determine L/R relationship
	if (vertices[0][Y] == vertices[1][Y]) {
		if (vertices[0][X] > vertices[1][X]) {
			SwapCoord(vertices[0], vertices[1]);
		}
	}
	else if (vertices[1][Y] == vertices[2][Y]) {
		if (vertices[1][X] < vertices[2][X]) {
			SwapCoord(vertices[1], vertices[2]);
		}
	}
	else if ((vertices[0][X] - vertices[2][X])*(vertices[1][Y] - vertices[2][Y]) / (vertices[0][Y] - vertices[2][Y]) + vertices[2][X] > vertices[1][X]) {
		SwapCoord(vertices[1], vertices[2]);
	}
}

void LEE(GzRender* render, GzCoord* vertices) {
	//Get Bounding Box
	int Up = floor(vertices[0][Y]);
	int Down = ceil(vertices[1][Y] > vertices[2][Y] ? vertices[1][Y] : vertices[2][Y]);
	int Left = floor(min(min(vertices[0][X], vertices[1][X]), vertices[2][X]));
	int Right = ceil(max(max(vertices[0][X], vertices[1][X]), vertices[2][X]));
    
	//Determine the right edge
	bool E01Right;
	bool E12Right;
	bool E20Right;
    
	//Get Display Parameter
	GzIntensity red, green, blue, alpha;
	GzDepth fbZ;

	if (vertices[0][Y] == vertices[1][Y] || vertices[1][Y] < vertices[2][Y])
	{
		E01Right = true;
		E12Right = true;
		E20Right = false;
	}
	else
	{
		E01Right = true;
		E12Right = false;
		E20Right = false;
	}

	//Calculate Z Plane for interpolation
	float A, B, C, D;
	GetZPlane(vertices, &A, &B, &C, &D);

	for (int j = Up; j < Down; j++) {
		if (j < 0 || j > render->display->xres) {
			continue;
		}
		for (int i = Left; i < Right; i++) {
			if (i < 0 || i > render->display->yres) {
				continue;
			}
			float E01 = EdgeSide(vertices[0], vertices[1], i, j, E01Right);
			float E12 = EdgeSide(vertices[1], vertices[2], i, j, E12Right);
			float E20 = EdgeSide(vertices[2], vertices[0], i, j, E20Right);	
			if (E01 > 0 && E12 > 0 && E20 > 0 || E01 < 0 && E12 < 0 && E20 < 0) {
				float pointZ = (-A * i - B * j - D) / C;
				GzGetDisplay(render->display, i, j, &red, &green, &blue, &alpha, &fbZ);
				if (pointZ > 0 && pointZ < fbZ){
					red = ctoi(render->flatcolor[0]);
					green = ctoi(render->flatcolor[1]);
					blue = ctoi(render->flatcolor[2]);
					fbZ = pointZ;
					GzPutDisplay(render->display, i, j, red, green, blue, alpha, fbZ);
				}
			}
		}
	}
}

float EdgeSide(const float* start, const float* end, int x, int y, bool right) {
	return (end[Y] - start[Y]) * (x - start[X]) - (end[X] - start[X]) * (y - start[Y]);
}

void GetZPlane(const GzCoord* vertices, float* A, float* B, float* C, float* D){
	GzCoord E01;
	E01[X] = vertices[1][X] - vertices[0][X];
	E01[Y] = vertices[1][Y] - vertices[0][Y];
	E01[Z] = vertices[1][Z] - vertices[0][Z];
	GzCoord E12;
	E12[X] = vertices[2][X] - vertices[1][X];
	E12[Y] = vertices[2][Y] - vertices[1][Y];
	E12[Z] = vertices[2][Z] - vertices[1][Z];

	*A = E01[Y] * E12[Z] - E01[Z] * E12[Y];
	*B = E01[Z] * E12[X] - E01[X] * E12[Z];
	*C = E01[X] * E12[Y] - E01[Y] * E12[X];

	*D = - *A * (vertices[0][X]) - *B * (vertices[0][Y]) - *C * (vertices[0][Z]);
}
