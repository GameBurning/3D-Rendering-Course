/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

#define PI 3.1415926
int setupXsp(GzRender *render);
int setupXpi(GzRender *render);
int setupXiw(GzRender *render);
int initCamera(GzRender *render);
int normalized(GzCoord vector);
void SetupTri(GzCoord* vertices);
void SwapCoord(float* v1, float* v2);
void LEE(GzRender* render, GzCoord* vertices);
void GetZPlane(const GzCoord* vertices, float* A, float* B, float* C, float* D);
float EdgeSide(const float* start, const float* end, int x, int y, bool right);
void ToScreen(const GzCoord* vert_world, GzMatrix Xsw, GzCoord* vert_screen);
short ctoi(float color);


int GzRotXMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value
	if (mat == NULL)
		return GZ_FAILURE;
	float rad = degree * (PI / 180.0);

	mat[0][0] = 1;
	mat[0][1] = 0;
	mat[0][2] = 0;
	mat[0][3] = 0;

	mat[1][0] = 0;
	mat[1][1] = cos(rad);
	mat[1][2] = sin(rad);
	mat[1][3] = 0;

	mat[2][0] = 0;
	mat[2][1] = -sin(rad);
	mat[2][2] = cos(rad);
	mat[2][3] = 0;

	mat[3][0] = 0;
	mat[3][1] = 0;
	mat[3][2] = 0;
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRotYMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value
	if (mat == NULL)
		return GZ_FAILURE;
	float rad = degree * (PI / 180.0);
	mat[0][0] = cos(rad);
	mat[0][1] = 0;
	mat[0][2] = sin(rad);
	mat[0][3] = 0;

	mat[1][0] = 0;
	mat[1][1] = 1;
	mat[1][2] = 0;
	mat[1][3] = 0;

	mat[2][0] = -sin(rad);
	mat[2][1] = 0;
	mat[2][2] = cos(rad);
	mat[2][3] = 0;

	mat[3][0] = 0;
	mat[3][1] = 0;
	mat[3][2] = 0;
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRotZMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value
	if (mat == NULL)
		return GZ_FAILURE;
	float rad = degree * (PI / 180.0);
	mat[0][0] = cos(rad);
	mat[0][1] = -sin(rad);
	mat[0][2] = 0;
	mat[0][3] = 0;

	mat[1][0] = sin(rad);
	mat[1][1] = cos(rad);
	mat[1][2] = 0;
	mat[1][3] = 0;

	mat[2][0] = 0;
	mat[2][1] = 0;
	mat[2][2] = 1;
	mat[2][3] = 0;

	mat[3][0] = 0;
	mat[3][1] = 0;
	mat[3][2] = 0;
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzTrxMat(GzCoord translate, GzMatrix mat)
{
	// Create translation matrix
	// Pass back the matrix using mat value
	if (mat == NULL)
		return GZ_FAILURE;
	mat[0][0] = 1;
	mat[0][1] = 0;
	mat[0][2] = 0;
	mat[0][3] = translate[X];

	mat[1][0] = 0;
	mat[1][1] = 1;
	mat[1][2] = 0;
	mat[1][3] = translate[Y];

	mat[2][0] = 0;
	mat[2][1] = 0;
	mat[2][2] = 1;
	mat[2][3] = translate[Z];

	mat[3][0] = 0;
	mat[3][1] = 0;
	mat[3][2] = 0;
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzScaleMat(GzCoord scale, GzMatrix mat)
{
	// Create scaling matrix
	// Pass back the matrix using mat value
	if (mat == NULL)
		return GZ_FAILURE;
	mat[0][0] = scale[X];
	mat[0][1] = 0;
	mat[0][2] = 0;
	mat[0][3] = 0;

	mat[1][0] = 0;
	mat[1][1] = scale[Y];
	mat[1][2] = 0;
	mat[1][3] = 0;

	mat[2][0] = 0;
	mat[2][1] = 0;
	mat[2][2] = scale[Z];
	mat[2][3] = 0;

	mat[3][0] = 0;
	mat[3][1] = 0;
	mat[3][2] = 0;
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

//----------------------------------------------------------
// Begin main functions

int GzNewRender(GzRender **render, GzDisplay *display)
{
/*  
- malloc a renderer struct 
- setup Xsp and anything only done once 
- save the pointer to display 
- init default camera 
*/ 
	*render = (GzRender*)malloc(sizeof(GzRender));

	(*render)->display = display;
	(*render)->matlevel = 0;

	//if(setupXsp(*render) && initCamera(*render))
	setupXsp(*render);
	initCamera(*render);
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

int GzBeginRender(GzRender *render)
{
/*  
- setup for start of each frame - init frame buffer color,alpha,z
- compute Xiw and projection xform Xpi from camera definition 
- init Ximage - put Xsp at base of stack, push on Xpi and Xiw 
- now stack contains Xsw and app can push model Xforms when needed 
*/ 

	//init frame buffer color
	for (int i = 0; i < render->display->xres * render->display->yres; i++)
	{
		render->display->fbuf[i].alpha = 1;
		render->display->fbuf[i].blue = 100;
		render->display->fbuf[i].green = 100;
		render->display->fbuf[i].red = 100;
		render->display->fbuf[i].z = MAXINT;
	}

	//compute Xiw
	setupXiw(render);

	//compute xpi
	setupXpi(render);



	GzPushMatrix(render, render->Xsp);
	GzPushMatrix(render, render->camera.Xpi);
	GzPushMatrix(render, render->camera.Xiw);

	return GZ_SUCCESS;
}

int GzPutCamera(GzRender *render, GzCamera *camera)
{
/*
- overwrite renderer camera structure with new camera definition
*/
	render->camera.FOV = camera->FOV;

	render->camera.position[0] = camera->position[0];
	render->camera.position[1] = camera->position[1];
	render->camera.position[2] = camera->position[2];

	render->camera.lookat[0] = camera->lookat[0];
	render->camera.lookat[1] = camera->lookat[1];
	render->camera.lookat[2] = camera->lookat[2];

	render->camera.worldup[0] = camera->worldup[0];
	render->camera.worldup[1] = camera->worldup[1];
	render->camera.worldup[2] = camera->worldup[2];
	normalized(render->camera.worldup);

	render->Xsp[2][2] = 2147483647 * tan((render->camera.FOV / 2.0) * (PI / 180.0));

	return GZ_SUCCESS;	
}

int GzPushMatrix(GzRender *render, GzMatrix	matrix)
{
/*
- push a matrix onto the Ximage stack
- check for stack overflow
*/
	if (render->matlevel >= MATLEVELS)
		return GZ_FAILURE;

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			render->Ximage[render->matlevel][i][j] = 0;

	if (render->matlevel == 0)
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				render->Ximage[render->matlevel][i][j] = matrix[i][j];
	else
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				for (int m = 0; m < 4; m++)
					render->Ximage[render->matlevel][i][j] += render->Ximage[render->matlevel - 1][i][m] * matrix[m][j];
	render->matlevel++;
	return GZ_SUCCESS;
}

int GzPopMatrix(GzRender *render)
{
/*
- pop a matrix off the Ximage stack
- check for stack underflow
*/
	if (render->matlevel <= 0)
		return GZ_FAILURE;
	render->matlevel--;
	return GZ_SUCCESS;
}

int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer	*valueList) /* void** valuelist */
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

int GzPutTriangle(GzRender	*render, int numParts, GzToken *nameList, GzPointer	*valueList)
/* numParts : how many names and values */
{
/*  
- pass in a triangle description with tokens and values corresponding to 
      GZ_POSITION:3 vert positions in model space 
- Xform positions of verts using matrix on top of stack 
- Clip - just discard any triangle with any vert(s) behind view plane 
       - optional: test for triangles with all three verts off-screen (trivial frustum cull)
- invoke triangle rasterizer  
*/ 

	GzCoord* vertices_model;
	GzCoord* vertices_screen;
	vertices_screen = (GzCoord*)malloc(sizeof(GzCoord) * 3);
	for (int i = 0; i < numParts; i++) {
		if (nameList[i] == GZ_NULL_TOKEN) {
			continue;
		}
		if (nameList[i] == GZ_POSITION) {
			vertices_model = (GzCoord*)valueList[i];
			ToScreen(vertices_model, render->Ximage[render->matlevel - 1], vertices_screen);
		}
	}

	bool inScreen = FALSE;
	for (int i = 0; i < 3; i++) {
		if (vertices_screen[i][X] >= 0 && vertices_screen[i][X] < render->display->xres &&
			vertices_screen[i][Y] >= 0 && vertices_screen[i][Y] < render->display->yres) {
			inScreen = TRUE;
			break;
		}
	}

	if (vertices_screen[0][Z] > 0 && vertices_screen[1][Z] > 0 && vertices_screen[2][Z] > 0 && inScreen) {
		SetupTri(vertices_screen);
		LEE(render, vertices_screen);
	}
	return GZ_SUCCESS;
	
}

/* NOT part of API - just for general assistance */

int setupXsp(GzRender *render)
{
	float d;
	float radian;
	radian = render->camera.FOV / 180.0 * PI;
	d = 1 / (tan(radian / 2));
	render->Xsp[0][0] = render->display->xres / 2.0;
	render->Xsp[0][1] = 0;
	render->Xsp[0][2] = 0;
	render->Xsp[0][3] = render->display->xres / 2.0;
	
	render->Xsp[1][0] = 0;
	render->Xsp[1][1] = -render->display->yres / 2.0;
	render->Xsp[1][2] = 0;
	render->Xsp[1][3] = render->display->yres / 2.0;

	render->Xsp[2][0] = 0;
	render->Xsp[2][1] = 0;
	render->Xsp[2][2] = MAXINT;
	render->Xsp[2][3] = 0;

	render->Xsp[3][0] = 0;
	render->Xsp[3][1] = 0;
	render->Xsp[3][2] = 0;
	render->Xsp[3][3] = 1;
	return GZ_SUCCESS;
}

int setupXpi(GzRender *render)
{
	float rad = (render->camera.FOV / 2.0) * (PI / 180.0);

	render->camera.Xpi[0][0] = 1;
	render->camera.Xpi[0][1] = 0;
	render->camera.Xpi[0][2] = 0;
	render->camera.Xpi[0][3] = 0;

	render->camera.Xpi[1][0] = 0;
	render->camera.Xpi[1][1] = 1;
	render->camera.Xpi[1][2] = 0;
	render->camera.Xpi[1][3] = 0;

	render->camera.Xpi[2][0] = 0;
	render->camera.Xpi[2][1] = 0;
	//ToDo Check
	render->camera.Xpi[2][2] = tan(rad);
	render->camera.Xpi[2][3] = 0;

	render->camera.Xpi[3][0] = 0;
	render->camera.Xpi[3][1] = 0;
	render->camera.Xpi[3][2] = tan(rad);
	render->camera.Xpi[3][3] = 1;

	return GZ_SUCCESS;
}


int setupXiw(GzRender *render)
{
	GzCoord cl, camZ;
	cl[X] = render->camera.lookat[X] - render->camera.position[X];
	cl[Y] = render->camera.lookat[Y] - render->camera.position[Y];
	cl[Z] = render->camera.lookat[Z] - render->camera.position[Z];
	normalized(cl);
	camZ[X] = cl[X];
	camZ[Y] = cl[Y];
	camZ[Z] = cl[Z];
	normalized(camZ);

	GzCoord camUp, camY;
	float upDotZ = render->camera.worldup[X] * camZ[X] + render->camera.worldup[Y] * camZ[Y] +
		render->camera.worldup[Z] * camZ[Z];
	camUp[X] = render->camera.worldup[X] - upDotZ*camZ[X];
	camUp[Y] = render->camera.worldup[Y] - upDotZ*camZ[Y];
	camUp[Z] = render->camera.worldup[Z] - upDotZ*camZ[Z];
	normalized(camUp);
	camY[X] = camUp[X];
	camY[Y] = camUp[Y];
	camY[Z] = camUp[Z];
	normalized(camY);

	GzCoord camX;
	camX[X] = camY[Y] * camZ[Z] - camY[Z] * camZ[Y];
	camX[Y] = camY[Z] * camZ[X] - camY[X] * camZ[Z];
	camX[Z] = camY[X] * camZ[Y] - camY[Y] * camZ[X];
	normalized(camX);

	render->camera.Xiw[0][0] = camX[X];
	render->camera.Xiw[0][1] = camX[Y];
	render->camera.Xiw[0][2] = camX[Z];
	render->camera.Xiw[0][3] = -(camX[X] * render->camera.position[X]
		+ camX[Y] * render->camera.position[Y]
		+ camX[Z] * render->camera.position[Z]);

	render->camera.Xiw[1][0] = camY[X];
	render->camera.Xiw[1][1] = camY[Y];
	render->camera.Xiw[1][2] = camY[Z];
	render->camera.Xiw[1][3] = -(camY[X] * render->camera.position[X]
		+ camY[Y] * render->camera.position[Y]
		+ camY[Z] * render->camera.position[Z]);

	render->camera.Xiw[2][0] = camZ[X];
	render->camera.Xiw[2][1] = camZ[Y];
	render->camera.Xiw[2][2] = camZ[Z];
	render->camera.Xiw[2][3] = -(camZ[X] * render->camera.position[X]
		+ camZ[Y] * render->camera.position[Y]
		+ camZ[Z] * render->camera.position[Z]);

	render->camera.Xiw[3][0] = 0;
	render->camera.Xiw[3][1] = 0;
	render->camera.Xiw[3][2] = 0;
	render->camera.Xiw[3][3] = 1;

	return GZ_SUCCESS;
}

int initCamera(GzRender *render)
{
	render->camera.FOV = DEFAULT_FOV;

	render->camera.lookat[0] = 0;
	render->camera.lookat[1] = 0;
	render->camera.lookat[2] = 0;

	render->camera.position[0] = DEFAULT_IM_X;
	render->camera.position[1] = DEFAULT_IM_Y;
	render->camera.position[2] = DEFAULT_IM_Z;

	render->camera.worldup[0] = 0;
	render->camera.worldup[1] = 1;
	render->camera.worldup[2] = 0;

	return GZ_SUCCESS;
}

int normalized(GzCoord vector) {
	float length = sqrt(vector[X] * vector[X] + vector[Y] * vector[Y] + vector[Z] * vector[Z]);
	vector[X] /= length;
	vector[Y] /= length;
	vector[Z] /= length;

	return GZ_SUCCESS;
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
				if (pointZ > 0 && pointZ < fbZ) {
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

void GetZPlane(const GzCoord* vertices, float* A, float* B, float* C, float* D) {
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

	*D = -*A * (vertices[0][X]) - *B * (vertices[0][Y]) - *C * (vertices[0][Z]);
}

float EdgeSide(const float* start, const float* end, int x, int y, bool right) {
	return (end[Y] - start[Y]) * (x - start[X]) - (end[X] - start[X]) * (y - start[Y]);
}

void ToScreen(const GzCoord* vert_world, GzMatrix Xsw, GzCoord* vert_screen)
{
	float tri_world[4][3];
	float tri_screen[4][3];

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 3; j++)
			tri_screen[i][j] = 0;

	for (int i = 0; i < 3; i++)
	{
		tri_world[X][i] = vert_world[i][X];
		tri_world[Y][i] = vert_world[i][Y];
		tri_world[Z][i] = vert_world[i][Z];
		tri_world[3][i] = 1;
	}

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 3; j++)
			for (int m = 0; m < 4; m++)
				tri_screen[i][j] += Xsw[i][m] * tri_world[m][j];

	for (int i = 0; i < 3; i++)
	{
		vert_screen[i][X] = tri_screen[X][i] / tri_screen[3][i];
		vert_screen[i][Y] = tri_screen[Y][i] / tri_screen[3][i];
		vert_screen[i][Z] = tri_screen[Z][i] / tri_screen[3][i];
	}
}

short	ctoi(float color)		/* convert float color to GzIntensity short */
{
  return(short)((int)(color * ((1 << 12) - 1)));
}