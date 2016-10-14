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

void LEE(GzRender* render, GzCoord* vertices, GzCoord* normals); //Hw4: add normals
void GetZPlane(const GzCoord* vertices, float* A, float* B, float* C, float* D);
float EdgeSide(const float* start, const float* end, int x, int y, bool right);
void ToScreen(const GzCoord* vert_world, GzMatrix Xsw, GzCoord* vert_screen);
void NormalToScreen(const GzCoord* normal_world, GzMatrix Xsw, GzCoord* normal_screen);
short ctoi(float color);
float GzTriangleArea(GzCoord v0, GzCoord v1, GzCoord v2);
int GzShadingEquation(GzRender *render, GzColor color, GzCoord norm);

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
	mat[1][2] = -sin(rad);
	mat[1][3] = 0;

	mat[2][0] = 0;
	mat[2][1] = sin(rad);
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
	(*render)->interp_mode = GZ_RGB_COLOR;
	(*render)->numlights = 0;

	GzColor Ka = DEFAULT_AMBIENT;
	GzColor Kd = DEFAULT_DIFFUSE;
	GzColor Ks = DEFAULT_SPECULAR;
	//initialize Ka
	(*render)->Ka[0] = Ka[0];
	(*render)->Ka[1] = Ka[1];
	(*render)->Ka[2] = Ka[2];
	//initialize Kd
	(*render)->Kd[0] = Kd[0];
	(*render)->Kd[1] = Kd[1];
	(*render)->Kd[2] = Kd[2];
	//initialize Ks
	(*render)->Ks[0] = Ks[0];
	(*render)->Ks[1] = Ks[1];
	(*render)->Ks[2] = Ks[2];

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
		for (int j = 0; j < 4; j++) {
			render->Ximage[render->matlevel][i][j] = 0;
			render->Xnorm[render->matlevel][i][j] = 0;
		}

	//HW4: add xnorm here
	if (render->matlevel == 0) {
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++) {
				render->Ximage[render->matlevel][i][j] = matrix[i][j];
				render->Xnorm[render->matlevel][i][j] = 0;
			}
		render->Xnorm[render->matlevel][0][0] = 1;
		render->Xnorm[render->matlevel][1][1] = 1;
		render->Xnorm[render->matlevel][2][2] = 1;
		render->Xnorm[render->matlevel][3][3] = 1;
	}
	else {
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				for (int m = 0; m < 4; m++)
					render->Ximage[render->matlevel][i][j] += render->Ximage[render->matlevel - 1][i][m] * matrix[m][j];
		if (render->matlevel == 1)
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					render->Xnorm[render->matlevel][i][j] = 0;
				}
			}
			render->Xnorm[render->matlevel][0][0] = 1;
			render->Xnorm[render->matlevel][1][1] = 1;
			render->Xnorm[render->matlevel][2][2] = 1;
			render->Xnorm[render->matlevel][3][3] = 1;
		}
		else
		{
			float k;
			GzMatrix R;
			//K = 1 / (a^2 + b^2 + c^2)^1/2
			k = 1 / sqrt(matrix[0][0] * matrix[0][0] + matrix[1][0] * matrix[1][0] + matrix[2][0] * matrix[2][0]);
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					R[i][j] = matrix[i][j] * k;
				}
				R[i][3] = 0;
				R[3][i] = 0;
			}
			R[3][3] = 1;
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++)
					for (int m = 0; m < 4; m++)
						render->Xnorm[render->matlevel][i][j] += render->Xnorm[render->matlevel-1][i][m] * matrix[m][j];
		}

	}
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

int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, GzPointer *valueList) /* void** valuelist */
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
		else if (nameList[i] == GZ_DIRECTIONAL_LIGHT) {
			// so this kinda assumes that the dir lights are all in one go, 
			if (render->numlights >= MAX_LIGHTS)
			{
				return GZ_FAILURE;
			}
			GzLight* dirLite = (GzLight*)valueList[i];
			render->lights[render->numlights] = *dirLite;
			render->numlights++;
		}
		else if (nameList[i] == GZ_AMBIENT_LIGHT) {
			GzLight* ambLite = (GzLight*)valueList[i];
			render->ambientlight = *ambLite;
		}
		else if (nameList[i] == GZ_DIFFUSE_COEFFICIENT) {
			GzColor* diffColor = (GzColor*)valueList[i];

			float diffR = diffColor[0][0];
			float diffG = diffColor[0][1];
			float diffB = diffColor[0][2];

			render->Kd[0] = diffR;
			render->Kd[1] = diffG;
			render->Kd[2] = diffB;
		}
		else if (nameList[i] == GZ_INTERPOLATE) {
			int* mode = (int*)valueList[i];
			render->interp_mode = *mode;
		}
		else if (nameList[i] == GZ_AMBIENT_COEFFICIENT) {
			GzColor* ambColor = (GzColor*)valueList[i];

			float ambR = ambColor[0][0];
			float ambG = ambColor[0][1];
			float ambB = ambColor[0][2];

			render->Ka[0] = ambR;
			render->Ka[1] = ambG;
			render->Ka[2] = ambB;
		}
		else if (nameList[i] == GZ_SPECULAR_COEFFICIENT) {
			GzColor* specColor = (GzColor*)valueList[i];

			float specR = specColor[0][0];
			float specG = specColor[0][1];
			float specB = specColor[0][2];

			render->Ks[0] = specR;
			render->Ks[1] = specG;
			render->Ks[2] = specB;
		}
		else if (nameList[i] == GZ_DISTRIBUTION_COEFFICIENT) {
			float* specCoeff = (float*)valueList[i]; // ugh why, int is fine!
			render->spec = *specCoeff;
		}
	}
	return GZ_SUCCESS;
}

int GzPutTriangle(GzRender	*render, int numParts, GzToken *nameList, GzPointer *valueList)
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
	
	GzCoord* vertices_screen;
	GzCoord* vertices_normal;
	vertices_screen = (GzCoord*)malloc(sizeof(GzCoord) * 3);
	vertices_normal = (GzCoord*)malloc(sizeof(GzCoord) * 3);;
	for (int i = 0; i < numParts; i++) {
		if (nameList[i] == GZ_NULL_TOKEN) {
			continue;
		}
		if (nameList[i] == GZ_POSITION) {
			GzCoord* vertices_model = (GzCoord*)valueList[i];
			ToScreen(vertices_model, render->Ximage[render->matlevel - 1], vertices_screen);
		}
		//HW4: get normals
		if (nameList[i] == GZ_NORMAL) {
			GzCoord* t_normal = (GzCoord*)valueList[i];
			normalized(t_normal[0]);
			normalized(t_normal[1]);
			normalized(t_normal[2]);
			NormalToScreen(t_normal, render->Xnorm[render->matlevel - 1], vertices_normal);
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
		LEE(render, vertices_screen, vertices_normal);
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

void LEE(GzRender* render, GzCoord* vertices, GzCoord* normals) {
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
			
			float triA = GzTriangleArea(vertices[0], vertices[1], vertices[2]);
			if (E01 > 0 && E12 > 0 && E20 > 0 || E01 < 0 && E12 < 0 && E20 < 0) {
				float pointZ = (-A * i - B * j - D) / C;
				GzGetDisplay(render->display, i, j, &red, &green, &blue, &alpha, &fbZ);
				//Hw4: ByLinear
				GzCoord p = { i, j, 1 };
				if (pointZ > 0 && pointZ < fbZ) {
					if (render->interp_mode == GZ_FLAT) {
						red = ctoi(render->flatcolor[0]);
						green = ctoi(render->flatcolor[1]);
						blue = ctoi(render->flatcolor[2]);
						fbZ = pointZ;
					}
					else if (render->interp_mode == GZ_COLOR) {
						GzColor colorV0, colorV1, colorV2;
						GzShadingEquation(render, colorV0, normals[0]);
						GzShadingEquation(render, colorV1, normals[1]);
						GzShadingEquation(render, colorV2, normals[2]);
						// Barycentric Interpolation
						// areas of each inner tris
						float A0 = GzTriangleArea(vertices[1], p, vertices[2]);
						float A1 = GzTriangleArea(vertices[0], p, vertices[2]);
						float A2 = GzTriangleArea(vertices[0], p, vertices[1]);

						// interpolate color
						float rf = (A0*colorV0[0] + A1*colorV1[0] + A2*colorV2[0]) / triA;
						float gf = (A0*colorV0[1] + A1*colorV1[1] + A2*colorV2[1]) / triA;
						float bf = (A0*colorV0[2] + A1*colorV1[2] + A2*colorV2[2]) / triA;
						if (rf > 1.0) rf = 1.0;
						if (gf > 1.0) gf = 1.0;
						if (bf > 1.0) bf = 1.0;
						red = (GzIntensity)ctoi(rf);
						green = (GzIntensity)ctoi(gf);
						blue = (GzIntensity)ctoi(bf);
						fbZ = pointZ;
					}
					else if (render->interp_mode == GZ_NORMALS) {
						// Barycentric Interpolation
						// areas of inner tris
						float A0 = GzTriangleArea(vertices[1], p, vertices[2]);
						float A1 = GzTriangleArea(vertices[0], p, vertices[2]);
						float A2 = GzTriangleArea(vertices[0], p, vertices[1]);

						// interpolate Normal of this point
						GzCoord pN;
						pN[X] = (A0*normals[0][X] + A1*normals[1][X] + A2*normals[2][X]) / triA;
						pN[Y] = (A0*normals[0][Y] + A1*normals[1][Y] + A2*normals[2][Y]) / triA;
						pN[Z] = (A0*normals[0][Z] + A1*normals[1][Z] + A2*normals[2][Z]) / triA;
						normalized(pN);

						// calculate color
						GzColor color;
						GzShadingEquation(render, color, pN);

						red = (GzIntensity)ctoi(color[0]);
						green = (GzIntensity)ctoi(color[1]);
						blue = (GzIntensity)ctoi(color[2]);
						fbZ = pointZ;
					}
					
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

void NormalToScreen(const GzCoord* normal_world, GzMatrix Xsw, GzCoord* normal_screen) {
	/*float outVect[4][3];
	float transVect[4][3];

	for (int i = 0; i < 3; i++){
		transVect[0][i] = normal_world[i][0];
		transVect[1][i] = normal_world[i][1];
		transVect[2][i] = normal_world[i][2];
		transVect[3][i] = 1;
	}

	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 3; j++){
			outVect[i][j] = 0;
		}
	}

	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 4; k++){
				outVect[i][j] += Xsw[i][k] * transVect[k][j];
			}
		}
	}

	for (int i = 0; i < 3; i++){
		normal_screen[i][0] = outVect[0][i] / outVect[3][i];
		normal_screen[i][1] = outVect[1][i] / outVect[3][i];
		normal_screen[i][2] = outVect[2][i] / outVect[3][i];
	}*/

	/*for (int i = 0; i < 3; i++){
		normalized(normal_screen[i]);
	}*/
	for (int j = 0; j < 3; ++j) {
		normal_screen[j][X] = Xsw[0][0] * normal_world[j][X] + Xsw[0][1] * normal_world[j][Y] + Xsw[0][2] * normal_world[j][Z];
		normal_screen[j][Y] = Xsw[1][0] * normal_world[j][X] + Xsw[1][1] * normal_world[j][Y] + Xsw[1][2] * normal_world[j][Z];
		normal_screen[j][Z] = Xsw[2][0] * normal_world[j][X] + Xsw[2][1] * normal_world[j][Y] + Xsw[2][2] * normal_world[j][Z];
	}
}


short ctoi(float color)		/* convert float color to GzIntensity short */
{
  return(short)((int)(color * ((1 << 12) - 1)));
}

int getColor(GzRender *render,GzColor color,GzCoord norm)
{
	
	return GZ_SUCCESS;
}

//Hw4:
float GzTriangleArea(GzCoord v0, GzCoord v1, GzCoord v2) {
	float a = .5 * (v0[X] * v1[Y] + v0[Y] * v2[X] + v1[X] * v2[Y] - v1[Y] * v2[X] - v0[Y] * v1[X] - v0[X] * v2[Y]);
	if (a < 0) return -a;
	else return a;
}

//Hw4:
float dotProduct(GzCoord v1, GzCoord v2) {
	return v1[X] * v2[X] + v1[Y] * v2[Y] + v1[Z] * v2[Z];
}

int GzShadingEquation(GzRender *render, GzColor color, GzCoord norm) {
	// computer color at each vertex
	normalized(norm);
	// E should just be camera lookat reversed? Actually no goddamn it
	GzCoord E;
	E[X] = 0;
	E[Y] = 0;
	E[Z] = -1;
	normalized(E);
	// calculate Rs for each point
	GzCoord* R = new GzCoord[render->numlights];
	// vertex 0
	float NdotL;
	int* liteCases = new int[render->numlights];
	// check N dot L and N dot E
	// if both positive fine liteCases[i] = 1;
	// if both negative flip normal, liteCases[i] = -1;
	// if different signs, skip, liteCases = 0;

	float NdotE = dotProduct(norm, E);
	for (int i = 0; i < render->numlights; ++i) {
		NdotL = dotProduct(norm, render->lights[i].direction);
		if (NdotL >= 0 && NdotE >= 0) {
			liteCases[i] = 1;
			R[i][X] = 2 * NdotL*norm[X] - render->lights[i].direction[X];
			R[i][Y] = 2 * NdotL*norm[Y] - render->lights[i].direction[Y];
			R[i][Z] = 2 * NdotL*norm[Z] - render->lights[i].direction[Z];
			normalized(R[i]);
		}
		else if (NdotL < 0 && NdotE < 0) {
			liteCases[i] = -1;
			R[i][X] = 2 * NdotL*(-norm[X]) - render->lights[i].direction[X];
			R[i][Y] = 2 * NdotL*(-norm[Y]) - render->lights[i].direction[Y];
			R[i][Z] = 2 * NdotL*(-norm[Z]) - render->lights[i].direction[Z];
			normalized(R[i]);
		}
		else {
			liteCases[i] = 0;
			continue;
		}
	}
	// check N dot L and N dot E, if both positi


	// sum all lights for Specular
	// Ks * sigma (le * (R dot E)^s) 
	GzColor specLightSum = { 0, 0, 0 };
	for (int i = 0; i < render->numlights; ++i) {
		if (liteCases[i] == 0) continue;
		float RdotE = dotProduct(R[i], E);
		if (RdotE < 0) RdotE = 0;
		if (RdotE > 1) RdotE = 1;
		// R
		specLightSum[0] += render->lights[i].color[0] * pow(RdotE, render->spec);
		// G
		specLightSum[1] += render->lights[i].color[1] * pow(RdotE, render->spec);
		// B
		specLightSum[2] += render->lights[i].color[2] * pow(RdotE, render->spec);
	}
	GzColor specComp;
	specComp[0] = render->Ks[0] * specLightSum[0]; // R
	specComp[1] = render->Ks[1] * specLightSum[1]; // G
	specComp[2] = render->Ks[2] * specLightSum[2]; // B

	// sum all lights for Diffuse
	// Kd * sigma (le * N dot L)
	
	GzColor diffLightSum = { 0, 0, 0 };
	for (int i = 0; i < render->numlights; ++i) {
		if (liteCases[i] == 0) continue;
		if (liteCases[i] == 1) {
			// R
			diffLightSum[0] += render->lights[i].color[0] * dotProduct(norm, render->lights[i].direction);
			// G
			diffLightSum[1] += render->lights[i].color[1] * dotProduct(norm, render->lights[i].direction);
			// B
			diffLightSum[2] += render->lights[i].color[2] * dotProduct(norm, render->lights[i].direction);
		}
		else if (liteCases[i] == -1) {
			GzCoord negN = { -norm[X], -norm[Y], -norm[Z] };
			// R
			diffLightSum[0] += render->lights[i].color[0] * dotProduct(negN, render->lights[i].direction);
			// G
			diffLightSum[1] += render->lights[i].color[1] * dotProduct(negN, render->lights[i].direction);
			// B
			diffLightSum[2] += render->lights[i].color[2] * dotProduct(negN, render->lights[i].direction);
		}
	}
	GzColor diffComp;
	diffComp[0] = render->Kd[0] * diffLightSum[0]; // R
	diffComp[1] = render->Kd[1] * diffLightSum[1]; // G
	diffComp[2] = render->Kd[2] * diffLightSum[2]; // B

												   // Ambient Component
	GzColor ambComp;
	ambComp[0] = render->Ka[0] * render->ambientlight.color[0]; // R
	ambComp[1] = render->Ka[1] * render->ambientlight.color[1]; // G
	ambComp[2] = render->Ka[2] * render->ambientlight.color[2]; // B

	color[0] = specComp[0] + diffComp[0] + ambComp[0]; // R
	color[1] = specComp[1] + diffComp[1] + ambComp[1]; // G
	color[2] = specComp[2] + diffComp[2] + ambComp[2]; // B

	return GZ_SUCCESS;
}