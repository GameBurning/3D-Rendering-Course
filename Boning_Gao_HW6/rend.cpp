/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
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

	//HW5
	//init texture function
	(*render)->tex_fun = NULL;

	(*render)->shiftX = (*render)->shiftY = 0;

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
						render->Xnorm[render->matlevel][i][j] += render->Xnorm[render->matlevel-1][i][m] * R[m][j];
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
		//HW5
		else if (nameList[i] == GZ_TEXTURE_MAP)
		{
			GzTexture texfunction = (GzTexture)(valueList[i]);
			render->tex_fun = texfunction;
		}
		else if (nameList[i] == GZ_AASHIFTX)
		{
			float *x_shift = (float*)(valueList[i]);
			render->shiftX = *x_shift;
		}
		else if (nameList[i] == GZ_AASHIFTY)
		{
			float *y_shift = (float*)(valueList[i]);
			render->shiftY = *y_shift;
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
	GzTextureIndex* uv_coord;
	GzCoord* vertices_model;
	vertices_screen = (GzCoord*)malloc(sizeof(GzCoord) * 3);
	vertices_normal = (GzCoord*)malloc(sizeof(GzCoord) * 3);
	uv_coord = (GzTextureIndex*)malloc(sizeof(GzTextureIndex) * 2);
	for (int i = 0; i < numParts; i++) {
		if (nameList[i] == GZ_NULL_TOKEN) {
			continue;
		}
		if (nameList[i] == GZ_POSITION) {
			vertices_model = (GzCoord*)valueList[i];
			ToScreen(vertices_model, render->Ximage[render->matlevel - 1], vertices_screen);
			for (int i = 0; i < 3; i++) {
				vertices_screen[i][X] -= render->shiftX;
				vertices_screen[i][Y] -= render->shiftY;
			}
		}
		//HW4: get normals
		if (nameList[i] == GZ_NORMAL) {
			GzCoord* t_normal = (GzCoord*)valueList[i];
			normalized(t_normal[0]);
			normalized(t_normal[1]);
			normalized(t_normal[2]);
			NormalToScreen(t_normal, render->Xnorm[render->matlevel - 1], vertices_normal);
		}
		//HW5: get uv
		if (nameList[i] == GZ_TEXTURE_INDEX){
			GzTextureIndex* uv = (GzTextureIndex*)valueList[i];
			for (int i = 0; i < 3; i++) {
				//transform uv into perspective space UV for each vertex
				float vz = vertices_screen[i][Z] / (INT_MAX - vertices_screen[i][Z]);
				uv_coord[i][0] = uv[i][0] / (vz + 1);
				uv_coord[i][1] = uv[i][1] / (vz + 1);
				//uv_coord[i][0] = uv[i][0];
				//uv_coord[i][1] = uv[i][1];
			}
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
		SetupTri(vertices_screen, vertices_normal, uv_coord);
		LEE(render, vertices_screen, vertices_normal, uv_coord);
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
	render->camera.Xpi[2][2] = 1;
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
//HW5: adds uv
void SetupTri(GzCoord* vertices, GzCoord* normals, GzTextureIndex* uvs) {
	//Determine Top/Bot relationship
	for (int i = 2; i > 0; i--) {
		if (vertices[i][Y] < vertices[i - 1][Y]) {
			SwapCoord(vertices[i], vertices[i - 1]);
			SwapCoord(normals[i], normals[i - 1]);
			SwapUV(uvs[i], uvs[i - 1]);
		}
	}
	if (vertices[2][Y] < vertices[1][Y]) {
		SwapCoord(vertices[2], vertices[1]);
		SwapCoord(normals[2], normals[1]);
		SwapUV(uvs[2], uvs[1]);
	}
	//Determine L/R relationship
	if (vertices[0][Y] == vertices[1][Y]) {
		if (vertices[0][X] > vertices[1][X]) {
			SwapCoord(vertices[0], vertices[1]);
			SwapCoord(normals[0], normals[1]);
			SwapUV(uvs[0], uvs[1]);
		}
	}
	else if (vertices[1][Y] == vertices[2][Y]) {
		if (vertices[1][X] < vertices[2][X]) {
			SwapCoord(vertices[1], vertices[2]);
			SwapCoord(normals[1], normals[2]);
			SwapUV(uvs[1], uvs[2]);
		}
	}
	else if ((vertices[0][X] - vertices[2][X])*(vertices[1][Y] - vertices[2][Y]) / (vertices[0][Y] - vertices[2][Y]) + vertices[2][X] > vertices[1][X]) {
		SwapCoord(vertices[1], vertices[2]);
		SwapCoord(normals[1], normals[2]);
		SwapUV(uvs[1], uvs[2]);
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

void SwapUV(float* v1, float* v2)
{
	float tempU = v1[0];
	float tempV = v1[1];
	v1[0] = v2[0];
	v1[1] = v2[1];
	v2[0] = tempU;
	v2[1] = tempV;
}

void LEE(GzRender* render, GzCoord* vertices, GzCoord* normals, GzTextureIndex* uvs) {
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
	//This perspective not work
	//float A, B, C, D;
	//GetZPlane(vertices, &A, &B, &C, &D);

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
				GzCoord p = { i, j, 1 };
				// Barycentric Interpolation
				// areas of each inner tris
				float A0 = triangleArea(vertices[1], p, vertices[2]);
				float A1 = triangleArea(vertices[0], p, vertices[2]);
				float A2 = triangleArea(vertices[0], p, vertices[1]);

				float x0 = vertices[0][0];
				float x1 = vertices[1][0];
				float x2 = vertices[2][0];
				float y0 = vertices[0][1];
				float y1 = vertices[1][1];
				float y2 = vertices[2][1];
				float z0 = vertices[0][2];
				float z1 = vertices[1][2];
				float z2 = vertices[2][2];

				float div = ((y1 - y2)*(x0 - x2) + (x2 - x1)*(y0 - y2)); // calculating the denominator in order to use for the calculation of coefficients, alpha,beta,gamma.
				float coeff1 = ((y1 - y2)*(i - x2) + (x2 - x1)*(j - y2)) / div;
				float coeff2 = ((y2 - y0)*(i - x2) + (x0 - x2)*(j - y2)) / div;
				float c_val = 1 - coeff1 - coeff2;

				float triA = triangleArea(vertices[0], vertices[1], vertices[2]);
				float pointZ = (A0 * vertices[0][Z] + A1 * vertices[1][Z] + A2*vertices[2][Z]) / triA;
				//float pointZ = z0 * coeff1 + z1 * coeff2 + z2 * c_val;

				GzGetDisplay(render->display, i, j, &red, &green, &blue, &alpha, &fbZ);
				
				if (pointZ > 0 && pointZ < fbZ) {
					if (render->interp_mode == GZ_FLAT) {
						red = ctoi(render->flatcolor[0]);
						green = ctoi(render->flatcolor[1]);
						blue = ctoi(render->flatcolor[2]);
						fbZ = pointZ;
					}

					else if (render->interp_mode == GZ_COLOR) {
						GzColor colorV0, colorV1, colorV2;
						CalculateColorTexture(render, colorV0, normals[0]);
						CalculateColorTexture(render, colorV1, normals[1]);
						CalculateColorTexture(render, colorV2, normals[2]);

						// HW5
						// Interpolate UV
						GzTextureIndex UV;
						UV[0] = (A0*uvs[0][0] + A1*uvs[1][0] + A2*uvs[2][0]) / triA;
						UV[1] = (A0*uvs[0][1] + A1*uvs[1][1] + A2*uvs[2][1]) / triA;

						//UV[0] = (uvs[0][0] + uvs[1][0] + uvs[2][0]) / 3;
						//UV[1] = (uvs[0][1] + uvs[1][1] + uvs[2][1]) / 3;

						//float atPointz0 = z0 / (INT_MAX - z0);	// Warping as per z0
						//float atPointz1 = z1 / (INT_MAX - z1);	// Warping as per z1
						//float atPointz2 = z2 / (INT_MAX - z2);	// Warping as per z2
						//float atPointzz = pointZ / (INT_MAX - pointZ);

						//float u0 = uvs[0][0] / (atPointz0 + 1);
						//float u1 = uvs[1][0] / (atPointz1 + 1);
						//float u2 = uvs[2][0] / (atPointz2 + 1);

						//float v0 = uvs[0][1] / (atPointz0 + 1);
						//float v1 = uvs[1][1] / (atPointz1 + 1);
						//float v2 = uvs[2][1] / (atPointz2 + 1);


						//UV[0] = coeff1 * u0 + coeff2 * u1 + c_val * u2;
						//UV[1] = coeff1 * v0 + coeff2 * v1 + c_val * v2;

						GzTextureIndex uv;
						float vz = pointZ / (INT_MAX - pointZ);

						uv[0] = UV[0] * (vz + 1);
						uv[1] = UV[1] * (vz + 1);

						//uv[0] = UV[0];
						//uv[1] = UV[1];

						// Get Texture Color
						GzColor textureColor;
						if (render->tex_fun != NULL) {
							render->tex_fun(uv[0], uv[1], textureColor);
						}

						// interpolate color
						float interp_r = (A0*colorV0[0] + A1*colorV1[0] + A2*colorV2[0]) / triA;
						float interp_g = (A0*colorV0[1] + A1*colorV1[1] + A2*colorV2[1]) / triA;
						float interp_b = (A0*colorV0[2] + A1*colorV1[2] + A2*colorV2[2]) / triA;
						if (render->tex_fun != NULL) {
							interp_r *= textureColor[RED];
							interp_b *= textureColor[BLUE];
							interp_g *= textureColor[GREEN];
						}

						if (interp_r > 1.0) interp_r = 1.0;
						if (interp_g > 1.0) interp_g = 1.0;
						if (interp_b > 1.0) interp_b = 1.0;
						red = (GzIntensity)ctoi(interp_r);
						green = (GzIntensity)ctoi(interp_g);
						blue = (GzIntensity)ctoi(interp_b);
						fbZ = pointZ;
					}
					
					else if (render->interp_mode == GZ_NORMALS) {
						// byinterpolate
						float A0 = triangleArea(vertices[1], p, vertices[2]);
						float A1 = triangleArea(vertices[0], p, vertices[2]);
						float A2 = triangleArea(vertices[0], p, vertices[1]);

						GzCoord interp_N;
						interp_N[X] = (A0*normals[0][X] + A1*normals[1][X] + A2*normals[2][X]) / triA;
						interp_N[Y] = (A0*normals[0][Y] + A1*normals[1][Y] + A2*normals[2][Y]) / triA;
						interp_N[Z] = (A0*normals[0][Z] + A1*normals[1][Z] + A2*normals[2][Z]) / triA;
						normalized(interp_N);

						// HW5
						// Interpolate UV
						GzTextureIndex UV;
						UV[0] = (A0*uvs[0][0] + A1*uvs[1][0] + A2*uvs[2][0]) / triA;
						UV[1] = (A0*uvs[0][1] + A1*uvs[1][1] + A2*uvs[2][1]) / triA;

						GzTextureIndex uv;
						float uv_z = pointZ / (INT_MAX - pointZ);
						uv[0] = UV[0] * (uv_z + 1);
						uv[1] = UV[1] * (uv_z + 1);
						//uv[0] = UV[0];
						//uv[1] = UV[1];

						// Get Texture Color
						GzColor textureColor;
						if (render->tex_fun != NULL) {
							render->tex_fun(uv[0], uv[1], textureColor);
						}

						if (render->tex_fun != NULL) {
							render->Ka[RED] = render->Kd[RED] = textureColor[RED];
							render->Ka[GREEN] = render->Kd[GREEN] = textureColor[GREEN];
							render->Ka[BLUE] = render->Kd[BLUE] = textureColor[BLUE];
						}

						// calculate color
						GzColor color;
						
						CalculateColor(render, color, interp_N);

						if (color[0] > 1.0) color[0] = 1.0;
						if (color[1] > 1.0) color[1] = 1.0;
						if (color[2] > 1.0) color[2] = 1.0;

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
	float tempNorm[4][3];
	float transNorm[4][3];

	for (int i = 0; i < 3; i++){
		transNorm[0][i] = normal_world[i][0];
		transNorm[1][i] = normal_world[i][1];
		transNorm[2][i] = normal_world[i][2];
		transNorm[3][i] = 1;
	}

	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 3; j++){
			tempNorm[i][j] = 0;
		}
	}

	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 4; k++){
				tempNorm[i][j] += Xsw[i][k] * transNorm[k][j];
			}
		}
	}

	for (int i = 0; i < 3; i++){
		normal_screen[i][0] = tempNorm[0][i] / tempNorm[3][i];
		normal_screen[i][1] = tempNorm[1][i] / tempNorm[3][i];
		normal_screen[i][2] = tempNorm[2][i] / tempNorm[3][i];
	}

	for (int i = 0; i < 3; i++){
		normalized(normal_screen[i]);
	}
}

short ctoi(float color)		/* convert float color to GzIntensity short */
{
  return(short)((int)(color * ((1 << 12) - 1)));
}

//Hw4:
float triangleArea(GzCoord v0, GzCoord v1, GzCoord v2) {
	float a = .5 * (v0[X] * v1[Y] + v0[Y] * v2[X] + v1[X] * v2[Y] - v1[Y] * v2[X] - v0[Y] * v1[X] - v0[X] * v2[Y]);
	if (a < 0) return -a;
	else return a;
}

//Hw4:
float dotProduct(GzCoord v1, GzCoord v2) {
	return v1[X] * v2[X] + v1[Y] * v2[Y] + v1[Z] * v2[Z];
}

int CalculateColor(GzRender *render, GzColor color, GzCoord norm) {
	normalized(norm);
	GzCoord* reflect = (GzCoord*)malloc(sizeof(GzCoord) * render->numlights);
	float NdotL;
	int* checkLights = (int*)malloc(sizeof(int)*render->numlights);
	GzCoord eye;
	eye[X] = 0;
	eye[Y] = 0;
	eye[Z] = -1;
	normalized(eye);
	float NdotE = dotProduct(norm, eye);
	for (int i = 0; i < render->numlights; i++) {
		NdotL = dotProduct(norm, render->lights[i].direction);
		if (NdotL >= 0 && NdotE >= 0) {
			checkLights[i] = 1;
			reflect[i][X] = 2 * NdotL*norm[X] - render->lights[i].direction[X];
			reflect[i][Y] = 2 * NdotL*norm[Y] - render->lights[i].direction[Y];
			reflect[i][Z] = 2 * NdotL*norm[Z] - render->lights[i].direction[Z];
			normalized(reflect[i]);
		}
		else if (NdotL < 0 && NdotE < 0) {
			checkLights[i] = -1;
			reflect[i][X] = 2 * NdotL*(-norm[X]) - render->lights[i].direction[X];
			reflect[i][Y] = 2 * NdotL*(-norm[Y]) - render->lights[i].direction[Y];
			reflect[i][Z] = 2 * NdotL*(-norm[Z]) - render->lights[i].direction[Z];
			normalized(reflect[i]);
		}
		else {
			checkLights[i] = 0;
			continue;
		}
	}
	GzColor ambientIntensity;
	ambientIntensity[0] = render->Ka[0] * render->ambientlight.color[0];
	ambientIntensity[1] = render->Ka[1] * render->ambientlight.color[1]; // G
	ambientIntensity[2] = render->Ka[2] * render->ambientlight.color[2]; // B
	GzColor diffuseSum = { 0, 0, 0 };
	for (int i = 0; i < render->numlights; ++i) {
		if (checkLights[i] == 0) continue;
		if (checkLights[i] == 1) {
			// R
			diffuseSum[0] += render->lights[i].color[0] * dotProduct(norm, render->lights[i].direction);
			// G
			diffuseSum[1] += render->lights[i].color[1] * dotProduct(norm, render->lights[i].direction);
			// B
			diffuseSum[2] += render->lights[i].color[2] * dotProduct(norm, render->lights[i].direction);
		}
		else if (checkLights[i] == -1) {
			GzCoord negN = { -norm[X], -norm[Y], -norm[Z] };
			// R
			diffuseSum[0] += render->lights[i].color[0] * dotProduct(negN, render->lights[i].direction);
			// G
			diffuseSum[1] += render->lights[i].color[1] * dotProduct(negN, render->lights[i].direction);
			// B
			diffuseSum[2] += render->lights[i].color[2] * dotProduct(negN, render->lights[i].direction);
		}
	}
	GzColor diffuseIntensity;
	diffuseIntensity[0] = render->Kd[0] * diffuseSum[0]; // R
	diffuseIntensity[1] = render->Kd[1] * diffuseSum[1]; // G
	diffuseIntensity[2] = render->Kd[2] * diffuseSum[2]; // B

	GzColor specularSum = { 0, 0, 0 };
	for (int i = 0; i < render->numlights; ++i) {
		if (checkLights[i] == 0) continue;
		float RdotE = dotProduct(reflect[i], eye);
		if (RdotE < 0) RdotE = 0;
		if (RdotE > 1) RdotE = 1;
		// R
		specularSum[0] += render->lights[i].color[0] * pow(RdotE, render->spec);
		// G
		specularSum[1] += render->lights[i].color[1] * pow(RdotE, render->spec);
		// B
		specularSum[2] += render->lights[i].color[2] * pow(RdotE, render->spec);
	}
	GzColor specularIntensity;
	specularIntensity[0] = render->Ks[0] * specularSum[0]; // R
	specularIntensity[1] = render->Ks[1] * specularSum[1]; // G
	specularIntensity[2] = render->Ks[2] * specularSum[2]; // B
	
	color[0] = specularIntensity[0] + diffuseIntensity[0] + ambientIntensity[0];
	color[1] = specularIntensity[1] + diffuseIntensity[1] + ambientIntensity[1];
	color[2] = specularIntensity[2] + diffuseIntensity[2] + ambientIntensity[2];

	return GZ_SUCCESS;
}

int CalculateColorTexture(GzRender *render, GzColor color, GzCoord norm) {
	normalized(norm);
	GzCoord* reflect = (GzCoord*)malloc(sizeof(GzCoord) * render->numlights);
	float NdotL;
	int* checkLights = (int*)malloc(sizeof(int)*render->numlights);
	GzCoord eye;
	eye[X] = 0;
	eye[Y] = 0;
	eye[Z] = -1;
	normalized(eye);
	float NdotE = dotProduct(norm, eye);
	for (int i = 0; i < render->numlights; i++) {
		NdotL = dotProduct(norm, render->lights[i].direction);
		if (NdotL >= 0 && NdotE >= 0) {
			checkLights[i] = 1;
			reflect[i][X] = 2 * NdotL*norm[X] - render->lights[i].direction[X];
			reflect[i][Y] = 2 * NdotL*norm[Y] - render->lights[i].direction[Y];
			reflect[i][Z] = 2 * NdotL*norm[Z] - render->lights[i].direction[Z];
			normalized(reflect[i]);
		}
		else if (NdotL < 0 && NdotE < 0) {
			checkLights[i] = -1;
			reflect[i][X] = 2 * NdotL*(-norm[X]) - render->lights[i].direction[X];
			reflect[i][Y] = 2 * NdotL*(-norm[Y]) - render->lights[i].direction[Y];
			reflect[i][Z] = 2 * NdotL*(-norm[Z]) - render->lights[i].direction[Z];
			normalized(reflect[i]);
		}
		else {
			checkLights[i] = 0;
			continue;
		}
	}

	GzColor diffuseSum = { 0, 0, 0 };
	for (int i = 0; i < render->numlights; ++i) {
		if (checkLights[i] == 0) continue;
		if (checkLights[i] == 1) {
			// R
			diffuseSum[0] += render->lights[i].color[0] * dotProduct(norm, render->lights[i].direction);
			// G
			diffuseSum[1] += render->lights[i].color[1] * dotProduct(norm, render->lights[i].direction);
			// B
			diffuseSum[2] += render->lights[i].color[2] * dotProduct(norm, render->lights[i].direction);
		}
		else if (checkLights[i] == -1) {
			GzCoord negN = { -norm[X], -norm[Y], -norm[Z] };
			// R
			diffuseSum[0] += render->lights[i].color[0] * dotProduct(negN, render->lights[i].direction);
			// G
			diffuseSum[1] += render->lights[i].color[1] * dotProduct(negN, render->lights[i].direction);
			// B
			diffuseSum[2] += render->lights[i].color[2] * dotProduct(negN, render->lights[i].direction);
		}
	}

	GzColor specularSum = { 0, 0, 0 };
	for (int i = 0; i < render->numlights; ++i) {
		if (checkLights[i] == 0) continue;
		float RdotE = dotProduct(reflect[i], eye);
		if (RdotE < 0) RdotE = 0;
		if (RdotE > 1) RdotE = 1;
		// R
		specularSum[0] += render->lights[i].color[0] * pow(RdotE, render->spec);
		// G
		specularSum[1] += render->lights[i].color[1] * pow(RdotE, render->spec);
		// B
		specularSum[2] += render->lights[i].color[2] * pow(RdotE, render->spec);
	}

	color[0] = specularSum[0] + diffuseSum[0];
	color[1] = specularSum[1] + diffuseSum[1];
	color[2] = specularSum[2] + diffuseSum[2];

	return GZ_SUCCESS;
}