#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"

float *image_h = NULL;
int xs_h, ys_h;
int reset_height = 1;
extern float *heightMap;
extern int normXs, normYs;

/* Image texture function */
float height_fun(float u, float v, float normal)
{
	//unsigned char		pixel[3];
	//unsigned char     dummy;
	//char  		foo[8];
	//int   		i, j;
	//FILE			*fd;

	//if (reset_height) {          /* open and load texture file */
	//	fd = fopen("rock_h.ppm", "rb");
	//	if (fd == NULL) {
	//		fprintf(stderr, "normal file not found\n");
	//		exit(-1);
	//	}
	//	fscanf(fd, "%s %d %d %c", foo, &xs_h, &ys_h, &dummy);
	//	image_h = (GzColor*)malloc(sizeof(GzColor)*(xs_h + 1)*(ys_h + 1));
	//	if (image_h == NULL) {
	//		fprintf(stderr, "malloc for texture image failed\n");
	//		exit(-1);
	//	}

	//	for (i = 0; i < xs_h*ys_h; i++) {	/* create array of GzColor values */
	//		fread(pixel, sizeof(pixel), 1, fd);
	//		image_h[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
	//		image_h[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
	//		image_h[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
	//	}

	//	reset_height = 0;          /* init is done */
	//	fclose(fd);
	//}

	if (reset_height) {
		xs_h = normXs;
		ys_h = normYs;
		image_h = (float*)malloc(sizeof(float)*(xs_h + 1)*(ys_h + 1));
		for (int i = 0; i < xs_h*ys_h; i++) {	
			/* create array of GzColor values */
			image_h[i] = (heightMap[i]) * (1.0 / 255.0);
		}
		reset_height = 0;

#ifdef USEEXTERNHEIGHT
		unsigned char		pixel[3];
		unsigned char     dummy;
		char  		foo[8];
		int   		i, j;
		FILE			*fd;

		if (reset_height) {          /* open and load texture file */
			fd = fopen(HEIGHTMAPFILE, "rb");
			if (fd == NULL) {
				fprintf(stderr, "normal file not found\n");
				exit(-1);
			}
			fscanf(fd, "%s %d %d %c", foo, &xs_h, &ys_h, &dummy);
			image_h = (float*)malloc(sizeof(float)*(xs_h + 1)*(ys_h + 1));
			if (image_h == NULL) {
				fprintf(stderr, "malloc for texture image failed\n");
				exit(-1);
			}

			for (i = 0; i < xs_h*ys_h; i++) {	/* create array of GzColor values */
				fread(pixel, sizeof(pixel), 1, fd);
				image_h[i] = (float)((int)pixel[RED]) * (1.0 / 255.0);
			}

			reset_height = 0;          /* init is done */
			fclose(fd);
		}
#endif
	}


	/* bounds-test u,v to make sure nothing will overflow image array bounds */
	if (u > 1)  u = 1;
	if (u < 0)  u = 0;
	if (v > 1)  v = 1;
	if (v < 0)  v = 0;

	/* determine texture cell corner values and perform bilinear interpolation */
	/* set color to interpolated GzColor value and return */
	float px, py;
	float s, t;
	px = u * (xs_h - 1);
	py = v * (ys_h - 1);
	s = px - floor(px);
	t = py - floor(py);

	//GzCoord normalA, normalB, normalC, normalD;
	////color at A
	//for (i = 0; i < 3; i++)  normalA[i] = image_h[xs_h * (int)floor(py) + (int)floor(px)][i];
	////color at B
	//for (i = 0; i < 3; i++)  normalB[i] = image_h[xs_h * (int)floor(py) + (int)ceil(px)][i];
	////color at C
	//for (i = 0; i < 3; i++)  normalC[i] = image_h[xs_h * (int)ceil(py) + (int)ceil(px)][i];
	////color at D
	//for (i = 0; i < 3; i++)  normalD[i] = image_h[xs_h * (int)ceil(py) + (int)floor(px)][i];
	////color at P
	//for (i = 0; i < 3; i++)  normal[i] = s * t * normalC[i] + (1 - s) * t * normalD[i] + s * (1 - t) * normalB[i] + (1 - s) * (1 - t) * normalA[i];
	float normalA, normalB, normalC, normalD;
	//color at A
	normalA = image_h[xs_h * (int)floor(py) + (int)floor(px)];
	//color at B
	normalB = image_h[xs_h * (int)floor(py) + (int)ceil(px)];
	//color at C
	normalC = image_h[xs_h * (int)ceil(py) + (int)ceil(px)];
	//color at D
	normalD = image_h[xs_h * (int)ceil(py) + (int)floor(px)];
	//color at P
	normal = s * t * normalC + (1 - s) * t * normalD + s * (1 - t) * normalB + (1 - s) * (1 - t) * normalA;
	return normal;
	//return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeHeightTexture()
{
	if (image_h != NULL)
		free(image_h);
	return GZ_SUCCESS;
}

