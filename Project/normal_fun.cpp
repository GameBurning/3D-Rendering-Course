#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"

GzColor	*image_normal = NULL;
int xs_normal, ys_normal;
int reset_normal = 1;
extern GzColor *normMap;
extern int normXs, normYs;

extern int NormalMapCreater();



/* Image texture function */
int normal_fun(float u, float v, GzColor normal)
{
	unsigned char		pixel[3];
	unsigned char     dummy;
	char  		foo[8];
	int   		i, j;
	FILE			*fd;

	if (reset_normal) {          /* open and load texture file */

		//fd = fopen("rock_n.ppm", "rb");
		//if (fd == NULL) {
		//	fprintf(stderr, "normal file not found\n");
		//	exit(-1);
		//}
		//fscanf(fd, "%s %d %d %c", foo, &xs_normal, &ys_normal, &dummy);
		//image_normal = (GzColor*)malloc(sizeof(GzColor)*(xs_normal + 1)*(ys_normal + 1));
		//if (image_normal == NULL) {
		//	fprintf(stderr, "malloc for texture image failed\n");
		//	exit(-1);
		//}

		//for (i = 0; i < xs_normal*ys_normal; i++) {	/* create array of GzColor values */
		//	fread(pixel, sizeof(pixel), 1, fd);
		//	image_normal[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
		//	image_normal[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
		//	image_normal[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
		//}

		//reset_normal = 0;          /* init is done */
		//fclose(fd);


		NormalMapCreater();
		//image_normal = normMap;
		xs_normal = normXs;
		ys_normal = normYs;
		image_normal = (GzColor*)malloc(sizeof(GzColor)*(xs_normal + 1)*(ys_normal + 1));
		for (i = 0; i < xs_normal*ys_normal; i++) {	/* create array of GzColor values */
			image_normal[i][RED] = (float)(normMap[i][RED]) * (1.0 / 255.0);
			image_normal[i][GREEN] = (float)(normMap[i][GREEN]) * (1.0 / 255.0);
			image_normal[i][BLUE] = (float)(normMap[i][BLUE]) * (1.0 / 255.0);
		}
		reset_normal = 0;

#ifdef USEEXTERNNORMAL				
		fd = fopen(NORMALMAPFILE, "rb");
		if (fd == NULL) {
			fprintf(stderr, "normal file not found\n");
			exit(-1);
		}
		fscanf(fd, "%s %d %d %c", foo, &xs_normal, &ys_normal, &dummy);
		image_normal = (GzColor*)malloc(sizeof(GzColor)*(xs_normal + 1)*(ys_normal + 1));
		if (image_normal == NULL) {
			fprintf(stderr, "malloc for texture image failed\n");
			exit(-1);
		}

		for (i = 0; i < xs_normal*ys_normal; i++) {	/* create array of GzColor values */
			fread(pixel, sizeof(pixel), 1, fd);
			image_normal[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
			image_normal[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
			image_normal[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
		}

		reset_normal = 0;          /* init is done */
		fclose(fd);
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
	px = u * (xs_normal - 1);
	py = v * (ys_normal - 1);
	s = px - floor(px);
	t = py - floor(py);

	GzCoord normalA, normalB, normalC, normalD;
	//color at A
	for (i = 0; i < 3; i++)  normalA[i] = image_normal[xs_normal * (int)floor(py) + (int)floor(px)][i];
	//color at B
	for (i = 0; i < 3; i++)  normalB[i] = image_normal[xs_normal * (int)floor(py) + (int)ceil(px)][i];
	//color at C
	for (i = 0; i < 3; i++)  normalC[i] = image_normal[xs_normal * (int)ceil(py) + (int)ceil(px)][i];
	//color at D
	for (i = 0; i < 3; i++)  normalD[i] = image_normal[xs_normal * (int)ceil(py) + (int)floor(px)][i];
	//color at P
	for (i = 0; i < 3; i++)  normal[i] = s * t * normalC[i] + (1 - s) * t * normalD[i] + s * (1 - t) * normalB[i] + (1 - s) * (1 - t) * normalA[i];
	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeNormalTexture()
{
	if (image_normal != NULL)
		free(image_normal);
	return GZ_SUCCESS;
}

