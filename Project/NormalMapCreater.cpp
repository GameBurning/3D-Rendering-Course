/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"


#define OUTPUTGRAY "output_heightMap.ppm"
#define OUTPUTNORM "output_normalMap.ppm"
GzColor	*normMap = NULL;
float *heightMap = NULL;
int normXs, normYs;
int normReset = 1;
float strength = 1;

int createGreyMap();
int createNormalMap();
float intensity(unsigned char* color, float s);
int pos(int x, int y);
void outputNormalMap();
void outputGrayScaleMap();

/* Image texture function */
int NormalMapCreater() {
	createGreyMap();
	createNormalMap();
	return GZ_SUCCESS;
}

int createGreyMap()
{
	unsigned char		pixel[3];
	unsigned char     dummy;
	char  		foo[8];
	int   		i, j;
	FILE			*fd;

	if (normReset) {          /* open and load texture file */
		fd = fopen(TEXTUREFILE, "rb");
		if (fd == NULL) {
			fprintf(stderr, "texture file not found\n");
			exit(-1);
		}
		fscanf(fd, "%s %d %d %c", foo, &normXs, &normYs, &dummy);
		normMap = (GzColor*)malloc(sizeof(GzColor)*(normXs + 1)*(normYs + 1));
		heightMap = (float*)malloc(sizeof(float)*(normXs + 1)*(normYs + 1));
		if ((normMap == NULL) || (heightMap == NULL)) {
			fprintf(stderr, "malloc for normMap or heightMap failed\n");
			exit(-1);
		}

		for (i = 0; i < normXs*normYs; i++) {	/* create array of GzColor values */
			fread(pixel, sizeof(pixel), 1, fd);
			heightMap[i] = intensity(pixel,strength);
		}

		normReset = 0;          /* init is done */
		fclose(fd);
	}
	outputGrayScaleMap();
	return GZ_SUCCESS;
}

int createNormalMap() {
	
	float topLeft, top, topRight, right, bottomRight, bottom, bottomLeft, left, dX, dY, dZ;
	for (int y = 0; y < normYs; y++) {
		for (int x = 0; x < normXs; x++ ) {
			topLeft = heightMap[pos(x - 1, y - 1)];
			top = heightMap[pos(x, y - 1)];
			topRight = heightMap[pos(x + 1, y - 1)];
			right = heightMap[pos(x + 1, y)];
			bottomRight = heightMap[pos(x + 1, y + 1)];
			bottom = heightMap[pos(x, y + 1)];
			bottomLeft = heightMap[pos(x - 1, y + 1)];
			left = heightMap[pos(x - 1, y)];
			// sobel operator
			// position.      Gx.            Gy
			// 1 2 3     |-1. 0. 1.|   |-1. -2. -1.|
			// 4 5 6     |-2. 0. 2.|   | 0.  0.  0.|
			// 7 8 9     |-1. 0. 1.|   | 1.  2.  1.|
			dX = (topRight + (2.0 * right) + bottomRight) - (topLeft + (2.0 * left) + bottomLeft);
			dY = (bottomLeft + (2.0 * bottom) + bottomRight) - (topLeft + (2.0 * top) + topRight);
			dZ = 255 / strength;
			dX = ((dX + 255.0) / 2.0);
			dY = ((dY + 255.0) / 2.0);
			dZ = ((dZ + 255.0) / 2.0);
			if (dX > 255)
				dX = 255;
			else if (dX < 0)
				dX = 0;
			if (dY > 255)
				dY = 255;
			else if (dY < 0)
				dY = 0;
			if (dZ > 255)
				dZ = 255;
			else if (dZ < 0)
				dZ = 0;
			normMap[pos(x, y)][RED] = dX;
			normMap[pos(x, y)][GREEN] = dY;
			normMap[pos(x, y)][BLUE] = dZ;
		}
	}
	outputNormalMap();
	return GZ_SUCCESS;
}

void outputNormalMap() {
	FILE *outfile;
	if ((outfile = fopen(OUTPUTNORM, "wb")) == NULL)
	{
		AfxMessageBox("The output file was not opened\n");
	}
	fprintf(outfile, "P6 %d %d 255\r", normXs, normYs);
	for (int i = 0; i < normXs * normYs; i++) {
		
		fprintf(outfile, "%c%c%c", (short)normMap[i][RED], (short)normMap[i][GREEN], (short)normMap[i][BLUE]);
	}
	fclose(outfile);
}

void outputGrayScaleMap() {
	FILE *outfile;
	if ((outfile = fopen(OUTPUTGRAY, "wb")) == NULL)
	{
		AfxMessageBox("The output file was not opened\n");
	}
	fprintf(outfile, "P6 %d %d 255\r", normXs, normYs);
	for (int i = 0; i < normXs * normYs; i++) {
		
		fprintf(outfile, "%c%c%c", (short)heightMap[i], (short)heightMap[i], (short)heightMap[i]);
	}
	fclose(outfile);
}

float intensity(unsigned char* color, float s) {
	float c = (float)(0.3*(float)((int)color[RED]) + 0.59*(float)((int)color[GREEN]) + 0.11*(float)((int)color[BLUE]))*(float)(s);
	if (c > 255)
		c = 255;
	else if (c < 0)
		c = 0;
	return c;
}

int pos(int x, int y) {
	int xx = x;
	int yy = y;
	if ((xx - 1) < 0) {
		xx = 0;
	}
	else if (xx + 1 >= normXs) {
		xx = normXs - 1;
	}
	if ((yy - 1) < 0) {
		yy = 0;
	}
	else if (yy + 1 >= normYs) {
		yy = normYs - 1;
	}
	return (xx + (yy*normXs));
}

/* Free texture memory */
int GzFreeMap()
{
	if (normMap != NULL)
		free(normMap);
	if (heightMap != NULL)
		free(heightMap);
	return GZ_SUCCESS;
}

