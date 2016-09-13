/*   CS580 HW1 display functions to be completed   */

#include   "stdafx.h"  
#include	"Gz.h"
#include	"disp.h"
#include	<fstream>
#include	<iostream>
#include	<string>
using namespace std;

int GzNewFrameBuffer(char** framebuffer, int width, int height)
{
/* HW1.1 create a framebuffer for MS Windows display:
 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
 -- pass back pointer 
 */
	if (width > MAXXRES || height > MAXYRES)
	{
		return GZ_FAILURE;
	}

	if ((*framebuffer = (char*)malloc(width * height * 3)) == NULL)
	{
		AfxMessageBox(_T("Frame width or height is too large.\n"));
		return GZ_FAILURE;
	}

	return GZ_SUCCESS;
}

int GzNewDisplay(GzDisplay	**display, int xRes, int yRes)
{
/* HW1.2 create a display:
  -- allocate memory for indicated resolution
  -- pass back pointer to GzDisplay object in display
*/
	if (xRes > MAXXRES || yRes > MAXYRES)
	{
		AfxMessageBox(_T("Diaplay size too large.\n"));
		return GZ_FAILURE;
	}

	if (((*display = (GzDisplay*)malloc(sizeof(GzDisplay))) == NULL) ||
		((*display)->fbuf = (GzPixel*)malloc(sizeof(GzPixel) * xRes * yRes)) == NULL)
	{
		return GZ_FAILURE;
	}
	
	(*display)->xres = xRes;
	(*display)->yres = yRes;
	return GZ_SUCCESS;
}

int GzFreeDisplay(GzDisplay	*display)
{
/* HW1.3 clean up, free memory */
	free(display->fbuf);
	free(display);
	return GZ_SUCCESS;
}


int GzGetDisplayParams(GzDisplay *display, int *xRes, int *yRes)
{
/* HW1.4 pass back values for a display */
	*xRes = display->xres;
	*yRes = display->yres;
	return GZ_SUCCESS;
}


int GzInitDisplay(GzDisplay	*display)
{
	//OutputDebugStringA("My output string.");
/* HW1.5 set everything to some default values - start a new frame */
	for (int i = 0; i < display->xres * display->yres; i++)
	{
		display->fbuf[i].alpha = 1;
		display->fbuf[i].blue = 100;
		display->fbuf[i].green = 100;
		display->fbuf[i].z = MAXINT;
		display->fbuf[i].red = 100;
		
	}
	return GZ_SUCCESS;
}

int clamp(int x, int min, int max)
{
	if(x < min){
		x = min;
	}
	else if (x > max) {
		x = max;
	}
	return x;
}

int GzPutDisplay(GzDisplay *display, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
/* HW1.6 write pixel values into the display */
	if (i >= 0 && i < display->xres && j >= 0 && j < display->yres) {
		int index = display->xres * j + i;
		display->fbuf[index].alpha = a;
		display->fbuf[index].blue = clamp(b,0,4095) >> 4;
		display->fbuf[index].green = clamp(g,0,4095) >> 4;
		display->fbuf[index].red = clamp(r, 0,4095) >> 4;
		display->fbuf[index].z = z;
	}
	return GZ_SUCCESS;
}

int GzGetDisplay(GzDisplay *display, int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
/* HW1.7 pass back a pixel value to the display */
	*r = display->fbuf[j * display->xres + i].red;
	*g = display->fbuf[j * display->xres + i].green;
	*b = display->fbuf[j * display->xres + i].blue;
	*a = display->fbuf[j * display->xres + i].alpha;
	*z = display->fbuf[j * display->xres + i].z;
	return GZ_SUCCESS;
}

int GzFlushDisplay2File(FILE* outfile, GzDisplay *display)
{

/* HW1.8 write pixels to ppm file -- "P6 %d %d 255\r" */
	fprintf(outfile, "P6 %d %d 255\r", display->xres, display->yres);
	for (int i = 0; i < display->xres * display->yres; i++) {
		fprintf(outfile, "%c%c%c", display->fbuf[i].red, display->fbuf[i].green, display->fbuf[i].blue);
	}
	return GZ_SUCCESS;
}

int GzFlushDisplay2FrameBuffer(char* framebuffer, GzDisplay *display)
{

/* HW1.9 write pixels to framebuffer: 
	- put the pixels into the frame buffer
	- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red 
	- NOT red, green, and blue !!!
*/
	for (int i = 0; i < display->xres * display->yres; i++)
	{
		framebuffer[3 * i] = display->fbuf[i].blue;
		framebuffer[3 * i + 1] = display->fbuf[i].green;
		framebuffer[3 * i + 2] = display->fbuf[i].red;
	}
	
	return GZ_SUCCESS;
}