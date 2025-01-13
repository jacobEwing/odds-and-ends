#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <SDL.h>
#include <SDL_image.h>
/* include local files */
#include "gfx_functions.h"

/***** global variables *****/
int x_resolution, y_resolution, screen_depth;
SDL_Surface *screen_image;
void *free_pointer[MAXFREEPOINTERS];
int numfreepointers;


/***
These next three functions are here to simplyfy memory allocation and management.
gfx_allocate_memory allocates it, returns a pointer to it, and puts it in a list
of pointers that get freed at the end of the program.
gfx_free_memory is called at the end and handles freeing everything in the list.
gfx_add_freepointer takes a pointer to already allocated memory and shoves it in
the list.
***/
void *gfx_allocate_memory(size_t size){
	void * returnval = NULL;
	if(numfreepointers < MAXFREEPOINTERS){
		returnval = malloc(size);
		if(returnval){
			free_pointer[numfreepointers++] = returnval;
		}
	}
	return returnval;
}

int gfx_add_freepointer(void *ptr){
	int returnval = 1;
	if(numfreepointers < MAXFREEPOINTERS){
		free_pointer[numfreepointers++] = ptr;
		returnval = 0;
	}
	return returnval;
}

void gfx_free_memory(void){
	int n;
	for(n = 0; n < numfreepointers; n++){
		free(free_pointer[n]);
	}
	numfreepointers = 0;
}

/***** calculates the x size and y size of the polygon defined bu "*points"
and assigns them to "*xsize" and "*ysize".  It assumes that the polygon
has two or more dimensions.
*****/
void gfx_poly_area(int *points, int numpoints, int *x1, int *y1, int *x2, int *y2){
	int n, minx, miny, maxx, maxy;
	minx = maxx = points[0];
	miny = maxy = points[numpoints];
	for(n = 1; n < numpoints; n++){
		if(points[n] < minx) minx = points[n];
		if(points[n] > maxx) maxx = points[n];
		if(points[n + numpoints] < miny) miny = points[n + numpoints];
		if(points[n + numpoints] > maxy) maxy = points[n + numpoints];
	}
	*x1 = minx;
	*x2 = maxx;
	*y1 = miny;
	*y2 = maxy;
}

/*****
This function takes a polygon on a source image, a polygon on a destination
image, and copies the former to the latter.  The polygons do not need to be the
same shape, but they do need to have the same number of points.  Also, they must
be convex.  If the target polygon is not convex, the wrong shape may be drawn.
If the source polygon is not convex, you get some really cool-looking warped
effects in the image.  The locations at which the points are drawn are relative
the corners of the image, which allows for any transformation you want on the
image simply by transforming the corners of it.  The corners must be listed in
a cyclic manner.
*****/
void gfx_draw_texture(SDL_Surface *target, int *destpoint, int *sourcepoint, int numpoints, SDL_Surface *source){
	unsigned char done = 0;
	int minpoint, current_pointnum[2], n;
	int t_deltax[2], t_deltay[2];
	int t_sgndeltax[2], t_sgndeltay[2];
	int t_deltaxabs[2], t_deltayabs[2], t_deltax_tally[2];
	int t_lastcornerx[2], t_lastcornery[2], t_nextcornerx[2], t_nextcornery[2];
	int t_borderx[2], t_bordery[2];
	int s_deltax[2], s_deltay[2];
	int s_sgndeltax[2], s_sgndeltay[2];
	int s_deltaxabs[2], s_deltayabs[2], s_deltax_tally[2], s_deltay_tally[2];
	int s_lastcornerx[2], s_lastcornery[2], s_nextcornerx[2], s_nextcornery[2];
	int s_borderx[2], s_bordery[2];
	int minx, miny, maxx, maxy;
	SDL_Surface *buffer;

	/*** To allow alpha to work properly, we need to add this buffer.
	This should also have the benefit of speeding things up when we are
	calling this functino to draw directly to the screen. ***/
	gfx_poly_area(destpoint, numpoints, &minx, &miny, &maxx, &maxy);
	if(source->format->Amask != 0){
		buffer = SDL_CreateRGBSurface(SDL_SWSURFACE, maxx - minx + 1, maxy - miny + 1, 32, RMASK, GMASK, BMASK, AMASK);
	}else{
		buffer = SDL_CreateRGBSurface(SDL_SWSURFACE, maxx - minx + 1, maxy - miny + 1, 32, BMASK, GMASK, RMASK, 0);
	}

	Slock(source);
	Slock(target);

	if(numpoints > 2){
		/***** find the starting point ******/
		minpoint = 0;
		for(n = 1; n < numpoints; n++){
			if(destpoint[n + numpoints] < destpoint[minpoint + numpoints])
				minpoint = n;
		}

		/******** initialize the starting points *********/
		for(n = 0; n < 2; n++){
			t_lastcornerx[n] = destpoint[minpoint];
			t_lastcornery[n] = destpoint[minpoint + numpoints];
			s_lastcornerx[n] = sourcepoint[minpoint];
			s_lastcornery[n] = sourcepoint[minpoint + numpoints];

			current_pointnum[n] = minpoint - (n << 1) + 1;
			if(current_pointnum[n] < 0) current_pointnum[n] += numpoints;
			else current_pointnum[n] %= numpoints;

			t_nextcornerx[n] = destpoint[current_pointnum[n]];
			t_nextcornery[n] = destpoint[current_pointnum[n] + numpoints];
			s_nextcornerx[n] = sourcepoint[current_pointnum[n]];
			s_nextcornery[n] = sourcepoint[current_pointnum[n] + numpoints];

			t_deltax[n] = t_nextcornerx[n] - t_lastcornerx[n];
			t_deltay[n] = t_nextcornery[n] - t_lastcornery[n];
			s_deltax[n] = s_nextcornerx[n] - s_lastcornerx[n];
			s_deltay[n] = s_nextcornery[n] - s_lastcornery[n];
			t_deltaxabs[n] = abs(t_deltax[n]);
			t_deltayabs[n] = abs(t_deltay[n]);
			s_deltaxabs[n] = abs(s_deltax[n]);
			s_deltayabs[n] = abs(s_deltay[n]);
			t_sgndeltax[n] = sgn(t_deltax[n]);
			t_sgndeltay[n] = sgn(t_deltay[n]);
			s_sgndeltax[n] = sgn(s_deltax[n]);
			s_sgndeltay[n] = sgn(s_deltay[n]);
			t_deltax_tally[n] = t_deltayabs[n] >> 1;
			s_deltax_tally[n] = 0;//t_deltayabs[n] >> 1;
			s_deltay_tally[n] = 0;//t_deltayabs[n] >> 1;
			t_borderx[n] = t_lastcornerx[n];
			t_bordery[n] = t_lastcornery[n];
			s_borderx[n] = s_lastcornerx[n];
			s_bordery[n] = s_lastcornery[n];
		}
		/******** draw the polygon *********/
		while(!done){
			gfx_copy_dline_to_hline(buffer, t_borderx[0] - minx, t_borderx[1] - minx, t_bordery[0] - miny, source,
						s_borderx[0], s_bordery[0], s_borderx[1], s_bordery[1]);
			for(n = 0; n < 2; n++){
				if(t_bordery[n] == t_nextcornery[n]){
					if(t_nextcornery[0] == t_nextcornery[1] && t_nextcornerx[0] == t_nextcornerx[1]){
						done = 1;
					}else{
						t_lastcornerx[n] = t_nextcornerx[n];
						t_lastcornery[n] = t_nextcornery[n];
						s_lastcornerx[n] = s_nextcornerx[n];
						s_lastcornery[n] = s_nextcornery[n];

						current_pointnum[n] = current_pointnum[n] - (n << 1) + 1;
						if(current_pointnum[n] < 0) current_pointnum[n] += numpoints;
						else current_pointnum[n] %= numpoints;

						t_nextcornerx[n] = destpoint[current_pointnum[n]];
						t_nextcornery[n] = destpoint[current_pointnum[n] + numpoints];
						s_nextcornerx[n] = sourcepoint[current_pointnum[n]];
						s_nextcornery[n] = sourcepoint[current_pointnum[n] + numpoints];

						t_deltax[n] = t_nextcornerx[n] - t_lastcornerx[n];
						t_deltay[n] = t_nextcornery[n] - t_lastcornery[n];
						s_deltax[n] = s_nextcornerx[n] - s_lastcornerx[n];
						s_deltay[n] = s_nextcornery[n] - s_lastcornery[n];

						t_deltaxabs[n] = abs(t_deltax[n]);
						t_deltayabs[n] = abs(t_deltay[n]);
						s_deltaxabs[n] = abs(s_deltax[n]);
						s_deltayabs[n] = abs(s_deltay[n]);

						t_sgndeltax[n] = sgn(t_deltax[n]);
						t_sgndeltay[n] = sgn(t_deltay[n]);
						s_sgndeltax[n] = sgn(s_deltax[n]);
						s_sgndeltay[n] = sgn(s_deltay[n]);

						t_deltax_tally[n] = t_deltayabs[n] >> 1;
						s_deltax_tally[n] = t_deltayabs[n] >> 1;
						s_deltay_tally[n] = t_deltayabs[n] >> 1;

						t_borderx[n] = t_lastcornerx[n];
						t_bordery[n] = t_lastcornery[n];
						s_borderx[n] = s_lastcornerx[n];
						s_bordery[n] = s_lastcornery[n];
					}
				}

				t_deltax_tally[n] += t_deltaxabs[n];
				while(t_deltax_tally[n] >= t_deltayabs[n]){
					t_deltax_tally[n] -= (t_deltayabs[n] + !t_deltayabs[n]);
					t_borderx[n] += t_sgndeltax[n];
				}
				s_deltax_tally[n] += s_deltaxabs[n];
				while(s_deltax_tally[n] >= t_deltayabs[n]){
					s_deltax_tally[n] -= (t_deltayabs[n] + !t_deltayabs[n]);
					s_borderx[n] += s_sgndeltax[n];
				}
				s_deltay_tally[n] += s_deltayabs[n];
				while(s_deltay_tally[n] >= t_deltayabs[n]){
					s_deltay_tally[n] -= (t_deltayabs[n] + !t_deltayabs[n]);
					s_bordery[n] += s_sgndeltay[n];
				}

				t_bordery[n] += t_sgndeltay[n];
			}
		}
	}
	gfx_blit(target, minx, miny, buffer);
	SDL_FreeSurface(buffer);
	Sulock(source);
	Sulock(target);
}

/***
This function is a very limited one that serves a single purpose which probably
will only ever be called from a single spot.  It takes the pixels along the line
(sx1, sy1)-(sx2, sy2) on the image "source" and copies them to the horizontal line
(drawx1, drawy)-(drawx2, drawy) on the image "target". It is called from the function
gfx_draw_textures and was separated for the sole purpose of simplifying the code
and making it a bit more humanly readable.
***/
void gfx_copy_dline_to_hline(SDL_Surface *target, int drawx1, int drawx2, int drawy, SDL_Surface *source, int sx1, int sy1, int sx2, int sy2){
	int deltadx, xtally = 0, ytally = 0, dtally1 = 0, dtally2 = 0;
	int drawx, readx, ready, sgndeltadx;
	int deltax, deltay, sgndeltax, sgndeltay, absdeltax, absdeltay;
	div_t division1, division2;

	deltadx = drawx2 - drawx1;
	if(deltadx < 0){
		deltadx *= -1;
		sgndeltadx = -1;
	}else{
		sgndeltadx = 1;
	}
	deltadx++;
	drawx = drawx1;
	deltax = sx2 - sx1;
	deltay = sy2 - sy1;
	absdeltax = abs(deltax);
	absdeltay = abs(deltay);
	sgndeltax = sgn(deltax);
	sgndeltay = sgn(deltay);

	/*** The order of this if structure is important.  The latter two
	conditions very much depend on the former two not being met. ***/
	if(deltadx >= absdeltax && deltadx >= absdeltay){
		readx = sx1;
		ready = sy1;
		while(drawx != drawx2){
			xtally += absdeltax;
			if(xtally > deltadx){
				xtally -= deltadx;
				readx += sgndeltax;
			}

			ytally += absdeltay;
			if(ytally > deltadx){
				ytally -= deltadx;
				ready += sgndeltay;
			}

			drawx += sgndeltadx;
			gfx_putrawpixel(target, drawx, drawy, gfx_getrawpixel(source, readx, ready));
		}
	}else if(absdeltax >= deltadx && absdeltay >= deltadx){
		division1 = div(absdeltax, deltadx);
		division2 = div(absdeltay, deltadx);
		division1.quot *= sgndeltax;
		division2.quot *= sgndeltay;
		readx = sx1 - division1.quot;
		ready = sy1 - division2.quot;
		while(drawx != drawx2){
			dtally2 += division2.rem;
			ready += division2.quot;
			if(dtally2 >= deltay){
				ready += sgndeltay;
				dtally2 -= deltadx;
			}

			dtally1 += division1.rem;
			readx += division1.quot;
			if(dtally1 >= deltax){
				readx += sgndeltax;
				dtally1 -= deltadx;
			}

			drawx += sgndeltadx;
			gfx_putrawpixel(target, drawx, drawy, gfx_getrawpixel(source, readx, ready));
		}
	}else if(absdeltax >= absdeltay){
		division1 = div(absdeltax, deltadx);
		division1.quot *= sgndeltax;
		readx = sx1 - division1.quot;
		ready = sy1;
		while(drawx != drawx2){
			ytally += absdeltay;
			if(ytally > deltadx){
				ytally -= deltadx;
				ready += sgndeltay;
			}

			dtally1 += division1.rem;
			readx += division1.quot;
			if(dtally1 >= deltax){
				readx += sgndeltax;
				dtally1 -= deltadx;
			}

			drawx += sgndeltadx;
			gfx_putrawpixel(target, drawx, drawy, gfx_getrawpixel(source, readx, ready));
		}
	}else if(absdeltay >= absdeltax){
		division1 = div(absdeltay, deltadx);
		division1.quot *= sgndeltay;
		readx = sx1;
		ready = sy1 - division1.quot;
		while(drawx != drawx2){
			xtally += absdeltax;
			if(xtally > deltadx){
				xtally -= deltadx;
				readx += sgndeltax;
			}

			dtally1 += division1.rem;
			ready += division1.quot;
			if(dtally1 >= deltay){
				ready += sgndeltay;
				dtally1 -= deltadx;
			}

			drawx += sgndeltadx;
			gfx_putrawpixel(target, drawx, drawy, gfx_getrawpixel(source, readx, ready));
		}
	}
}

/***
Draws the image "img" onto "target" at location (blitx, blity), and 
scaled to the dimensions of xblitsize by yblitsize.  This is done by using
a variation on Bresenham's line algorithm to handle the scaling.
***/
void gfx_scale_blit(SDL_Surface *target, int blitx, int blity, int xblitsize, int yblitsize, SDL_Surface *img){
	int real_xsize, real_ysize, xtally, ytally;
	int x1, y1, x2, y2, readx, ready, drawx, drawy;
	Uint8 R, G, B;

	real_xsize = img->w;
	real_ysize = img->h;
	if(xblitsize < 0) xblitsize *= -1;
	if(yblitsize < 0) yblitsize *= -1;
	x1 = blitx;
	y1 = blity;
	x2 = blitx + xblitsize - 1;
	y2 = blity + yblitsize - 1;
	ready = 0;
	ytally = 1;
	for(drawy = y1; drawy <= y2; drawy++){
		readx = 0;
		xtally = 1;
		for(drawx = x1; drawx <= x2; drawx++){
			gfx_getpixel(img, readx, ready, &R, &G, &B);
			gfx_putpixel(target, drawx, drawy, R, G, B);

			xtally += real_xsize;
			while(xtally > xblitsize && readx < real_xsize - 1){
				xtally -= xblitsize;
				readx ++;
			}
		}

		ytally += real_ysize;
		while(ytally > yblitsize && ready < real_ysize){
			ytally -= yblitsize;
			ready ++;
		}
	}
}

int gfx_init(int xres, int yres, int bdepth, Uint32 flags){
	int returnval = 0;

	x_resolution = xres;
	y_resolution = yres;
	screen_depth = bdepth;

	screen_image = SDL_SetVideoMode(x_resolution, y_resolution, screen_depth, flags);

	if(screen_image == NULL){
		fprintf(stderr, "Unable to set video: %s\n", SDL_GetError());
		returnval = 1;
	}

	return returnval;
}

SDL_Surface *gfx_getimage(SDL_Surface *src, int x1, int y1, int x2, int y2){
	SDL_Rect source, dest;
	SDL_Surface *returnval;

	if(x1 > x2){
		x1 ^= x2;
		x2 ^= x1;
		x1 ^= x2;
	}
	if(y1 > y2){
		y1 ^= y2;
		y2 ^= y1;
		y1 ^= y2;
	}
	
	source.x = x1;
	source.y = y1;
	dest.x = 0;
	dest.y = 0;
	source.w = dest.w = x2 - x1 + 1;
	source.h = dest.h = y2 - y1 + 1;

	returnval = SDL_CreateRGBSurface(SDL_SWSURFACE, source.w, source.h, 32, RMASK, GMASK, BMASK, AMASK);
	if(returnval != NULL)
		SDL_BlitSurface(src, &source, returnval, &dest);

	return returnval;
}

void gfx_blit(SDL_Surface *target, int cx, int cy, SDL_Surface *img){
	SDL_Rect dest;
	dest.x = cx;
	dest.y = cy;

	SDL_BlitSurface(img, NULL, target, &dest);
}

void gfx_blitarea(SDL_Surface *target, int cx, int cy, SDL_Surface *img, int rx1, int ry1, int rx2, int ry2){
	SDL_Rect dest, source;
	dest.x = cx;
	dest.y = cy;

	source.x = rx1;
	source.y = ry1;
	source.w = rx2 - rx1;
	source.h = ry2 - ry1;
	SDL_BlitSurface(img, &source, target, &dest);
}

/**** Blits an image onto the screen at the specified angle and scale.
All transformations are done relative to {rotx, roty} ****/
void gfx_rotblit(SDL_Surface *target, int x, int y, int rotx, int roty, float ang, float scale, SDL_Surface *img){
	float sine, cosine, x1, y1, x2, y2;
	int drawpoint[8], readpoint[8], w, h;

	/*** First, calculate the polygon to which we're drawing ***/
	sine = sin(ang);
	cosine = cos(ang);
	x1 = - rotx * scale;
	y1 = - roty * scale;
	x2 = x1 + img->w * scale;
	y2 = y1 + img->h * scale;
	
	drawpoint[0] = x + (int)(x1 * cosine - y1 * sine);
	drawpoint[4] = y + (int)(x1 * sine + y1 * cosine);

	drawpoint[1] = x + (int)(x2 * cosine - y1 * sine);
	drawpoint[5] = y + (int)(x2 * sine + y1 * cosine);

	drawpoint[2] = x + (int)(x2 * cosine - y2 * sine);
	drawpoint[6] = y + (int)(x2 * sine + y2 * cosine);

	drawpoint[3] = x + (int)(x1 * cosine - y2 * sine);
	drawpoint[7] = y + (int)(x1 * sine + y2 * cosine);

	/*** Now assign the polygon from which we're reading ***/
	w = img->w - 1;
	h = img->h - 1;
	readpoint[0] = 0;
	readpoint[4] = 0;
	readpoint[1] = w;
	readpoint[5] = 0;
	readpoint[2] = w;
	readpoint[6] = h;
	readpoint[3] = 0;
	readpoint[7] = h;

	/*** ready to roll, now call the function that actually does the real work ***/
	gfx_draw_texture(target, drawpoint, readpoint, 4, img);
}

/*** Calculates the smallest rectangle that surrounds the image once
the image has been scaled and rotated.  Saves the top left and bottom
right corners in the variables {minx, miny} and {maxx, maxy} ***/
void gfx_calcrotarea(int x, int y, int rotx, int roty, float ang, float scale, SDL_Surface *img, int *minx, int *miny, int *maxx, int *maxy){
	float sine, cosine, x1, y1, x2, y2;
	int drawpoint[8], n, xindex, yindex;

	/*** First, calculate the polygon to which we're drawing ***/
	sine = sin(ang);
	cosine = cos(ang);
	x1 = - rotx * scale;
	y1 = - roty * scale;
	x2 = x1 + img->w * scale;
	y2 = y1 + img->h * scale;
	
	drawpoint[0] = x + (int)(x1 * cosine - y1 * sine);
	drawpoint[4] = y + (int)(x1 * sine + y1 * cosine);

	drawpoint[1] = x + (int)(x2 * cosine - y1 * sine);
	drawpoint[5] = y + (int)(x2 * sine + y1 * cosine);

	drawpoint[2] = x + (int)(x2 * cosine - y2 * sine);
	drawpoint[6] = y + (int)(x2 * sine + y2 * cosine);

	drawpoint[3] = x + (int)(x1 * cosine - y2 * sine);
	drawpoint[7] = y + (int)(x1 * sine + y2 * cosine);

	*minx = drawpoint[0];
	*maxx = drawpoint[0];
	*miny = drawpoint[4];
	*maxy = drawpoint[4];

	for(n = 1; n < 4; n++){
		xindex = n;
		yindex = n + 4;
		if(drawpoint[xindex] < *minx) *minx = drawpoint[xindex];
		if(drawpoint[xindex] > *maxx) *maxx = drawpoint[xindex];
		if(drawpoint[yindex] < *miny) *miny = drawpoint[yindex];
		if(drawpoint[yindex] > *maxy) *maxy = drawpoint[yindex];
	}
}

/***
Reads the RGB values of the pixel at locaiton X, Y on image img, and assigns
them to the R, G and B variables.
***/
void gfx_getpixel(SDL_Surface *img, int x, int y, Uint8 *R, Uint8 *G, Uint8 *B){
	int bpp = img->format->BytesPerPixel;
	Uint8 *p = (Uint8 *)img->pixels + y * img->pitch + x * bpp;
	SDL_GetRGB(*p, img->format, R, G, B);
}

/*** returns the colour of the pixel at location X, Y on image img ***/
Uint32 gfx_getrawpixel(SDL_Surface *img, int x, int y){
	int bpp = img->format->BytesPerPixel;
	/* Here p is the address to the pixel we want to retrieve */
	Uint8 *p = (Uint8 *)img->pixels + y * img->pitch + x * bpp;

	switch(bpp) {
		case 1: return *p;

		case 2: return *(Uint16 *)p;

		case 3:
			if(SDL_BYTEORDER == SDL_BIG_ENDIAN)
				return p[0] | p[1] << 8 | p[2] << 16;
			else
				return p[0] << 16 | p[1] << 8 | p[2];

		case 4:
			return *(Uint32 *)p;

		default:
			return 0;       /* shouldn't happen, but avoids warnings */
	}
}

/**
Draws a pixel at location x, y on "surface" using the colour "pixel".
This is copied almost exactly from an example given in the SDL manual.
**/
void gfx_putrawpixel(SDL_Surface *surface, int x, int y, Uint32 pixel){
	int bpp = surface->format->BytesPerPixel;
	/* Here p is the address to the pixel we want to set */
	Uint8 *p = (Uint8 *)surface->pixels + y * surface->pitch + x * bpp;
	if(x >= 0 && x < surface->w && y >= 0 && y < surface->h){
		switch(bpp) {
			case 1:
				*p = pixel;
				break;

			case 2:
				*(Uint16 *)p = pixel;
				break;

			case 3:
				if(SDL_BYTEORDER == SDL_BIG_ENDIAN) {
					p[0] = (pixel >> 16) & 0xff;
					p[1] = (pixel >> 8) & 0xff;
					p[2] = pixel & 0xff;
				} else {
					p[0] = pixel & 0xff;
					p[1] = (pixel >> 8) & 0xff;
					p[2] = (pixel >> 16) & 0xff;
				}
				break;

			case 4:
				*(Uint32 *)p = pixel;
				break;
		}
	}
}

/***
Draws a pixel composed of the R, G and B colour bases at location x, y, on
the surface img.
***/
void gfx_putpixel(SDL_Surface *img, int x, int y, Uint8 R, Uint8 G, Uint8 B){
	Uint32 colour = SDL_MapRGB(img->format, R, G, B);
	Uint8 *bp8;
	Uint16 *bp16;
	Uint32 *bp32;
	if(x >= 0 && y >= 0 && x < img->w && y < img->h){
		switch(img->format->BytesPerPixel){
			case 1: // Assuming 8-bpp
				bp8 = (Uint8 *)img->pixels + y * img->pitch + x;
				*bp8 = colour;
				break;
			case 2: // Probably 15-bpp or 16-bpp
				bp16 = (Uint16 *)img->pixels + y * img->pitch / 2 + x;
				*bp16 = colour;
				break;
			case 3: // Slow 24-bpp mode, usually not used
				bp8 = (Uint8 *)img->pixels + y * img->pitch + x * 3;
				if(SDL_BYTEORDER == SDL_LIL_ENDIAN)
				{
					bp8[0] = colour;
					bp8[1] = colour >> 8;
					bp8[2] = colour >> 16;
				}else{
					bp8[2] = colour;
					bp8[1] = colour >> 8;
					bp8[0] = colour >> 16;
				}
				break;
			case 4: // Probably 32-bpp
				bp32 = (Uint32 *)img->pixels + y * (img->pitch >> 2) + x;
				*bp32 = colour;
				break;
		}
	}
}

void gfx_line(SDL_Surface *img, int x1, int y1, int x2, int y2, Uint8 R, Uint8 G, Uint8 B){
	//Uint32 colour = SDL_MapRGB(img->format, R, G, B);
	int n, deltax, deltay, sgndeltax, sgndeltay, deltaxabs, deltayabs, x, y, drawx, drawy;

	deltax = x2 - x1;
	deltay = y2 - y1;
	deltaxabs = abs(deltax);
	deltayabs = abs(deltay);
	sgndeltax = sgn(deltax);
	sgndeltay = sgn(deltay);
	x = deltayabs >> 1;
	y = deltaxabs >> 1;
	drawx = x1;
	drawy = y1;

	gfx_putpixel(img, drawx, drawy, R, G, B);

	if(deltaxabs >= deltayabs){
		for(n = 0; n < deltaxabs; n++){
			y += deltayabs;
			if(y >= deltaxabs){
				y -= deltaxabs;
				drawy += sgndeltay;
			}
			drawx += sgndeltax;
			gfx_putpixel(img, drawx, drawy, R, G, B);
		}
	}else{  
		for(n = 0; n < deltayabs; n++){
			x += deltaxabs;
			if(x >= deltayabs){
				x -= deltayabs;
				drawx += sgndeltax;
			}
			drawy += sgndeltay;
			gfx_putpixel(img, drawx, drawy, R, G, B);
		}
	}


}

void Slock(SDL_Surface *img){
	if(SDL_MUSTLOCK(img)){
		SDL_LockSurface(img);
/*		if(SDL_LockSurface(img) < 0){
			return;
		}
*/	}
}

void Sulock(SDL_Surface *img){
	if(SDL_MUSTLOCK(img))
		SDL_UnlockSurface(img);
}

int sgn(int v){
	if(v < 0) return -1;
	return v > 0;
}

void gfx_solid_box(SDL_Surface *img, int x1, int y1, int x2, int y2, Uint8 R, Uint8 G, Uint8 B, Uint8 A){
	int y, maxx, maxy;
	maxx = img->w;
	maxy = img->h;

	/*** account for limits ***/
	/* first, ensure that [x1, y1] is top left and [x2, y2] is bottom right. */
	if(x1 > x2) { x1 ^= x2; x2 ^= x1; x1 ^= x2; }
	if(y1 > y2) { y1 ^= y2; y2 ^= y1; y1 ^= y2; }

	/* Now make sure everythings's on the screen */
	if(x1 < maxx && x2 >= 0 && y1 < maxy && y2 >= 0){
		if(x1 < 0) x1 = 0;
		if(y1 < 0) y1 = 0;
		if(x2 >= maxx) x2 = maxx - 1;
		if(y2 >= maxy) y2 = maxy - 1;

		/*** everything is now limited to the screen, draw a box ***/
		for(y = y1; y <= y2; y++){
			gfx_hline(img, x1, x2, y, R, G, B, A);
		}
	}
	
}

/*** draws a horizontal line on the "img" surface.  Checks restrictions first ***/
void gfx_hline(SDL_Surface *img, int x1, int x2, int y, Uint8 R, Uint8 G, Uint8 B, Uint8 A){
	Uint32 colour = SDL_MapRGBA(img->format, R, G, B, A);
	Uint8 *bp8;
	Uint16 *bp16;
	Uint32 *bp32;
	int x;

	if(y >= 0 && y < img->h){
		if(x1 > x2){ /** make sure x1 is on the left and x2 on the right **/
			x1 ^= x2;
			x2 ^= x1;
			x1 ^= x2;
		}

		if(x2 >= 0 && x1 < img->w){ /** make sure the line is on the screen **/
			if(x2 > img->w) x2 = img->w - 1; /** clip any parts that are off the screen **/
			if(x1 < 0) x1 = 0;

			/** graphical limits are accounted for, now draw **/
			switch(img->format->BytesPerPixel){
				case 1: // Assuming 8-bpp
					bp8 = (Uint8 *)img->pixels + y * img->pitch + x1;
					for(x = x1; x <= x2; x++){
						*bp8 = colour;
						bp8++;
					}
					break;

				case 2: // Probably 15-bpp or 16-bpp
					bp16 = (Uint16 *)img->pixels + (y * img->pitch >> 1) + x1;
					for(x = x1; x <= x2; x++){
						*bp16 = colour;
						bp16++;
					}
					break;

				case 3: // Slow 24-bpp mode, usually not used
					bp8 = (Uint8 *)img->pixels + y * img->pitch + x1 * 3;
					if(SDL_BYTEORDER == SDL_LIL_ENDIAN)
					{
						for(x = x1; x <= x2; x++){
							bp8[0] = colour;
							bp8[1] = colour >> 8;
							bp8[2] = colour >> 16;
							bp8 += 3;
						}
					}else{
						for(x = x1; x <= x2; x++){
							bp8[2] = colour;
							bp8[1] = colour >> 8;
							bp8[0] = colour >> 16;
							bp8 += 3;
						}
					}
					break;

				case 4: // Probably 32-bpp
					bp32 = (Uint32 *)img->pixels + y * (img->pitch >> 2) + x1;
					for(x = x1; x <= x2; x++){
						*bp32 = colour;
						bp32 ++;
					}
					break;
			}
			
		}
	}
}

/***** draws a convex polygon *****/
void gfx_polygon(SDL_Surface *img, int *pointx, int *pointy, int numpoints, Uint8 R, Uint8 G, Uint8 B, Uint8 A){
	char done = 0;
	int minpoint, current_pointnum[2], n;
	int deltax[2], deltay[2];
	int sgndeltax[2], sgndeltay[2];
	int deltaxabs[2], deltayabs[2];
	int x[2], lastcornerx[2], lastcornery[2], nextcornerx[2], nextcornery[2];
	int drawx[2], drawy[2];

	if(numpoints > 2){
		/***** find the starting point ******/
		minpoint = 0;
		for(n = 1; n < numpoints; n++){
			if(pointy[n] < pointy[minpoint])
				minpoint = n;
		}

		/******** initialize the starting points *********/
		for(n = 0; n < 2; n++){
			lastcornerx[n] = pointx[minpoint];
			lastcornery[n] = pointy[minpoint];

			current_pointnum[n] = minpoint - (n << 1) + 1;
			if(current_pointnum[n] < 0) current_pointnum[n] += numpoints;
			else current_pointnum[n] %= numpoints;

			nextcornerx[n] = pointx[current_pointnum[n]];
			nextcornery[n] = pointy[current_pointnum[n]];

			deltax[n] = nextcornerx[n] - lastcornerx[n];
			deltay[n] = nextcornery[n] - lastcornery[n];
			deltaxabs[n] = abs(deltax[n]);
			deltayabs[n] = abs(deltay[n]);
			sgndeltax[n] = sgn(deltax[n]);
			sgndeltay[n] = sgn(deltay[n]);
			x[n] = deltayabs[n] >> 1;
			drawx[n] = lastcornerx[n];
			drawy[n] = lastcornery[n];
		}
		/******** draw the polygon *********/
		while(!done){

			gfx_hline(img, drawx[0], drawx[1], drawy[0], R, G, B, A);
			for(n = 0; n < 2; n++){
				if(drawy[n] == nextcornery[n]){
					if(nextcornery[0] == nextcornery[1] && nextcornerx[0] == nextcornerx[1]) done = 1;
					else{
						lastcornerx[n] = nextcornerx[n];
						lastcornery[n] = nextcornery[n];

						current_pointnum[n] = current_pointnum[n] - (n << 1) + 1;
						if(current_pointnum[n] < 0) current_pointnum[n] += numpoints;
						else current_pointnum[n] %= numpoints;

						nextcornerx[n] = pointx[current_pointnum[n]];
						nextcornery[n] = pointy[current_pointnum[n]];

						deltax[n] = nextcornerx[n] - lastcornerx[n];
						deltay[n] = nextcornery[n] - lastcornery[n];
						deltaxabs[n] = abs(deltax[n]);
						deltayabs[n] = abs(deltay[n]);
						sgndeltax[n] = sgn(deltax[n]);
						sgndeltay[n] = sgn(deltay[n]);
						x[n] = deltayabs[n] >> 1;
						drawx[n] = lastcornerx[n];
						drawy[n] = lastcornery[n];
					}
				}
				x[n] += deltaxabs[n];
				while(x[n] >= deltayabs[n]){
					x[n] -= (deltayabs[n] + !deltayabs[n]);
					drawx[n] += sgndeltax[n];
				}
				drawy[n] += sgndeltay[n];
			}
		}
	}
}

SDL_Surface *gfx_load_png(const char *filename){
	SDL_Surface *image;
	SDL_RWops *rwop;
	rwop = SDL_RWFromFile((const char *)filename, "rb");
	image = IMG_LoadPNG_RW(rwop);
	if(!image) {
		fprintf(stdout, "IMG_LoadPNG_RW: %s\n", IMG_GetError());
	}
	return image;
}
