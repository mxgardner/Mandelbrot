
#include "bitmap.h"

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>

int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );
void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max );

typedef struct {
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    int image_width;
    int image_height;
    int max_iterations;
    int start_row;
    int end_row;
    struct bitmap *bm; 
} thread_data_t;

void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}

void* compute_band(void* arg) {
    thread_data_t* data = (thread_data_t*) arg;

    double xmin = data->xmin;
    double xmax = data->xmax;
    double ymin = data->ymin;
    double ymax = data->ymax;
    int image_width = data->image_width;
    int image_height = data->image_height;
    int max_iterations = data->max_iterations;
    int start_row = data->start_row;
    int end_row = data->end_row;
    struct bitmap* bm = data->bm;

    // loop over the rows in this thread
    for (int j = start_row; j < end_row; j++) {
        for (int i = 0; i < image_width; i++) {
            // calculate the Mandelbrot point and store the result in the bitmap
            double x = xmin + i * (xmax - xmin) / image_width;
            double y = ymin + j * (ymax - ymin) / image_height;
            int iters = iterations_at_point(x, y, max_iterations);
            bitmap_set(bm, i, j, iters);
        }
    }

    return NULL;
}

int main( int argc, char *argv[] )
{
	char c;

	// These are the default configuration values used
	// if no command line arguments are given.

	const char *outfile = "mandel.bmp";
	double xcenter = 0;
	double ycenter = 0;
	double scale = 4;
	int    image_width = 500;
	int    image_height = 500;
	int    max = 1000;
	int num_threads = 1;

	// For each command line argument given,
	// override the appropriate configuration value.

	while((c = getopt(argc,argv,"x:y:s:W:H:m:o:h"))!=-1) {
		switch(c) {
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				scale = atof(optarg);
				break;
			case 'W':
				image_width = atoi(optarg);
				break;
			case 'H':
				image_height = atoi(optarg);
				break;
			case 'm':
				max = atoi(optarg);
				break;
			case 'o':
				outfile = optarg;
				break;
			case 'h':
				show_help();
				exit(1);
				break;
			case 'n':
            	num_threads = atoi(optarg);
            	break;
		}
	}

	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d outfile=%s\n",xcenter,ycenter,scale,max,outfile);

	// Create a bitmap of the appropriate size.
	struct bitmap *bm = bitmap_create(image_width,image_height);

	// Fill it with a dark blue, for debugging
	bitmap_reset(bm,MAKE_RGBA(0,0,255,0));

	// Compute the Mandelbrot image
	compute_image(bm,xcenter-scale,xcenter+scale,ycenter-scale,ycenter+scale,max);

	// Calculate the image bounds based on center and scale
	double xmin = xcenter - scale;
	double xmax = xcenter + scale;
	double ymin = ycenter - scale;
	double ymax = ycenter + scale;

	  // Create and start threads
    pthread_t threads[num_threads];
    thread_data_t thread_data[num_threads];

    // Divide the image into bands and assign each band to a thread
    int rows_per_thread = image_height / num_threads;
    for (int t = 0; t < num_threads; t++) {
        thread_data[t].xmin = xmin;
        thread_data[t].xmax = xmax;
        thread_data[t].ymin = ymin;
        thread_data[t].ymax = ymax;
        thread_data[t].image_width = image_width;
        thread_data[t].image_height = image_height;
        thread_data[t].max_iterations = max;
        thread_data[t].bm = bm;
        thread_data[t].start_row = t * rows_per_thread;
        thread_data[t].end_row = (t == num_threads - 1) ? image_height : (thread_data[t].start_row + rows_per_thread);

        // Create the thread
        if (pthread_create(&threads[t], NULL, compute_band, (void*)&thread_data[t])) {
            fprintf(stderr, "Error creating thread %d\n", t);
            return 1;
        }
    }

    // Wait for all threads to finish
    for (int t = 0; t < num_threads; t++) {
        if (pthread_join(threads[t], NULL)) {
            fprintf(stderr, "Error joining thread %d\n", t);
            return 1;
        }
    }

	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}

	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max )
{
	int i,j;

	int width = bitmap_width(bm);
	int height = bitmap_height(bm);

	// For every pixel in the image...

	for(j=0;j<height;j++) {

		for(i=0;i<width;i++) {

			// Determine the point in x,y space for that pixel.
			double x = xmin + i*(xmax-xmin)/width;
			double y = ymin + j*(ymax-ymin)/height;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x,y,max);

			// Set the pixel in the bitmap.
			bitmap_set(bm,i,j,iters);
		}
	}
}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point( double x, double y, int max )
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while( (x*x + y*y <= 4) && iter < max ) {

		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter,max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color( int i, int max )
{
	int gray = 255*i/max;
	return MAKE_RGBA(0,0,gray,0);
}
