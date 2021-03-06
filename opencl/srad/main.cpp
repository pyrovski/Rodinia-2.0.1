/*

  --------------------------------------------
  UPDATES
  ------------------------------------------------
  March 2006
  Rob Janiczek
  -- creation of prototype version

  March 2006
  Drew Gilliam
  -- rewriting of prototype version into current version
  -- got rid of multiple function calls, all code in a
  single function (for speed)
  -- code cleanup & commenting
  -- code optimization efforts

  April 2006
  Drew Gilliam
  -- added diffusion coefficent saturation on [0,1]

  December 2009
  Lukasz G. Szafaryn
  -- reading from image, command line inputs

  January 2010
  Lukasz G. Szafaryn
  -- comments

  October 2010
  W. Mike Wilkerson
  -- Ported to OpenCL from CUDA.


  ------------------------------------------------
  DESCRIPTION
  ------------------------------------------------
  The Heart Wall application tracks the movement of a mouse heart over a sequence of 104 609x590 ultrasound images to record response to the stimulus.
  In its initial stage, the program performs image processing operations on the first image to detect initial, partial shapes of inner and outer heart walls.
  These operations include: edge detection, SRAD despeckling (also part of Rodinia suite), morphological transformation and dilation. In order to reconstruct
  approximated full shapes of heart walls, the program generates ellipses that are superimposed over the image and sampled to mark points on the heart walls
  (Hough Search). In its final stage (Heart Wall Tracking presented here), program tracks movement of surfaces by detecting the movement of image areas under
  sample points as the shapes of the heart walls change throughout the sequence of images.

  SRAD is one of the first stages of the Heart Wall application. SRAD (Speckle Reducing Anisotropic Diffusion) is a diffusion method for ultrasonic and radar imaging
  applications based on partial differential equations (PDEs). It is used to remove locally correlated noise, known as speckles, without destroying important image
  features. SRAD consists of several pieces of work: image extraction, continuous iterations over the image (preparation, reduction, statistics, computation 1 and
  computation 2) and image compression. The sequential dependency between all of these stages requires synchronization after each stage (because each stage
  operates on the entire image).

  For more information, see "PAPERS" and "DOWNLOADS" below.


  ------------------------------------------------
  PAPERS:
  ------------------------------------------------
  L. G. Szafaryn, K. Skadron, and J. J. Saucerman. "Experiences Accelerating MATLAB Systems Biology Applications." In Proceedings of the Workshop on Biomedicine
  in Computing: Systems, Architectures, and Circuits (BiC) 2009, in conjunction with the 36th IEEE/ACM International Symposium on Computer Architecture (ISCA),
  June 2009. <http://www.cs.virginia.edu/~skadron/Papers/BiC09.pdf>

  Y. Yu, S. Acton, Speckle reducing anisotropic diffusion, IEEE Transactions on Image Processing 11(11)(2002) 1260-1270.
  <http://people.virginia.edu/~sc5nf/01097762.pdf>


  ------------------------------------------------
  DOWNLOADS:
  ------------------------------------------------
  Rodinia Benchmark Suite <https://www.cs.virginia.edu/~skadron/wiki/rodinia/index.php/Main_Page>


  ------------------------------------------------
  NOTE:
  ------------------------------------------------
  Input image is generated by expanding the original image (image.pgm) via concatenating its parts. The original image needs to be located in the same folder as source
  files.


  ------------------------------------------------
  IMPLEMENTATION-SPECIFIC DESCRIPTION (OpenCL):
  ------------------------------------------------
  This is the OpenCL version of SRAD code.

  In CUDA version of this application, each stage is a separate kernel (due to synchronization requirements) that operates on data already residing in GPU memory.
  In order to improve GPU performance, data was transferred to GPU at the beginning of the code and then transferred back to CPU after all of the computation stages
  were completed in GPU. Some of the kernels use GPU shared memory for additional improvement in performance. Speedup achievable with CUDA version depends on
  the image size (up to the point where GPU saturates). For 1000x1000 image, 100 iterations and 0.5 saturation coefficient, CUDA version achieves approximately
  10x speedup compared to single-threaded C version.

  Speedup numbers reported in the description of this application were obtained on the machine with: Intel Quad Core CPU, 4GB of RAM, Nvidia GTX280 GPU.


  ------------------------------------------------
  RUNNING THIS CODE:
  ------------------------------------------------
  The following are the command parameters to the application:

  1) Number of iterations. Needs to be integer > 0.
  2) Saturation coefficient. Needs to be float > 0.
  3) Number of rows in the input image. Needs to be integer > 0.
  4) Number of columns in the input image. Needs to be integer > 0.

  Example:
  a.out 100 0.5 502 458
	
	Running a.out without parameters will default to the parameters shown in the line above.

*/

// Defines and includes.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include "define.c"
#include "graphics.c"
#include "resize.c"
#include "timer_lin.c"  // use "timer_win.c" if compiling on Windows!
//#include "kernels.h"

// Include the OpenCL header files
#ifdef __APPLE__
#include <OpenCL/cl.h>
#include <OpenCL/cl_gl.h>
#include <OpenCL/cl_gl_ext.h>
#include <OpenCL/cl_ext.h>
#else
#include <CL/cl.h>
#include <CL/cl_gl.h>
#include <CL/cl_gl_ext.h>
#include <CL/cl_ext.h>
#endif

using namespace std;

void checkErr(int ret)
{
	if ((ret != CL_SUCCESS) && (VERBOSE))
		printf("\nError code %i!\n\n", ret);
}

// Program main.
int main(int argc, char* argv[])
{
	// Time...
	long long time0;
	long long time1;
	long long time2;
	long long time3;
	long long time4;
	long long time4b;
	long long time5;
	long long time6;
	long long time6c;
	long long time6d;
	long long time7;
	long long time8;
	long long time9;
	long long time10;
	long long time11;
	long long time12;

	time0 = get_time();

	// Inputs image, input paramenters.
	float* image_ori;
	int image_ori_rows;
	int image_ori_cols;
	long image_ori_elem;

	// ...image, input paramenters.
	float* image;       // input image.
	int Nr, Nc;         // IMAGE nbr of rows/cols/elements.
	long Ne;            // Total number of elements.

	// algorithm parameters
	int niter;          // nbr of iterations
	float lambda;       // update step size

	// size of IMAGE
	int r1, r2, c1, c2; // row/col coordinates of uniform ROI
	long NeROI;         // ROI nbr of elements

	// surrounding pixel indicies
	int *iN, *iS, *jE, *jW;

	// counters
	int iter;           // primary loop
	int i, j;           // image row/col

	// memory sizes
	int mem_size_i;
	int mem_size_j;
	int mem_size_single;

	/* GPU Variables. */
	size_t lwsize;
	int numBlocks;
	int blocks_x;

	/* Kernels */
	cl_kernel extractKernel;
	cl_kernel prepareKernel;
	cl_kernel sradKernel;
	cl_kernel srad2Kernel;
	cl_kernel compressKernel;

	/* memory sizes */
	int mem_size;       // matrix memory size

	/* HOST */
	float total;
	float total2;
	float meanROI;
	float meanROI2;
	float varROI;
	float q0sqr;

	// Misc...
	cl_int ret;

	time1 = get_time();

	/* Get input parameters. Default is 100 0.5 502 458*/
	if (argc != 5)
	{
		printf("Usage: %s ITERATIONS SATURATION Y-SIZE X-SIZE \n", argv[0]);
		printf("--> Using default parameters of 100 0.5 502 458... \n");
		niter = 100;
		lambda = 0.5;
		Nr = 502; // The default y-dimension.
		Nc = 458; // The default x-dimension.
	}
	else
	{
		niter = atoi(argv[1]);
		lambda = atof(argv[2]);
		Nr = atoi(argv[3]); // The y-dimension of the image from the command line.
		Nc = atoi(argv[4]); // The x-dimension of the image from the command line.
	}

	time2 = get_time();

	// Gather the number of pixels from command line for original image.
	image_ori_rows = Nr; // number of pixels in y-direction
	image_ori_cols = Nc; // number of pixels in x-direction
	image_ori_elem = image_ori_rows * image_ori_cols;

	image_ori = (float*) malloc(sizeof(float) * image_ori_elem);

	// Read in the original image file.
	read_graphics("../../data/srad/image.pgm", image_ori, image_ori_rows, image_ori_cols, 1);

	time3 = get_time();

	// Total number of elements in the image...
	Ne = Nr * Nc;

	// The new image that will result from applying SRAD.
	image = (float*) malloc(sizeof(float) * Ne);

	// Resize image. (Assuming column-major storage of image_orig)
	resize(image_ori, image_ori_rows, image_ori_cols, image, Nr, Nc, 1);

	time4 = get_time();

	/* The following code queries the number of platforms and devices, and
	 * lists the information about both.
	 */
	cl_platform_id platform_id[100];
	cl_uint platforms_n = 0;
	clGetPlatformIDs(100, platform_id, &platforms_n);
	if (DISPLAY_INFO)
	{
		printf("\n=== %d OpenCL platform(s) found: ===\n", platforms_n);
		for (int i = 0; i < platforms_n; i++)
		{
			char buffer[10240];
			printf("  -- %d --\n", i);
			clGetPlatformInfo(platform_id[i], CL_PLATFORM_PROFILE, 10240, buffer,
			                  NULL);
			printf("  PROFILE = %s\n", buffer);
			clGetPlatformInfo(platform_id[i], CL_PLATFORM_VERSION, 10240, buffer,
			                  NULL);
			printf("  VERSION = %s\n", buffer);
			clGetPlatformInfo(platform_id[i], CL_PLATFORM_NAME, 10240, buffer, NULL);
			printf("  NAME = %s\n", buffer);
			clGetPlatformInfo(platform_id[i], CL_PLATFORM_VENDOR, 10240, buffer, NULL);
			printf("  VENDOR = %s\n", buffer);
			clGetPlatformInfo(platform_id[i], CL_PLATFORM_EXTENSIONS, 10240, buffer,
			                  NULL);
			printf("  EXTENSIONS = %s\n", buffer);
		}
	}

	cl_device_id device_id[100];
	cl_uint devices_n = 0;
	//clGetDeviceIDs(platform_id[0], CL_DEVICE_TYPE_GPU, 100, device_id, &devices_n);
	clGetDeviceIDs(platform_id[0], CL_DEVICE_TYPE_GPU, 100, device_id, &devices_n);
	if (DISPLAY_INFO)
	{
		printf("=== %d OpenCL device(s) found on platform:\n", devices_n);
		for (int i = 0; i < devices_n; i++)
		{
			char buffer[10240];
			cl_uint buf_uint;
			cl_ulong buf_ulong;
			printf("  -- %d --\n", i);
			clGetDeviceInfo(device_id[i], CL_DEVICE_NAME, sizeof(buffer), buffer,
			                NULL);
			printf("  DEVICE_NAME = %s\n", buffer);
			clGetDeviceInfo(device_id[i], CL_DEVICE_VENDOR, sizeof(buffer), buffer,
			                NULL);
			printf("  DEVICE_VENDOR = %s\n", buffer);
			clGetDeviceInfo(device_id[i], CL_DEVICE_VERSION, sizeof(buffer), buffer,
			                NULL);
			printf("  DEVICE_VERSION = %s\n", buffer);
			clGetDeviceInfo(device_id[i], CL_DRIVER_VERSION, sizeof(buffer), buffer,
			                NULL);
			printf("  DRIVER_VERSION = %s\n", buffer);
			clGetDeviceInfo(device_id[i], CL_DEVICE_MAX_COMPUTE_UNITS,
			                sizeof(buf_uint), &buf_uint, NULL);
			printf("  DEVICE_MAX_COMPUTE_UNITS = %u\n", (unsigned int) buf_uint);
			clGetDeviceInfo(device_id[i], CL_DEVICE_MAX_CLOCK_FREQUENCY,
			                sizeof(buf_uint), &buf_uint, NULL);
			printf("  DEVICE_MAX_CLOCK_FREQUENCY = %u\n", (unsigned int) buf_uint);
			clGetDeviceInfo(device_id[i], CL_DEVICE_GLOBAL_MEM_SIZE,
			                sizeof(buf_ulong), &buf_ulong, NULL);
			printf("  DEVICE_GLOBAL_MEM_SIZE = %llu\n",
			       (unsigned long long) buf_ulong);
			clGetDeviceInfo(device_id[i], CL_DEVICE_LOCAL_MEM_SIZE,
			                sizeof(buf_ulong), &buf_ulong, NULL);
			printf("  CL_DEVICE_LOCAL_MEM_SIZE = %llu\n",
			       (unsigned long long) buf_ulong);
		}
		printf("\n");
	}

	/* Create an OpenCL context. */
	cl_context context = clCreateContext(NULL,       // Properties.
	                                     devices_n,  // The number of entries in device_id[].
	                                     device_id,  // The list of devices.
	                                     NULL,       // pfn_notify, a callback function (not used).
	                                     NULL,       // User data passed  when pfn_notify is called.
	                                     &ret);      // Error code to be returned.
	if ((ret != CL_SUCCESS) && (VERBOSE))
		printf("\nError at clCreateContext! Error code %i\n\n", ret);

	/* Create a command queue. */
	cl_command_queue command_queue = clCreateCommandQueue(context, device_id[0],
	                                                      0, &ret);
	if ((ret != CL_SUCCESS) && (VERBOSE))
		printf("\nError at clCreateCommandQueue! Error code %i\n\n", ret);

	r1 = 0;       // top row index of ROI
	r2 = Nr - 1;  // bottom row index of ROI
	c1 = 0;       // left column index of ROI
	c2 = Nc - 1;  // right column index of ROI

	time4b = get_time();

	/* ROI image size */
	NeROI = (r2 - r1 + 1) * (c2 - c1 + 1); // number of elements in ROI, ROI size

	/* Allocate variables for surrounding pixels */
	mem_size_i = Nr;
	iN = (int *) malloc(mem_size_i * sizeof(int)); // north surrounding element
	iS = (int *) malloc(mem_size_i * sizeof(int)); // south surrounding element
	mem_size_j = Nc;
	jW = (int *) malloc(mem_size_j * sizeof(int)); // west surrounding element
	jE = (int *) malloc(mem_size_j * sizeof(int)); // east surrounding element

	/* N/S/W/E indices of surrounding pixels (every element of IMAGE) */
	for (i = 0; i < Nr; i++)
	{
		iN[i] = i - 1; // holds index of IMAGE row above
		iS[i] = i + 1; // holds index of IMAGE row below
	}
	for (j = 0; j < Nc; j++)
	{
		jW[j] = j - 1; // holds index of IMAGE column on the left
		jE[j] = j + 1; // holds index of IMAGE column on the right
	}

	/* N/S/W/E boundary conditions, fix surrounding indices outside boundary of image */
	iN[0] = 0;           // changes IMAGE top row index from -1 to 0
	iS[Nr - 1] = Nr - 1; // changes IMAGE bottom row index from Nr to Nr-1
	jW[0] = 0;           // changes IMAGE leftmost column index from -1 to 0
	jE[Nc - 1] = Nc - 1; // changes IMAGE rightmost column index from Nc to Nc-1

	time5 = get_time();

	/* Allocate memory for entire IMAGE on DEVICE */
	mem_size = Ne; // get the size of float representation of input IMAGE

	/* The image on the device */
	cl_mem d_I = clCreateBuffer(context, CL_MEM_READ_WRITE, mem_size
	                            * sizeof(cl_float), NULL, &ret);
	if ((ret != CL_SUCCESS) && (VERBOSE))
		printf("\nError at clCreateBuffer! Error code %i \n", ret);

	// Copy from host to device.
	ret = clEnqueueWriteBuffer(command_queue, d_I, CL_TRUE, 0, mem_size
	                           * sizeof(cl_float), image, 0, NULL, NULL);
	if ((ret != CL_SUCCESS) && (VERBOSE))
		printf("\nError at clEnqueueWriteBuffer! Error code %i \n", ret);

	/* Allocate memory for coordinates on DEVICE */
	cl_mem d_iN = clCreateBuffer(context, CL_MEM_READ_WRITE, mem_size_i
	                             * sizeof(cl_int), NULL, NULL);
	clEnqueueWriteBuffer(command_queue, d_iN, CL_TRUE, 0, mem_size_i
	                     * sizeof(cl_int), iN, 0, NULL, NULL);
	cl_mem d_iS = clCreateBuffer(context, CL_MEM_READ_WRITE, mem_size_i
	                             * sizeof(cl_int), NULL, NULL);
	clEnqueueWriteBuffer(command_queue, d_iS, CL_TRUE, 0, mem_size_i
	                     * sizeof(cl_int), iS, 0, NULL, NULL);
	cl_mem d_jE = clCreateBuffer(context, CL_MEM_READ_WRITE, mem_size_j
	                             * sizeof(cl_int), NULL, NULL);
	clEnqueueWriteBuffer(command_queue, d_jE, CL_TRUE, 0, mem_size_j
	                     * sizeof(cl_int), jE, 0, NULL, NULL);
	cl_mem d_jW = clCreateBuffer(context, CL_MEM_READ_WRITE, mem_size_j
	                             * sizeof(cl_int), NULL, NULL);
	clEnqueueWriteBuffer(command_queue, d_jW, CL_TRUE, 0, mem_size_j
	                     * sizeof(cl_int), jW, 0, NULL, NULL);

	/* Allocate memory for partial sums on DEVICE */
	cl_mem d_sums = clCreateBuffer(context, CL_MEM_READ_WRITE, mem_size
	                               * sizeof(cl_float), NULL, NULL);
	cl_mem d_sums2 = clCreateBuffer(context, CL_MEM_READ_WRITE, mem_size
	                                * sizeof(cl_float), NULL, NULL);

	/* Allocate memory for derivatives */
	cl_mem d_dN = clCreateBuffer(context, CL_MEM_READ_WRITE, mem_size
	                             * sizeof(cl_float), NULL, NULL);
	cl_mem d_dS = clCreateBuffer(context, CL_MEM_READ_WRITE, mem_size
	                             * sizeof(cl_float), NULL, NULL);
	cl_mem d_dW = clCreateBuffer(context, CL_MEM_READ_WRITE, mem_size
	                             * sizeof(cl_float), NULL, NULL);
	cl_mem d_dE = clCreateBuffer(context, CL_MEM_READ_WRITE, mem_size
	                             * sizeof(cl_float), NULL, NULL);

	/* Allocate memory for coefficient on DEVICE */
	cl_mem d_c = clCreateBuffer(context, CL_MEM_READ_WRITE, mem_size
	                            * sizeof(cl_float), NULL, NULL);

	time6 = get_time();


	/* Load the source code for all of the kernels into the array source_str */
	FILE* theFile;
	char* source_str;
	size_t source_size;
	theFile = fopen("kernels.cl", "r");
	if (!theFile)
	{
		fprintf(stderr, "Failed to load kernel file.\n");
		exit(1);
	}
	// Obtain length of source file.
	fseek(theFile, 0, SEEK_END);
	source_size = ftell(theFile);
	rewind(theFile);
	// Read in the file.
	source_str = (char*) malloc(sizeof(char) * source_size);
	fread(source_str, 1, source_size, theFile);
	fclose(theFile);

	// Create a program from the kernel source.
	cl_program program = clCreateProgramWithSource(context, 1,
	                                               (const char **) &source_str, NULL, // Number of chars in kernel src. NULL means src is null-terminated.
	                                               &ret);                             // Return status message in the ret variable.

	// Build (compile) the program.
	ret = clBuildProgram(program, NULL, NULL, NULL, NULL, NULL);

	if ((ret != CL_SUCCESS) && (VERBOSE))
		printf("\nError at clBuildProgram! Error code %i\n\n", ret);

	/* Show error info from building the program. */
	if (VERBOSE)
	{
		cout << "\n*************************************************" << endl;
		cout << "***   OUTPUT FROM COMPILING THE KERNEL FILE   ***" << endl;
		cout << "*************************************************" << endl;
		// Shows the log
		char* build_log;
		size_t log_size;
		// First call to know the proper size
		clGetProgramBuildInfo(program, device_id[0], CL_PROGRAM_BUILD_LOG, 0, NULL,
		                      &log_size);
		build_log = new char[log_size + 1];
		// Second call to get the log
		clGetProgramBuildInfo(program, device_id[0], CL_PROGRAM_BUILD_LOG,
		                      log_size, build_log, NULL);
		build_log[log_size] = '\0';
		cout << build_log << endl;
		delete[] build_log;
		cout << "\n*************************************************" << endl;
		cout << "*** END OUTPUT FROM COMPILING THE KERNEL FILE ***" << endl;
		cout << "*************************************************\n\n" << endl;
	}

	time6c = get_time();

	// Set the global work size.
	size_t globalWorkSize[] = { Ne };

	// The "extract" kernel.
	extractKernel = clCreateKernel(program, "extract", NULL);
	clSetKernelArg(extractKernel, 0, sizeof(cl_long), (void *) &Ne);
	clSetKernelArg(extractKernel, 1, sizeof(cl_mem), (void *) &d_I);
	
	// Get the kernel work group size.
	clGetKernelWorkGroupInfo(extractKernel, device_id[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &lwsize, NULL);
	int howManyThreads = lwsize;
	while (Ne % howManyThreads != 0)
	{
		howManyThreads--;
	}
	if (VERBOSE)
	{
		cout << "Max local threads is " << lwsize << ". ";
		cout << "Using " << howManyThreads << "for local work size." << endl;
	}
	size_t localWorkSize[] = { howManyThreads };
	if (VERBOSE)
	{
		cout << "===========================" << endl;
		cout << "Using local work size " << howManyThreads << endl;
		cout << "===========================" << endl;
	}

	clFinish(command_queue);
	time6d = get_time();

	// Scale from 0-255 down to 0-1, and extract.
	ret = clEnqueueNDRangeKernel(command_queue, extractKernel, 1, NULL,
	                             globalWorkSize, localWorkSize, 0, NULL, NULL);
	if ((ret != CL_SUCCESS) && (VERBOSE))
		printf("\nError code %i when attempting to launch extractKernel\n\n", ret);

	clFinish(command_queue);
	time7 = get_time();

	/* MAIN COMPUTATION LOOP */
	for (iter = 0; iter < niter; iter++)
	{
		// The "prepare" kernel.
		prepareKernel = clCreateKernel(program, "prepare", NULL);
		checkErr(clSetKernelArg(prepareKernel, 0, sizeof(cl_long), (void *) &Ne));
		checkErr(clSetKernelArg(prepareKernel, 1, sizeof(cl_mem), (void *) &d_I));
		checkErr(clSetKernelArg(prepareKernel, 2, sizeof(cl_mem), (void *) &d_sums));
		checkErr(
			clSetKernelArg(prepareKernel, 3, sizeof(cl_mem), (void *) &d_sums2));

		// Launch the "prepare" kernel.
		ret = clEnqueueNDRangeKernel(command_queue, prepareKernel, 1, NULL,
		                             globalWorkSize, localWorkSize, 0, NULL, NULL);
		if ((ret != CL_SUCCESS) && (VERBOSE))
			printf("\nError code %i when attempting to launch prepareKernel\n\n", ret);

		// Figure out the number of thread blocks that will be used.
		numBlocks = Ne / (howManyThreads);
		if (Ne % (howManyThreads) != 0)
		{
			printf(
				"Error: number of local threads must be a number which divides evenly into the global number of threads.\n");
			printf("Local work size: %i\nGlobal work size: %i\n", 
			       (int)localWorkSize[0],
			       (int)globalWorkSize[0]);
		}

		// q0sqr is basically a threshold on edge detection. If set too low, it will create the spots in your examples. Too high, it will over-smooth.
		q0sqr = 114957;                   

		// Make sure all prior kernels have finished.
		clFinish(command_queue);

		// Set up the "srad" kernel and set parameters.
		sradKernel = clCreateKernel(program, "srad", NULL);
		clSetKernelArg(sradKernel, 0, sizeof(cl_float), (void*) &lambda);        // SRAD coefficient
		clSetKernelArg(sradKernel, 1, sizeof(cl_int), (void*) &Nr);              // # of rows in input image
		clSetKernelArg(sradKernel, 2, sizeof(cl_int), (void*) &Nc);              // # of columns in input image
		clSetKernelArg(sradKernel, 3, sizeof(cl_long), (void*) &Ne);             // # of elements in input image
		clSetKernelArg(sradKernel, 4, sizeof(cl_mem), (void*) &d_iN);            // indices of North surrounding pixels
		clSetKernelArg(sradKernel, 5, sizeof(cl_mem), (void*) &d_iS);            // indices of South surrounding pixels
		clSetKernelArg(sradKernel, 6, sizeof(cl_mem), (void*) &d_jE);            // indices of East surrounding pixels
		clSetKernelArg(sradKernel, 7, sizeof(cl_mem), (void*) &d_jW);            // indices of West surrounding pixels
		clSetKernelArg(sradKernel, 8, sizeof(cl_mem), (void*) &d_dN);            // North derivative
		clSetKernelArg(sradKernel, 9, sizeof(cl_mem), (void*) &d_dS);            // South derivative
		clSetKernelArg(sradKernel, 10, sizeof(cl_mem), (void*) &d_dW);           // West derivative
		clSetKernelArg(sradKernel, 11, sizeof(cl_mem), (void*) &d_dE);           // East derivative
		clSetKernelArg(sradKernel, 12, sizeof(cl_float), (void*) &q0sqr);        // standard deviation of ROI
		clSetKernelArg(sradKernel, 13, sizeof(cl_mem), (void*) &d_c);            // diffusion coefficient
		clSetKernelArg(sradKernel, 14, sizeof(cl_mem), (void*) &d_I);            // output image
		clSetKernelArg(sradKernel, 15, sizeof(cl_int), (void*) &howManyThreads); // Thread block size.

		// Set up the "srad2" kernel and set parameters.
		srad2Kernel = clCreateKernel(program, "srad2", NULL);
		clSetKernelArg(srad2Kernel, 0, sizeof(cl_float), (void*) &lambda);        // SRAD coefficient
		clSetKernelArg(srad2Kernel, 1, sizeof(cl_int), (void*) &Nr);              // # of rows in input image
		clSetKernelArg(srad2Kernel, 2, sizeof(cl_int), (void*) &Nc);              // # of columns in input image
		clSetKernelArg(srad2Kernel, 3, sizeof(cl_long), (void*) &Ne);             // # of elements in input image
		clSetKernelArg(srad2Kernel, 4, sizeof(cl_mem), (void*) &d_iN);            // indices of North surrounding pixels
		clSetKernelArg(srad2Kernel, 5, sizeof(cl_mem), (void*) &d_iS);            // indices of South surrounding pixels
		clSetKernelArg(srad2Kernel, 6, sizeof(cl_mem), (void*) &d_jE);            // indices of East surrounding pixels
		clSetKernelArg(srad2Kernel, 7, sizeof(cl_mem), (void*) &d_jW);            // indices of West surrounding pixels
		clSetKernelArg(srad2Kernel, 8, sizeof(cl_mem), (void*) &d_dN);            // North derivative
		clSetKernelArg(srad2Kernel, 9, sizeof(cl_mem), (void*) &d_dS);            // South derivative
		clSetKernelArg(srad2Kernel, 10, sizeof(cl_mem), (void*) &d_dW);           // West derivative
		clSetKernelArg(srad2Kernel, 11, sizeof(cl_mem), (void*) &d_dE);           // East derivative
		clSetKernelArg(srad2Kernel, 12, sizeof(cl_mem), (void*) &d_c);            // diffusion coefficient
		clSetKernelArg(srad2Kernel, 13, sizeof(cl_mem), (void*) &d_I);            // output image
		clSetKernelArg(srad2Kernel, 14, sizeof(cl_int), (void*) &howManyThreads); // Thread block size.

		// execute srad kernel
		ret = clEnqueueNDRangeKernel(command_queue, sradKernel, 1, NULL,
		                             globalWorkSize, localWorkSize, 0, NULL, NULL);
		if ((ret != CL_SUCCESS) && (VERBOSE))
			printf("\nError code %i when attempting to launch sradKernel\n\n", ret);

		// Wait for the srad to finish...
		clFinish(command_queue);

		// execute srad2 kernel
		ret = clEnqueueNDRangeKernel(command_queue, srad2Kernel, 1, NULL,
		                             globalWorkSize, localWorkSize, 0, NULL, NULL);
		if ((ret != CL_SUCCESS) && (VERBOSE))
			printf("\nError code %i when attempting to launch srad2Kernel\n\n", ret);

	}

	// Make sure all prior kernels have finished.
	clFinish(command_queue);

	// Get the time.
	time8 = get_time();

	// The "compress" kernel.
	compressKernel = clCreateKernel(program, "compress", NULL);
	clSetKernelArg(compressKernel, 0, sizeof(cl_long), (void *) &Ne);
	clSetKernelArg(compressKernel, 1, sizeof(cl_mem), (void *) &d_I);
	clSetKernelArg(compressKernel, 2, sizeof(cl_int), (void *) &howManyThreads);

	// Scale image up from 0-1 to 0-255 and compress.
	ret = clEnqueueNDRangeKernel(command_queue, compressKernel, 1, NULL,
	                             globalWorkSize, localWorkSize, 0, NULL, NULL);
	if ((ret != CL_SUCCESS) && (VERBOSE))
		printf("\nError code %i when attempting to launch compressKernel\n\n", ret);

	clFinish(command_queue);

	// Get the time.
	time9 = get_time();

	// Copy results back to host.
	ret = clEnqueueReadBuffer(command_queue,               // The command queue.
	                          d_I,                         // The image on the device.
	                          CL_TRUE,                     // Blocking? (ie. Wait at this line until read has finished?)
	                          0,                           // Offset. None in this case.
	                          mem_size * sizeof(cl_float), // Size to copy.
	                          image,                       // The pointer to the image on the host.
	                          0,                           // Number of events in wait list. Not used.
	                          NULL,                        // Event wait list. Not used.
	                          NULL);                       // Event object for determining status. Not used.
	if ((ret != CL_SUCCESS) && (VERBOSE))
		printf(
			"\nError code %i when attempting to copy d_I from device to host!\n\n",
			ret);

	// Make sure all prior kernels have finished.
	clFinish(command_queue);

	// Get the time.
	time10 = get_time();

	// Write image after processing.
	write_graphics("image_out.pgm", image, Nr, Nc, 1, 255);

	time11 = get_time();

	// Clean up memory.
	free(image_ori);
	free(image);
	free(iN);
	free(iS);
	free(jW);
	free(jE);

	/* Clean up OpenCL stuff. */
	clFlush(command_queue);
	clFinish(command_queue);
	// Release kernels...
	clReleaseKernel(extractKernel);
	clReleaseKernel(prepareKernel);
	clReleaseKernel(sradKernel);
	clReleaseKernel(srad2Kernel);
	clReleaseKernel(compressKernel);
	// Now the program...
	clReleaseProgram(program);
	// Clean up the device memory...
	clReleaseMemObject(d_I);
	clReleaseMemObject(d_c);
	clReleaseMemObject(d_iN);
	clReleaseMemObject(d_iS);
	clReleaseMemObject(d_jE);
	clReleaseMemObject(d_jW);
	clReleaseMemObject(d_dN);
	clReleaseMemObject(d_dS);
	clReleaseMemObject(d_dE);
	clReleaseMemObject(d_dW);
	clReleaseMemObject(d_sums);
	clReleaseMemObject(d_sums2);
	// ...and finally, the queue and context.
	clReleaseCommandQueue(command_queue);
	clReleaseContext(context);

	time12 = get_time();

	// display the timing results.
	printf("Time spent in different stages of the application:\n");

	printf("%.12f s, %.12f %% : SETUP VARIABLES\n", 
	       (float)(time1 - time0) / 1000000,
	       (float)(time1 - time0) / (float)(time12 - time0) * 100);

	printf("%.12f s, %.12f %% : READ COMMAND LINE PARAMETERS\n",
	       (float)(time2 - time1) / 1000000,
	       (float)(time2 - time1) / (float)(time12 - time0) * 100);

	printf("%.12f s, %.12f %% : READ IMAGE FROM FILE\n", (float) (time3 - time2)
	       / 1000000, (float) (time3 - time2) / (float) (time12 - time0) * 100);
	printf("%.12f s, %.12f %% : RESIZE IMAGE\n",
	       (float) (time4 - time3) / 1000000, (float) (time4 - time3)
	       / (float) (time12 - time0) * 100);
	printf("%.12f s, %.12f %% : GET DEVICES, CREATE CONTEXT AND QUEUE\n",
	       (float) (time4b - time4) / 1000000, (float) (time4b - time4)
	       / (float) (time12 - time0) * 100);
	printf("%.12f s, %.12f %% : CPU/GPU SETUP, CPU/GPU MEMORY ALLOCATION\n",
	       (float) (time5 - time4b) / 1000000, (float) (time5 - time4)
	       / (float) (time12 - time0) * 100);
	printf("%.12f s, %.12f %% : COPY DATA TO CPU->GPU\n", (float) (time6 - time5)
	       / 1000000, (float) (time6 - time5) / (float) (time12 - time0) * 100);
	printf("%.12f s, %.12f %% : LOAD KERNEL FROM EXTERNAL FILE\n", (float) (time6c
		  - time6) / 1000000, (float) (time6c - time6) / (float) (time12 - time0)
	       * 100);
	printf("%.12f s, %.12f %% : CREATE KERNEL, SET ARGUMENTS AND WORK SIZE.\n",
	       (float) (time6d - time6c) / 1000000, (float) (time6d - time6c)
	       / (float) (time12 - time0) * 100);
	printf("%.12f s, %.12f %% : EXTRACT IMAGE\n", (float) (time7 - time6)
	       / 1000000, (float) (time7 - time6c) / (float) (time12 - time0) * 100);
	printf("%.12f s, %.12f %% : COMPUTE\n", (float) (time8 - time7) / 1000000,
	       (float) (time8 - time7) / (float) (time12 - time0) * 100);
	printf("%.12f s, %.12f %% : COMPRESS IMAGE\n", (float) (time9 - time8)
	       / 1000000, (float) (time9 - time8) / (float) (time12 - time0) * 100);
	printf("%.12f s, %.12f %% : COPY DATA TO GPU->CPU\n", (float) (time10 - time9)
	       / 1000000, (float) (time10 - time9) / (float) (time12 - time0) * 100);
	printf("%.12f s, %.12f %% : SAVE IMAGE INTO FILE\n", (float) (time11 - time10)
	       / 1000000, (float) (time11 - time10) / (float) (time12 - time0) * 100);
	printf("%.12f s, %.12f %% : FREE MEMORY\n", (float) (time12 - time11)
	       / 1000000, (float) (time12 - time11) / (float) (time12 - time0) * 100);
	printf("Total time:\n");
	printf("%.12f s\n", (float) (time12 - time0) / 1000000);
}
