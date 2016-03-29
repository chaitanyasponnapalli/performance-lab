//Chaitanya Ponnapalli
//006948176
#include <stdio.h>
#include "cs1300bmp.h"
#include <iostream>
#include <stdint.h>
#include <fstream>
#include "Filter.h"
#include <omp.h>

using namespace std;

#include "rtdsc.h"

//
// Forward declare the functions
//
Filter * readFilter(string filename);
double applyFilter(Filter *filter, cs1300bmp *input, cs1300bmp *output);

int
main(int argc, char **argv)
{

  if ( argc < 2) {
    fprintf(stderr,"Usage: %s filter inputfile1 inputfile2 .... \n", argv[0]);
  }

  //
  // Convert to C++ strings to simplify manipulation
  //
  string filtername = argv[1];

  //
  // remove any ".filter" in the filtername
  //
  string filterOutputName = filtername;
  string::size_type loc = filterOutputName.find(".filter");
  if (loc != string::npos) {
    //
    // Remove the ".filter" name, which should occur on all the provided filters
    //
    filterOutputName = filtername.substr(0, loc);
  }

  Filter *filter = readFilter(filtername);

  double sum = 0.0;
  int samples = 0;

  for (int inNum = 2; inNum < argc; inNum++) {
    string inputFilename = argv[inNum];
    string outputFilename = "filtered-" + filterOutputName + "-" + inputFilename;
    struct cs1300bmp *input = new struct cs1300bmp;
    struct cs1300bmp *output = new struct cs1300bmp;
    int ok = cs1300bmp_readfile( (char *) inputFilename.c_str(), input);

    if ( ok ) {
      double sample = applyFilter(filter, input, output);
      sum += sample;
      samples++;
      cs1300bmp_writefile((char *) outputFilename.c_str(), output);
    }
    delete input;
    delete output;
  }
  fprintf(stdout, "Average cycles per sample is %f\n", sum / samples);

}

struct Filter *
readFilter(string filename)
{
  ifstream input(filename.c_str());

  if ( ! input.bad() ) {
    int size = 0;
    input >> size;
    Filter *filter = new Filter(size);
    int div;
    input >> div;
    filter -> setDivisor(div);
    for (int i=0; i < size; i++) {
      for (int j=0; j < size; j++) {
	int value;
	input >> value;
	filter -> set(i,j,value);
      }
    }
    return filter;
  }
}

#if defined(__arm__)
static inline unsigned int get_cyclecount (void)
{
 unsigned int value;
 // Read CCNT Register
 asm volatile ("MRC p15, 0, %0, c9, c13, 0\t\n": "=r"(value)); 
 return value;
}

static inline void init_perfcounters (int32_t do_reset, int32_t enable_divider)
{
 // in general enable all counters (including cycle counter)
 int32_t value = 1;

 // peform reset: 
 if (do_reset)
 {
   value |= 2;     // reset all counters to zero.
   value |= 4;     // reset cycle counter to zero.
 }

 if (enable_divider)
   value |= 8;     // enable "by 64" divider for CCNT.

 value |= 16;

 // program the performance-counter control-register:
 asm volatile ("MCR p15, 0, %0, c9, c12, 0\t\n" :: "r"(value)); 

 // enable all counters: 
 asm volatile ("MCR p15, 0, %0, c9, c12, 1\t\n" :: "r"(0x8000000f)); 

 // clear overflows:
 asm volatile ("MCR p15, 0, %0, c9, c12, 3\t\n" :: "r"(0x8000000f));
}



#endif

double
applyFilter(struct Filter *filter, cs1300bmp *input, cs1300bmp *output)
{
	#if defined(__arm__)
	init_perfcounters (1, 1);
	#endif

	long long cycStart, cycStop;
	double start,stop;
	#if defined(__arm__)
	cycStart = get_cyclecount();
	#else
	cycStart = rdtscll();
	#endif
	//variables kept out of the function and removed from the loops, so they are initialized only once to remove redundancy of calculations
	int width = input -> width;
	int height = input -> height;
	int in_width = (input -> width) - 1;
	int in_height = (input -> height) - 1;
	output -> width = width;
	output -> height = height;	
	//function call is performed outside the loop, so that it is called once instead of for every iteration
	int filterDivisor = filter -> getDivisor();
	//store each filtered pixel values after unrolling the two loops used in filtering the image
	int filterCol[9];
	//declared the row, col and plane as short instead of int
	short row, col, plane;
	int colorPlane;
	//Move the call to 'get' outside and unroll the loop, so that boundary conditions need not be checked, during assigning of the values
	int filterPlaneVal[3][3];
	filterPlaneVal[0][0] = filter -> get(0, 0);
	filterPlaneVal[0][1] = filter -> get(0, 1);
	filterPlaneVal[0][2] = filter -> get(0, 2);
	filterPlaneVal[1][0] = filter -> get(1, 0);
	filterPlaneVal[1][1] = filter -> get(1, 1);
	filterPlaneVal[1][2] = filter -> get(1, 2);
	filterPlaneVal[2][0] = filter -> get(2, 0);
	filterPlaneVal[2][1] = filter -> get(2, 1);
	filterPlaneVal[2][2] = filter -> get(2, 2);
	
	//#pragma omp parallel for //- This distorted the images, so removed it
	//To increase the spatial locality, the iterations are performed over each plane(inside which the row and col are traversed)
	for(plane = 0; plane < 3; plane++) {
		//#pragma omp parallel for - However significantly slowed down the filtering process
		//loops are interchanged so that the iteration is in row major order, supported by c. It is more efficient to access the elements in the order in which they are stored
		for(row = 1; row < in_height ; row++) {			
				//performed loop unrolling, to reduce the redundancy of boundary tests and time taken during iterations. These instructions can also be executed in parallel 
				col = 1;//tried to reduce the number of loop iterations and thereby calls to struct member input -> color[plane][row][col], to speed up the process. This sped up the process by a factor of 40 clock cycles. Also, this can reduce the calls as the values can be re-assigned, for the loop					
				
				filterCol[0] = input -> color[plane][row-1][col-1];
				filterCol[1] = input -> color[plane][row-1][col];
				filterCol[2] = input -> color[plane][row-1][col+1];
				filterCol[3]= input -> color[plane][row][col-1];
				filterCol[4]= input -> color[plane][row][col];
				filterCol[5]= input -> color[plane][row][col+1];
				filterCol[6]= input -> color[plane][row+1][col-1];
				filterCol[7]= input -> color[plane][row+1][col];
				filterCol[8]= input -> color[plane][row+1][col+1];
			
				//stored the values in local variables to reduce the access time. The array declared above is used as an accumulator, so that calculations can be performed in parallel, without much dependency. 
				colorPlane = filterCol[0] * filterPlaneVal[0][0]+ filterCol[1] * filterPlaneVal[0][1]+ filterCol[2] * filterPlaneVal[0][2]+filterCol[3] * filterPlaneVal[1][0]+filterCol[4] * filterPlaneVal[1][1] +filterCol[5] * filterPlaneVal[1][2] +filterCol[6] * filterPlaneVal[2][0]+filterCol[7] * filterPlaneVal[2][1] +filterCol[8] * filterPlaneVal[2][2];

				//reduce the unnecessary clock cycles for division, if divisor is equals 1
				if ( filterDivisor > 1){
					colorPlane = colorPlane/filterDivisor;
				}				
				//used ternary operator to remove unnecessary if-else branching
				output -> color[plane][row][col] = colorPlane < 0   ? 0   : (colorPlane > 255 ? 255 : colorPlane);									
			//permits re-usability of the values and reduces calls to input -> color[plane][row][col]			
			for(col = 2; col < in_width; col++) {
				filterCol[0] = filterCol[1];//assigns the calculated value of col value 1
				filterCol[1] = filterCol[2];
				filterCol[2] = input -> color[plane][row-1][col+1];
				filterCol[3]= filterCol[4];
				filterCol[4]= filterCol[5];
				filterCol[5]= input -> color[plane][row][col+1];
				filterCol[6]= filterCol[7];
				filterCol[7]= filterCol[8];
				filterCol[8]=input -> color[plane][row+1][col+1];
				colorPlane = filterCol[0] * filterPlaneVal[0][0]+ filterCol[1] * filterPlaneVal[0][1]+ filterCol[2] * filterPlaneVal[0][2]+filterCol[3] * filterPlaneVal[1][0]+filterCol[4] * filterPlaneVal[1][1] +filterCol[5] * filterPlaneVal[1][2] +filterCol[6] * filterPlaneVal[2][0]+filterCol[7] * filterPlaneVal[2][1] +filterCol[8] * filterPlaneVal[2][2];
				if ( filterDivisor > 1){
					colorPlane = colorPlane / filterDivisor;
				}				
				output -> color[plane][row][col] = colorPlane;				
				//placed this as this would be the more probable branch that would be taken. This prevents unnecessary boundary checking conditions
				if ( colorPlane > 0 && colorPlane < 255)
					continue;			
				output -> color[plane][row][col] = colorPlane < 0   ? 0   : (colorPlane > 255 ? 255 : colorPlane);				
			}
		}
	}			
	#if defined(__arm__)
	cycStop = get_cyclecount();
	#else
	cycStop = rdtscll();
	#endif

	double diff = cycStop-cycStart;
	diff = diff * 64;
	double diffPerPixel = diff / (output -> width * output -> height);
	fprintf(stderr, "Took %f cycles to process, or %f cycles per pixel\n",
	  diff, diff / (output -> width * output -> height));
	return diffPerPixel;
}
