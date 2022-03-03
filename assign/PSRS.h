#ifndef PSRS_H
#define PSRS_H



/* making life easier since I compile with C11, and it doesn't know uint64_t, 
 * yet "unsigned long" is quite a pain to maintain */ 
typedef unsigned long uint_64;


/* general purpose auxiliary functions */
void assert_other(int, char[]);



/* MIN and MAX macros. Used for their superficial destiny */
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))



/**********************************************************************************************
 *********** MAIN MECHANISM (INTERFACE FOR PSRS (Parallel Sorting Regular Sampling)) **********
 **********************************************************************************************/
 
 
/* Given an array and an amount of processors, We will distribute parts of the array into the different threads
 * that will be spawnedd, and let them sort each local part of the array that they control.
 * Once done, we will gather an amount of <processors - 1> samples from each local sorted array part, 
 * Into one big samples array. After done, We will gather <processors - 1> of pivots from the samples array (kind like Hoare's median of three rule).
 * Now, we will make each thread partition its local part using the pivots (well, each local part is already sorted, so we just need to find the borders of each partition)
 * We have now come to the result of having <processors> partitions for each local part thrad. Now gets the fun part:
 * For the thread of rank <thread_num>, we will copy each and every thread's local part's <thread_num>-th partition, into a new location.
 * We have now accumulated <processors> amount of partitions, and we sort them, using regular quicksort (any algorithm is welcomed).
 * Now, we have got that the thread of rank <thread_num>, hold a sorted accumulated array of every number in the original array,
 * that is smaller than the <thread_num>-th pivot. So, all there is left to do, is to link every sorted accumulated array of each thread
 * To the thread behind and next to it, but that part of the algorithm is done throughout the running (Link means to just make their heads and tails touch). */
void parl_hyperQuickSort(uint_64*, int, int);


/* The Parallel Region of the algorithm */
void Parallel_Region (uint_64* original_array,
					  int original_array_size,
				      int part_size, 
					  int local_sample_distance, 
					  int samples_amount, 
					  uint_64** local_parts_array, 
					  uint_64* samples, 
					  int* partition_borders, 
					  int* cummulative_partition_sizes, 
					  int* cummulative_partition_locations, 
					  uint_64* pivots,
					  int processors);
				
			
			
			
/************************************************************************** 
 *********************** MAIN MECHANISM'S AUXILIARY ***********************
 **************************************************************************/
 
/* Calculate the optimal amount of threads */
void Calc_Processors_Amount(int* processors, int array_size);

/* get crucial information about the running thread: rank, indicies of its borders at the original array, size, etc */
void get_specs(int* start,
			   int* end,
			   int* size,
			   int* thread_num, 
			   int part_size,
			   int origin_length);

/* Extract a part of the original array. Extract <size_copy> elements from the original array into the new one, and store the new local thread's array into the
 * an array that maintains pointers to all of the local arrays.
 * We then sort the array that we parted */ 
void Extract_Part(uint_64** local_part,
				  uint_64* original_arr, 
				  int start_copy,
				  int size_copy, 
				  uint_64** local_parts_array, 
				  int thread_num);	
				  
				  
				  
					  
					  

/********************************************************************************************************************************************** 
 *********************** GATHER SAMPLES FROM THREADS' LOCAL PARTS, SORT IT, EXTRACT GLOBAL PIVOTS *********************************************
 **********************************************************************************************************************************************/		
 
 /* This function extracts at a maximum, <parts> samples from the array_send, into the array_recv. Each sample must reside at least <sample_dist> indices away */
 void Extract_Samples(uint_64* array_recv,
					 uint_64* array_send,
					 int thread_num,
					 int end,
					 int sample_dist,
					 int parts);
					 	  
/* given the array of sample we've gathered - Sort it, And extract <parts - 1> pivots from it */
void Sort_Samples(uint_64* samples, int sample_size, uint_64* pivots, int parts);
					
					
					
					
					
					
/******************************************************************************************************************************************************
 ********************************** GIVEN <PROCESSORS-1> PIVOTS: Calculae The Borders induced from those partitions ***********************************
 ******************************************************************************************************************************************************/

/* this function returns the first value in the array that is greater than the pivot */
int BinarySearch(uint_64* array, uint_64 pivot, int left, int right);
	
	
/* calculate the borders of a local part array, for the pivots that are within the range of [first_pivot, last_pivot].
 * We are guaranteed that those pivots' borders are within the range [start, ..., end] of the local_part_array.
 * first_pivot and last_pivot are just the index of the pivot when the pivots are put into a sorted array.
 * Done in a recursive manner. Recursion depth is O(log(end - start)), and each iteration costs O(log(end - start)).
 * For optimization manners, since this part is done in a recursive Divide-and-Conquer manner,
 * Once a bottom_border is found for a partition of a pivot, the recursion that seraches for pivots higher than the current pivot,
 * is going to search for borders in the range [bottom_border, ..., end] of the array.
 * the recursion that seraches for pivots lower than the current pivot, is gonig to search for borders in the range [start, ..., bottom_border - 1] of the array  */				 	  
void Extract_Partition_Borders(uint_64* local_part_array,
                            int start,			
                            int end, 
                            int* partition_borders,
                            int offset_in_arr, 
                            uint_64* pivots,
                            int first_pivot,
                            int last_pivot);
               





/******************************************************************************************************************************************************************** 
 ********************************** LET THE THREAD GATHER THE PARTITIONS HE MANAGES, FROM OTHER THREADs' LOCAL PARTS ************************************************
 ********************************************************************************************************************************************************************/



/* Gathers all of the partitions of all of the threads, that correlate to the partition that this thread was destined to manage. 
 * It then, copies that array into the place of the original array, into the fitting place (relative to the lower partitions' size that they occupy) */
void Accumulate_Partitions(uint_64** cummulative_partition,
						   uint_64* original_array, 
						   int thread_num,
						   int* partition_borders,
						   int start_copy,
						   uint_64** local_parts_array,
						   int parts);
						   
						   
/* This function gets a thread_num id, and then sums up all of the sizes of the <thread_num>-th partitions, in each thread's local array
 * That would be equal to a thread's "cummulative partition"'s length */			   
void Partition_Length(int* length, 
					  int thread_num,
					  int* partition_borders,
					  int all_partition_borders,
					  int parts);



/************************************************************************************************************************************************************************* 
 ********************************** GIVEN ACCUMULATED PARTITIONS, FIND THEIR ABSOLUTE LOCATION INSIDE THE NEW ARRAY ******************************************************
 *************************************************************************************************************************************************************************/

/* This function gets as an argument, the size of each "cummulative partition", and like that,
 * Calculates the absolute location of such i-th, where the absolute location is simply where the "cummulative partition" would start at the resulted array */
void Calc_Locations(int* partitions_locations, int* cummulative_partition_sizes, int parts);





/****************************************************************** 
 *********************** INTERFACE for QUICKSORT ******************
 ******************************************************************/

/* quick-sort regular algorithm */
void quickSort(uint_64*, int, int);

/* hoare's partition */
int partition(uint_64*, int, int);

/* swap for quick sort */
void swap(uint_64 *, uint_64 *);






#endif
