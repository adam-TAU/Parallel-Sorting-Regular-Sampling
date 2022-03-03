#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "hyperquicksort_parallel.h"




/***************** MAIN MECHANISM (INTERFACE FOR PSRS (Parallel Sorting Regular Sampling)) **********************/



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
void parl_hyperQuickSort(uint_64 *a, int processors, int n) {

	/* base cases - that a regular sort proved to be better */
	if (processors <= 1 || n <= 1000) {
		quickSort(a, 0, n - 1);
		return;
	}
	
	Calc_Processors_Amount(&processors, n);

	/* declarations of global shared variables */
  int part_size, local_sample_distance, samples_amount;
  uint_64 *samples, *pivots;
  int *partition_borders, *cummulative_partition_sizes, *cummulative_partition_locations;
  uint_64 **local_parts_array;



	/* definitions of global shared variables */
  part_size  = (n + processors - 1) / processors; // calculate size of local arrays (must add  (p-1)  in order to round up things
  local_sample_distance = (part_size + processors - 1) / processors; // in a local array, <local_sample_distance> is the distance between each sample
  samples_amount = processors * (processors - 1); // the size of the sampling array ( (p-1) pivots from each thread. There are p threads)

  local_parts_array = (uint_64**) malloc(processors * sizeof(uint_64 *)); // array of pointers to local arrays
  samples = (uint_64*) malloc(samples_amount * sizeof(uint_64)); // array of the gathered samples
  partition_borders = (int*) malloc(processors * (processors + 1) * sizeof(int)); // for (p-1) pivots, there will be (p+1) borders for each local array, making p * (p+1) borders in total. This global shared array holds those border indices.
  cummulative_partition_sizes = (int*) malloc(processors * sizeof(int));
  cummulative_partition_locations = (int*) malloc(processors * sizeof(int)); // offsets, within the origin array a, that tell where each "cummulative partition" starts at. That helps when extracting the i-th partition of each thread's local array into one big array
  pivots = (uint_64*) malloc((processors - 1) * sizeof(uint_64)); // the (p - 1) sorted pivots chosen from "sample"



	/* The parallel Region */
	Parallel_Region(a,
				   n,
				   part_size,
				   local_sample_distance, 
				   samples_amount, 
				   local_parts_array, 
				   samples, 
				   partition_borders, 
				   cummulative_partition_sizes, 
				   cummulative_partition_locations, 
				   pivots,
				   processors);


	/* free program */
  free(local_parts_array);
  free(samples);
  free(partition_borders);
  free(cummulative_partition_sizes);
  free(cummulative_partition_locations);
  free(pivots);

}







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
					  int processors)

{	
	/* Parallel Region - spawned with <processors> threads */
  #pragma omp parallel
  {
    int thread_num, start, end, local_part_size, offset, cummulative_partition_size;
    uint_64 *local_part, *cummulative_partition;

	get_specs(&start, &end, &local_part_size, &thread_num, part_size, original_array_size); // fetch local size, start, and end
	Extract_Part(&local_part, original_array, start, local_part_size, local_parts_array, thread_num); // extract a part of the array, and sort it
	
	/* populating the samples array, with p * (p - 1) samples. Populate by the threads_num */
	end = end % part_size; // relative end location (rather than origin end location) of the local array
	Extract_Samples(samples, local_part, thread_num, end, local_sample_distance, processors);
    #pragma omp barrier

	/* sort the sample array, and choose (p - 1) pivots from it, and populate it into the "pivots" array */
    #pragma omp single
    Sort_Samples(samples, samples_amount, pivots, processors); // sort the gathered samples and derive <processors - 1> pivots from it
    // implied Barrier (OpenMP)
	
	
	
	/* calculate the borders of each partition of the thread's local array, and inject them into the "partition_borders" shared array */
    // enter the obvious borders = the start and the end
    offset = thread_num * (processors + 1);
    partition_borders[offset] = 0;
    partition_borders[offset + processors] = end + 1;
    // calculate the borders of the thread's local part array
    Extract_Partition_Borders(local_part, 0, local_part_size-1, partition_borders, offset, pivots, 1, processors-1);
    #pragma omp barrier
    

	/* calculate how long each "cummulative partition" is going to be */
	Partition_Length(&cummulative_partition_sizes[thread_num], thread_num, partition_borders, processors * (processors + 1), processors);
    #pragma omp barrier


	/* calculate the relative position of this thread's "cummulative partition", in regards to the start of the origin array 'a' */
    #pragma omp single
	Calc_Locations(cummulative_partition_locations, cummulative_partition_sizes, processors);
	// implied Barrier (OpenMP)
	
	
	/* calculate the length of this "cummulative partition", and then accumulate the  */
	cummulative_partition_size = (thread_num == (processors-1)) ? (original_array_size - cummulative_partition_locations[thread_num]) : (cummulative_partition_sizes[thread_num]);    
	
	/* accumulate the <thread_num>-th partition of each part, and then gather them into their place in the original array */
	Accumulate_Partitions(&cummulative_partition, original_array, thread_num, partition_borders, cummulative_partition_locations[thread_num], local_parts_array, processors);

	/* Sort the accumulated partition, and free the memory that this thread allocated */
    quickSort(cummulative_partition, 0, cummulative_partition_size - 1);
	#pragma omp barrier
    free(local_part); // free local array
  }
}






/**************************************************** MAIN MECHANISM'S AUXILIARY **********************************************************/

/* Calculate the optimal amount of threads */
void Calc_Processors_Amount(int* processors, int array_size) {
	int max_processors = MIN((*processors), omp_get_max_threads()); // making sure we don't overload our threads with too many partitionings, that will make the CPUs inefficient


  // processors^3 <= n - We need this to hold true
  *processors = pow(*processors, 3);
  if(*processors > array_size){
	*processors = floor(pow(array_size,0.33));
	*processors -= (*processors % 2);
  }else{
   *processors = max_processors;
   *processors -= (*processors % 2);
  }
  omp_set_num_threads(*processors);
}



/* get crucial information about the running thread: rank, indicies of its borders at the original array, size, etc */
void get_specs(int* start,
			   int* end,
			   int* size,
			   int* thread_num, 
			   int part_size,
			   int origin_length) 	   
{
    *thread_num = omp_get_thread_num(); // rank of thread
    *start = (*thread_num) * part_size; // traditional start
    *end = *start + part_size - 1; // traditional end 
    *end = ((*end) >= origin_length) ? (origin_length - 1) : *end; // if (start + size) > n: this is the last thread. Move end to the end of the origin array
    *size = (*end - *start + 1); // calculate the local size of the array we're working on, through end and start (since "size" isn't accurate at last thread
}






/* Extract a part of the original array. Extract <size_copy> elements from the original array into the new one, and store the new local thread's array into the
 * an array that maintains pointers to all of the local arrays.
 * We then sort the array that we parted */ 
void Extract_Part(uint_64** local_part,
				  uint_64* original_arr, 
				  int start_copy,
				  int size_copy, 
				  uint_64** local_parts_array, 
				  int thread_num)
{
	/* create a local array, and memcpy into it from the origin array, from an offset of "start", and size "local_part_size" */
    *local_part = (uint_64*) malloc(size_copy * sizeof(uint_64));
    memcpy(*local_part, original_arr + start_copy, size_copy * sizeof(uint_64));
    local_parts_array[thread_num] = *local_part; // saving a pointer to the thread's local array in the shared "pointers-array"

	/* sort the local array */
	quickSort(*local_part, 0, size_copy - 1);
}
/*************************************************************************************************************************************************/












/******************************************* GATHER SAMPLES FROM THREADS' LOCAL PARTS, SORT IT, EXTRACT GLOBAL PIVOTS ***********************************************/

/* This function extracts at a maximum, <parts> samples from the array_send, into the array_recv. Each sample must reside at least <sample_dist> indices away */
void Extract_Samples(uint_64* array_recv,
					 uint_64* array_send,
					 int thread_num,
					 int end,
					 int sample_dist,
					 int parts)
					 
{

    int offset = thread_num * (parts - 1) - 1; // where the thread will write its samples (into the collective sampling array)
    
    for(int i = 1; i < parts; i++) { // parts, 2 * parts, ..., (parts-1) * parts
      
      if (i * sample_dist - 1 <= end) array_recv[offset + i] = array_send[i * sample_dist - 1]; // i * p sample
      else array_recv[offset + i] = array_send[end]; // if the local array was the last thread's, hence a shortended array, there aren't enough samples ( < (p-1) )
    }
}




/* given the array of sample we've gathered - Sort it, And extract <parts - 1> pivots from it */
void Sort_Samples(uint_64* samples, int sample_size, uint_64* pivots, int parts) {

      quickSort(samples, 0, sample_size - 1); /* the sorting algorithm of the sample array */
      for(int i = 0; i < parts - 1; i++) {
        // the amount of processors, on average, is going to exclude p/2 numbers into the last thread's local array.
        // Therefore, we'll have on average p/2 zeros at the start of the array. 0 is bad pivot.
        pivots[i] = samples[i * parts + parts / 2]; 
      }
}
/*********************************************************************************************************************************************************************/













/**************************************** GIVEN <PROCESSORS-1> PIVOTS: Calculae The Borders induced from those partitions ******************************************/

/* this function returns the first value in the array that is greater than the pivot */
int BinarySearch(uint_64* array, uint_64 pivot, int left, int right) {

	int center;
  
  while(left <= right) {
    center = (left + right) / 2;
    
    if (array[center] > pivot) right = center - 1;
    else left = center + 1;
  }
  
  return left;
}



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
                            int* partition_borders,			// the partition_borders array
                            int offset_in_arr,               // this thread's local array's start at the partition_borders array
                            uint_64* pivots,
                            int first_pivot,         
                            int last_pivot)          
{
  int mid, bottom_border;
  uint_64 curr_pivot;

  mid = (first_pivot + last_pivot) / 2; // get the middle pivot rating between [1, ... , (p-1)]
  curr_pivot = pivots[mid-1]; // get the value of the current pivot
  
  /* binary search for the first value in the array that is greater than the pivot */
	bottom_border = BinarySearch(local_part_array, curr_pivot, start, end);
  
  // enter the index of the first value in the local array that is greater than the pivot, as the border of the pivot, from its left (lowerbound) */
  partition_borders[offset_in_arr + mid] = bottom_border; 

	/* recursion */
  if(first_pivot < mid) { // search for borders implied by pivots who are *lower* than the current pivots[mid]. Since they're lower, we can reduce the end to the (lowerbound - 1)
    Extract_Partition_Borders(local_part_array, start, bottom_border - 1, partition_borders, offset_in_arr, pivots, first_pivot, mid - 1);
  }
  if(mid < last_pivot) { // search for borders implied by pivots who are *bigger* than the current pivots[mid]. Since they're bigger, we can reduce the start to (lowerbound)
    Extract_Partition_Borders(local_part_array, bottom_border, end, partition_borders, offset_in_arr, pivots, mid + 1, last_pivot);
  }
}
/***************************************************************************************************************************************/














/******************************************* LET THE THREAD GATHER THE PARTITIONS HE MANAGES, FROM OTHER THREADs' LOCAL PARTS ******************************************************/

/* Gathers all of the partitions of all of the threads, that correlate to the partition that this thread was destined to manage. 
 * It then, copies that array into the place of the original array, into the fitting place (relative to the lower partitions' size that they occupy) */
void Accumulate_Partitions(uint_64** cummulative_partition,
						   uint_64* original_array, 
						   int thread_num,
						   int* partition_borders,
						   int start_copy,
						   uint_64** local_parts_array,
						   int parts) 
{
    /* Extract the partitions in each thread, that are supposed to be handled by this thread (so all of the <thread_num> partitions, of all threads, shall be merged into one, here */
	*cummulative_partition = original_array + start_copy;
	
    for(int i = 0, j = 0; i < parts; i++) {
      int bottom, top, thread_partition_size, offset;
      
      // get the i-th thread's, <thread_num> partition borders. While the 0-th partition is starts at index 0.
      offset = i * (parts + 1) + thread_num; 
      bottom = partition_borders[offset];
      top = partition_borders[offset+1];
      thread_partition_size = (top - bottom);
      
      // copy the partition (using its borders - relative to the i-th thread's local array) into the final position of the cummulative <thread_num> partition (0th partition starts at index 0)
      if(thread_partition_size > 0) {
        memcpy(*cummulative_partition+j, ((local_parts_array)[i] + bottom), thread_partition_size * sizeof(uint_64));
        j += thread_partition_size;
      }
    }
}





/* This function gets a thread_num id, and then sums up all of the sizes of the <thread_num>-th partitions, in each thread's local array
 * That would be equal to a thread's "cummulative partition"'s length */
void Partition_Length(int* length, 
					  int thread_num,
					  int* partition_borders,
					  int all_partition_borders,
					  int parts)
{
	*length = 0; // initialize current "cummulative partition" to be 0
    for(int i = thread_num; i < all_partition_borders; i += (parts + 1)) { // every <thread_num>-s partition of some local array of a thread, is going to be owned by this "cummulative partition"
      *length += partition_borders[i + 1] - partition_borders[i]; // calculate the local partition size and sum it up
    }
}
/********************************************************************************************************************************************************************/












/******************************************************* GIVEN ACCUMULATED PARTITIONS, FIND THEIR ABSOLUTE LOCATION INSIDE THE NEW ARRAY ******************************************************/

/* This function gets as an argument, the size of each "cummulative partition", and like that,
 * Calculates the absolute location of such i-th, where the absolute location is simply where the "cummulative partition" would start at the resulted array */ 
void Calc_Locations(int* partitions_locations, int* cummulative_partition_sizes, int parts) {
	  partitions_locations[0] = 0;
	  for(int i = 1; i < parts; i++) {
		partitions_locations[i] = cummulative_partition_sizes[i-1] + partitions_locations[i-1];
	  }
}
/****************************************************************************************************************************************************************************************************/











/*********************** INTERFACE for QUICKSORT ******************/


/* quick-sort regular algorithm */
void quickSort(uint_64* arr, int low, int high)
{

    if (low < high) {
        /* pi is partitioning index, arr[p] is now
         * at right place */
        uint_64 pi = partition(arr, low, high);
 
        /* Separately sort elements before
         * partition and after partition */
        quickSort(arr, low, pi);
        quickSort(arr, pi + 1, high);
    }
}



/* hoare's partition */
int partition(uint_64* arr, int low, int high)
{
    uint_64 pivot = arr[low];
    int i = low - 1;
    int j = high + 1;
    while (1)
    {
        do {
            i++;
        } while (arr[i] < pivot);
 
        do {
            j--;
        } while (arr[j] > pivot);
 
        if (i >= j) {
            return j;
        }
 
        swap(&arr[i], &arr[j]);
    }
}




/* swap for quick sort */
void swap(uint_64 *xp, uint_64 *yp)
{
    uint_64 temp = *xp;
    *xp = *yp;
    *yp = temp;
}
/******************************************************/
