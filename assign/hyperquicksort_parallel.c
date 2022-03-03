#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "hyperquicksort_parallel.h"




/***************** MAIN MECHANISM (INTERFACE FOR PSRS (Parallel Sorting Regular Sampling)) **********************/



/* sort an array in non-descending order */
void parl_hyperQuickSort(uint_64 *a, int processors, int n) {

	/* base cases - that a regular sort proved to be better */
	if (processors <= 1 || n <= 1000) {
		quickSort(a, 0, n - 1);
		return;
	}
	
	processors = MIN(processors, omp_get_max_threads()); // making sure we don't overload our threads with too many partitionings, that will make the CPUs inefficient

  int p, size, rsize, sample_size;
  uint_64 *sample, *pivots;
  int *partition_borders, *bucket_sizes, *result_positions;
  uint_64 **loc_a_ptrs;

  // Determine the appropriate number of threads to use
  // p^3 <= n - We need this to hold true
  p = processors;
  p = p*p*p;
  if(p > n){
	p = floor(pow(n,0.33));
	p-=p%2;
  }else{
   p = processors;
   p-=p%2;
  }
  omp_set_num_threads(p);

	/* control panel */
  size  = (n + p - 1) / p; // calculate size of local arrays (must add  (p-1)  in order to round up things
  rsize = (size + p - 1) / p; // in a local array, the distance between each sample
  sample_size = p * (p - 1); // the size of the sampling array ( (p-1) pivots from each thread. There are p threads)

  loc_a_ptrs = (uint_64**) malloc(p * sizeof(uint_64 *)); // array of pointers to local arrays
  sample = (uint_64*) malloc(sample_size * sizeof(uint_64)); // array of the gathered samples
  partition_borders = (int*) malloc(p * (p + 1) * sizeof(int)); // for (p-1) pivots, there will be (p+1) borders for each local array, making p * (p+1) borders in total
  bucket_sizes = (int*) malloc(p * sizeof(int));
  result_positions = (int*) malloc(p * sizeof(int)); // offsets, within the origin array a, that tell where each "cummulative partition" starts at. That helps when extracting the i-th partition of each thread's local array into one big array
  pivots = (uint_64*) malloc((p - 1) * sizeof(uint_64)); // the (p - 1) sorted pivots chosen from "sample"


	/* Parallel Region - spawned with p threads */
  #pragma omp parallel
  {
    int i, j, max, thread_num, start, end, loc_size, offset, cummulative_partition_size;
    uint_64 *loc_a, *cummulative_partition, *current_a;

	get_specs(&start, &end, &loc_size, &thread_num, size, n); // fetch local size, start, and end
	Extract_Part(&loc_a, a, start, loc_size, loc_a_ptrs, thread_num); // extract a part of the array, and sort it
	
	/* populating the sample array, with p * (p - 1) samples. Populate by the threads_num */
	end = end % size; // relative end location (rather than origin end location) of the local array
	Extract_Samples(sample, loc_a, thread_num, end, rsize, p);
    #pragma omp barrier

	/* sort the sample array, and choose (p - 1) pivots from it, and populate it into the "pivots" array */
    #pragma omp single
    Sort_Samples(sample, sample_size, pivots, p);
    // implied Barrier (OpenMP)
	
	
	
	/* calculate the borders of each partition of the thread's local array, and inject them into the "partition_borders" shared array */
    // enter the obvious borders = the start and the end
    offset = thread_num * (p + 1);
    partition_borders[offset] = 0;
    partition_borders[offset + p] = end + 1;
    // calculate the "inner" borders of the thread's local array
    Extract_Partition_Borders(loc_a, 0, loc_size-1, partition_borders, offset, pivots, 1, p-1);
    #pragma omp barrier
    

	/* calculate how long each "cummulative partition" is going to be */
	Partition_Length(&bucket_sizes[thread_num], thread_num, partition_borders, p * (p + 1), p);
    #pragma omp barrier


	/* calculate the relative position of this thread's "cummulative partition", in regards to the start of the origin array 'a' */
    #pragma omp single
	Calc_Locations(result_positions, bucket_sizes, p);
	// implied Barrier (OpenMP)
	
	
	/* calculate the length of this "cummulative partition", and then accumulate the  */
	cummulative_partition_size = (thread_num == (p-1)) ? (n - result_positions[thread_num]) : (bucket_sizes[thread_num]);    
	
	/* accumulate the <thread_num>-th partition of each part, and then gather them into their place in the original array */
	Accumulate_Partitions(&cummulative_partition, a, thread_num, partition_borders, result_positions[thread_num], loc_a_ptrs, p);

	/* Sort the accumulated partition, and free the memory that this thread allocated */
    quickSort(cummulative_partition, 0, cummulative_partition_size - 1);
	#pragma omp barrier
    free(loc_a); // free local array
  }

	/* free program */
  free(loc_a_ptrs);
  free(sample);
  free(partition_borders);
  free(bucket_sizes);
  free(result_positions);
  free(pivots);

}







/************************** MAIN MECHANISM'S AUXILIARY ********************************/






/* calculate the borders of a local array, with a given array of pivots, the start and the end of the portion of the array
 * that you want to determine its borders with, the "ranks" of the pivots you wish to determine lower borders for
 * Done in a recursive manner. Recursion depth is O(log(end - start)), and each iteration costs O(log(end - start)) */
void Extract_Partition_Borders(uint_64 array[],    // array being sorted
                            int start,			
                            int end,              // separate the array into current process range
                            int partition_borders[],			// the partition_borders array
                            int offset_in_arr,               // this thread's local array's start at the partition_borders array
                            uint_64 pivots[],   // the pivot values
                            int first_pivot,         // first 
                            int last_pivot)          // last pivot
{
  int mid, bottom_border;
  uint_64 curr_pivot;

  mid = (first_pivot + last_pivot) / 2; // get the middle pivot rating between [1, ... , (p-1)]
  curr_pivot = pivots[mid-1]; // get the value of the pivot
  
  /* binary search for the first value in the array that is greater than the pivot */
	bottom_border = BinarySearch(array, curr_pivot, start, end);
  
  // enter the index of the first value in the local array that is greater than the pivot, as the border of the pivot, from its left (lowerbound) */
  partition_borders[offset_in_arr + mid] = bottom_border; 

	/* recursion */
  if(first_pivot < mid) { // search for borders implied by pivots who are *lower* than the current pivots[mid]. Since they're lower, we can reduce the end to the (lowerbound - 1)
    Extract_Partition_Borders(array, start, bottom_border - 1, partition_borders, offset_in_arr, pivots, first_pivot, mid - 1);
  }
  if(mid < last_pivot) { // search for borders implied by pivots who are *bigger* than the current pivots[mid]. Since they're bigger, we can reduce the start to (lowerbound)
    Extract_Partition_Borders(array, bottom_border, end, partition_borders, offset_in_arr, pivots, mid + 1, last_pivot);
  }
}






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
	/* create a local array, and memcpy into it from the origin array, from an offset of "start", and size "loc_size" */
    *local_part = (uint_64*) malloc(size_copy * sizeof(uint_64));
    memcpy(*local_part, original_arr + start_copy, size_copy * sizeof(uint_64));
    local_parts_array[thread_num] = *local_part; // saving a pointer to the thread's local array in the shared "pointers-array"

	/* sort the local array */
	quickSort(*local_part, 0, size_copy - 1);


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








/* This function gets as an argument, the size of each "cummulative partition", and like that,
 * Calculates the absolute location of such i-th, where the absolute location is simply where the "cummulative partition" would start at the resulted array */ 
void Calc_Locations(int* partitions_locations, int* cummulative_partition_sizes, int parts) {
	  partitions_locations[0] = 0;
	  for(int i = 1; i < parts; i++) {
		partitions_locations[i] = cummulative_partition_sizes[i-1] + partitions_locations[i-1];
	  }
}
/****************************************************************************************/







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
    uint_64
     temp = *xp;
    *xp = *yp;
    *yp = temp;
}
/******************************************************/
