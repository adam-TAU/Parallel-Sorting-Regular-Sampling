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

  #pragma omp parallel
  {
    int i, j, max, thread_num, start, end, loc_size, offset, this_result_size;
    uint_64 *loc_a, *this_result, *current_a;

	/* fetch local size, start, and end */
    thread_num = omp_get_thread_num(); // rank of thread
    start = thread_num * size; // traditional start
    end = start + size - 1; // traditional end 
    end = (end >= n) ? (n - 1) : end; // if (start + size) > n: this is the last thread. Move end to the end of the origin array
    loc_size = (end - start + 1); // calculate the local size of the array we're working on, through end and start (since "size" isn't accurate at last thread
    end = end % size; // relative end location (rather than origin end location) of the local array



	/* create a local array, and memcpy into it from the origin array, from an offset of "start", and size "loc_size" */
    loc_a = (uint_64*) malloc(loc_size * sizeof(uint_64));
    memcpy(loc_a, a + start, loc_size * sizeof(uint_64));
    loc_a_ptrs[thread_num] = loc_a; // saving a pointer to the thread's local array in the shared "pointers-array"

	/* sort the local array */
	quickSort(loc_a, 0, loc_size - 1);
	
	
	/* populating the sample array, with p * (p - 1) samples. Populate by the threads_num */
	Extract_Samples(sample, loc_a, thread_num, end, rsize, p);
    #pragma omp barrier


	/* sort the sample array, and choose (p - 1) pivots from it, and populate it into the "pivots" array */
    #pragma omp single
    {
      merge_sort(sample, sample_size); /* the sorting algorithm of the sample array */
      for(i = 0; i < p - 1; i++) {
        pivots[i] = sample[i * p + p / 2]; // the amount of processors, on average, is going to exclude p/2 numbers into the last thread's local array. Therefore, we'll have on average p/2 zeros at the start of the array. 0 is bad pivot.
      }
    }
    // there's already an implied barrier here (thanks OpenMP)
	
	
	
	/* calculate the borders of each partition of the thread's local array, and inject them into the "partition_borders" shared array */
    // enter the obvious borders = the start and the end
    offset = thread_num * (p + 1);
    partition_borders[offset] = 0;
    partition_borders[offset + p] = end + 1;
    // calculate the "inner" borders of the thread's local array
    calc_partition_borders(loc_a, 0, loc_size-1, partition_borders, offset, pivots, 1, p-1);
    #pragma omp barrier



	/* calculate how long each "cummulative partition" is going to be */
    max = p * (p + 1);
    bucket_sizes[thread_num] = 0; // initialize current "cummulative partition" to be 0
    for(i = thread_num; i < max; i += p + 1) { // every <thread_num>-s partition of some local array of a thread, is going to be owned by this "cummulative partition"
      bucket_sizes[thread_num] += partition_borders[i + 1] - partition_borders[i]; // calculate the local partition size and sum it up
    }
    #pragma omp barrier




    #pragma omp single
    {
      result_positions[0] = 0;
      for(i = 1; i < p; i++) {
        result_positions[i] = bucket_sizes[i-1] + result_positions[i-1];
      }
    }
	// there's already an implied barrier here (thanks OpenMP)	
	
	
	

	/* calculate the length of this "cummulative partition" */
	this_result_size = (thread_num == (p-1)) ? (n - result_positions[thread_num]) : (result_positions[thread_num+1] - result_positions[thread_num]);    
    


    /* Extract the partitions in each thread, that are supposed to be handled by this thread (so all of the <thread_num> partitions, of all threads, shall be merged into one, here */
    this_result = a + result_positions[thread_num];

    for(i = 0, j = 0; i < p; i++) {
      int low, high, partition_size;
      
      // get the i-th thread's, <thread_num> partition borders. While the 0-th partition is starts at index 0.
      offset = i * (p + 1) + thread_num; 
      low = partition_borders[offset];
      high = partition_borders[offset+1];
      partition_size = (high - low);
      
      // copy the partition (using its borders - relative to the i-th thread's local array) into the final position of the cummulative <thread_num> partition (0th partition starts at index 0)
      if(partition_size > 0) {
        memcpy(this_result+j, &(loc_a_ptrs[i][low]), partition_size * sizeof(uint_64));
        j += partition_size;
      }
    }


    // sort the local array inplace, at its new placement (the origin array)
    quickSort(this_result, 0, this_result_size - 1);
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


/* calculate the borders of a local array, with a given array of pivots, the start and the end of the portion of the array
 * that you want to determine its borders, the "ranks" of the pivots you wish to determine lower borders for
 * Done in a recursive manner. Recursion depth is O(log(end - start)), and each iteration costs O(log(end - start)) */
void calc_partition_borders(uint_64 array[],    // array being sorted
                            int start,			
                            int end,              // separate the array into current process range
                            int result[],			// the partition_borders array
                            int at,               // this thread's local array's start at the partition_borders array
                            uint_64 pivots[],   // the pivot values
                            int first_pv,         // first 
                            int last_pv)          // last pivot
{
  int mid, lowerbound;
  uint_64 pv;

  mid = (first_pv + last_pv) / 2; // get the middle pivot rating between [1, ... , (p-1)]
  pv = pivots[mid-1]; // get the value of the pivot
  
  /* binary search for the first value in the array that is greater than the pivot */
	lowerbound = BinarySearch(array, pv, start, end);
  
  // enter the index of the first value in the local array that is greater than the pivot, as the border of the pivot, from its left (lowerbound) */
  result[at + mid] = lowerbound; 

	/* recursion */
  if(first_pv < mid) { // search for borders implied by pivots who are *lower* than the current pivots[mid]. Since they're lower, we can reduce the end to the (lowerbound - 1)
    calc_partition_borders(array, start, lowerbound - 1, result, at, pivots, first_pv, mid - 1);
  }
  if(mid < last_pv) { // search for borders implied by pivots who are *bigger* than the current pivots[mid]. Since they're bigger, we can reduce the start to (lowerbound)
    calc_partition_borders(array, lowerbound, end, result, at, pivots, mid + 1, last_pv);
  }
}






/* this function returns the first value in the array that is greater than the pivot */
int BinarySearch(uint_64* array, uint_64 pivot, int left, int right) {

	int center;
  
  while(left <= right) {
    center = (left + right) / 2;
    if(array[center] > pivot) {
      right = center - 1;
    } else {
      left = center + 1;
    }
  }
  
  return left;
}



/* This function extracts at a maximum, <parts> samples from the array_send, into the array_recv. Each sample must reside at least <sample_dist> indices away */
void Extract_Samples(uint_64* array_recv, uint_64* array_send, int thread_num, int end, int sample_dist, int parts) {

    int offset = thread_num * (parts - 1) - 1; // where the thread will write its samples (into the collective sampling array)
    
    for(int i = 1; i < parts; i++) { // parts, 2 * parts, ..., (parts-1) * parts
      if(i * sample_dist <= end) {
        array_recv[offset + i] = array_send[i * sample_dist - 1]; // i * p sample
      } else {
        array_recv[offset + i] = array_send[end]; // if the local array was the last thread's, hence a shortended array, there aren't enough samples ( < (p-1) )
      }
    }
}

/******************************************************************/







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










/*********************** INTERFACE for MERGESORT ******************/

/* merge-sort regular algorithm */
uint_64 *merge_sort(uint_64 * arr, int size){
	// Arrays shorter than 1 are already sorted
	if(size > 1){
		int middle = size / 2, i; 
		uint_64 *left, *right;
		left = arr;
		right = arr + middle; 
		
		left = merge_sort(left, middle);
		right = merge_sort(right, size-middle);
		return merge(left, right, middle,size-middle);
	}else { return arr; }
}

/* merging two sorted lists */
uint_64 *merge(uint_64 * left, uint_64 * right, int l_end, int r_end){
	int temp_off, l_off, r_off, size = l_end+r_end;
	uint_64 *temp = (uint_64*) malloc(sizeof(uint_64) * l_end);

	// Copy lower half into temp buffer
	for(l_off=0, temp_off=0; left+l_off != right; l_off++, temp_off++){
		*(temp + temp_off) = *(left + l_off);
	}
	
	temp_off=0; l_off=0; r_off=0;

	while(l_off < size){
		if(temp_off < l_end){
			if(r_off < r_end){
				if(*(temp+temp_off) < *(right+r_off)){
					*(left+l_off) = *(temp+temp_off);
					temp_off++;
				}else{
					*(left+l_off) = *(right+r_off);
					r_off++;
				}
			}else{
				*(left+l_off) = *(temp+temp_off);
				temp_off++;
			}
		}else{
			if(r_off < r_end) {
				*(left + l_off) = *(right + r_off);
				r_off++;
			}else{
				printf("\nERROR - merging loop going too far\n");
			}
		}
		l_off++;
	}
	free(temp);
	return left;
}



/******************************************************/
