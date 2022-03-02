#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "hyperquicksort_parallel.h"



/* headers */
void insertion_sort(uint_64 *arr, int n);
void calc_partition_borders(uint_64 array[],
               int start,
               int end,
               int sublist_sizes[],
               int at,
               uint_64 pivots[],
               int first_p,
               int last_p);
void sortll(uint_64 *a, int len);
uint_64 *merge(uint_64 * left, uint_64 * right, int l_end, int r_end);
uint_64 *merge_sort(uint_64 * arr, int size);

/* auxiliary funcs declaration */
void quickSort(uint_64*, int, int);
int partition(uint_64*, int, int);
void swap(uint_64 *, uint_64 *);





/***************** MAIN MECHANISM **********************/



/* sort an array in non-descending order */
void parl_hyperQuickSort(uint_64 *a, int p, int n) {

  int size, rsize, sample_size;
  uint_64 *sample, *pivots;
  int *partition_borders, *bucket_sizes, *result_positions;
  uint_64 **loc_a_ptrs;

  // Determine the appropriate number of threads to use
  // p^3 <= n - We need this to hold true
  p = p*p*p;
  if(p > n){
	p = floor(pow(n,0.33));
	p-=p%2;
  }else{
   p = omp_get_max_threads();
   p-=p%2;
  }
  omp_set_num_threads(p);

  size  = (n + p - 1) / p;
  rsize = (size + p - 1) / p;
  sample_size = p * (p - 1);

  loc_a_ptrs = malloc(p * sizeof(uint_64 *));
  sample = malloc(sample_size * sizeof(uint_64));
  partition_borders = malloc(p * (p + 1) * sizeof(int));
  bucket_sizes = malloc(p * sizeof(int));
  result_positions = malloc(p * sizeof(int));
  pivots = malloc((p - 1) * sizeof(uint_64));

  #pragma omp parallel
  {
    int i, j, max, thread_num, start, end, loc_size, offset, this_result_size;
    uint_64 *loc_a, *this_result, *current_a;

    thread_num = omp_get_thread_num();
    start = thread_num * size;
    end = start + size - 1;
    if(end >= n) end = n - 1;
    loc_size = (end - start + 1);
    end = end % size;

    loc_a = malloc(loc_size * sizeof(uint_64));
    memcpy(loc_a, a + start, loc_size * sizeof(uint_64));
    loc_a_ptrs[thread_num] = loc_a;

    sortll(loc_a, loc_size); // Testing shows that this sequential sort is quickest in this instance
	
    offset = thread_num * (p - 1) - 1;

    for(i = 1; i < p; i++) {
      if(i * rsize <= end) {
        sample[offset + i] = loc_a[i * rsize - 1];
      } else {
        sample[offset + i] = loc_a[end];
      }
    }

    #pragma omp barrier

    #pragma omp single
    {
      merge_sort(sample, sample_size); // Testing shows that this sequential sort is quickest in this instance
      for(i = 0; i < p - 1; i++) {
        pivots[i] = sample[i * p + p / 2];
      }
    }

    #pragma omp barrier

    offset = thread_num * (p + 1);
    partition_borders[offset] = 0;
    partition_borders[offset + p] = end + 1;
    calc_partition_borders(loc_a, 0, loc_size-1, partition_borders, offset, pivots, 1, p-1);

    #pragma omp barrier

    max = p * (p + 1);
    bucket_sizes[thread_num] = 0;
    for(i = thread_num; i < max; i += p + 1) {
      bucket_sizes[thread_num] += partition_borders[i + 1] - partition_borders[i];
    }

    #pragma omp barrier

    #pragma omp single
    {
      result_positions[0] = 0;
      for(i = 1; i < p; i++) {
        result_positions[i] = bucket_sizes[i-1] + result_positions[i-1];
      }
    }

    #pragma omp barrier

    this_result = a + result_positions[thread_num];

    if(thread_num == p-1) {
      this_result_size = n - result_positions[thread_num];
    } else {
      this_result_size = result_positions[thread_num+1] - result_positions[thread_num];
    }

    // pluck this threads sublist from each of the local arrays
    this_result = a + result_positions[thread_num];

    for(i = 0, j = 0; i < p; i++) {
      int low, high, partition_size;
      offset = i * (p + 1) + thread_num;
      low = partition_borders[offset];
      high = partition_borders[offset+1];
      partition_size = (high - low);
      if(partition_size > 0) {
        memcpy(this_result+j, &(loc_a_ptrs[i][low]), partition_size * sizeof(uint_64));
        j += partition_size;
      }
    }

    // sort p local sorted arrays
    sortll(this_result, this_result_size); // Testing shows that this sequential sort is quickest in this instance
	
    #pragma omp barrier
    free(loc_a);
  }

  free(loc_a_ptrs);
  free(sample);
  free(partition_borders);
  free(bucket_sizes);
  free(result_positions);
  free(pivots);

}


/* determine the boundaries for the sublists of an local array */
void calc_partition_borders(uint_64 array[],    // array being sorted
                            int start,
                            int end,              // separate the array into current process range
                            int result[],
                            int at,               // this process start point in result
                            uint_64 pivots[],   // the pivot values
                            int first_pv,         // first pivot
                            int last_pv)          // last pivot
{
  int mid, lowerbound, upperbound, center;
  uint_64 pv;

  mid = (first_pv + last_pv) / 2;
  pv = pivots[mid-1];
  lowerbound = start;
  upperbound = end;
  while(lowerbound <= upperbound) {
    center = (lowerbound + upperbound) / 2;
    if(array[center] > pv) {
      upperbound = center - 1;
    } else {
      lowerbound = center + 1;
    }
  }
  result[at + mid] = lowerbound;

  if(first_pv < mid) {
    calc_partition_borders(array, start, lowerbound - 1, result, at, pivots, first_pv, mid - 1);
  }
  if(mid < last_pv) {
    calc_partition_borders(array, lowerbound, end, result, at, pivots, mid + 1, last_pv);
  }
}




/*
  Sort a portion of an array
*/
void sortll(uint_64 *a, int len)
{
	quickSort(a, 0, len);
  // qsort(a, len, sizeof(uint_64), lcompare);
}


/******************************************************************/


/*********************** INTERFACE for QUICKSORT ******************/


/* sorting the given file using multi-core parallelism */
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



/*
  Standard merge sort
*/
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

uint_64 *merge(uint_64 * left, uint_64 * right, int l_end, int r_end){
	int temp_off, l_off, r_off, size = l_end+r_end;
	uint_64 *temp = malloc(sizeof(uint_64) * l_end);

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
