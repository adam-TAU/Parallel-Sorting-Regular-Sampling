#ifndef QUICKSORT_PARALLEL_H
#define QUICKSORT_PARALLEL_H



/* making life easier since I compile with C11, and it doesn't know uint64_t, 
 * yet "unsigned long" is quite a pain to maintain */ 
typedef unsigned long uint_64;




/* Main Mechanism's function declarations */
void parl_hyperQuickSort(uint_64*, int, int);
void calc_partition_borders(uint_64 array[],
               int start,
               int end,
               int sublist_sizes[],
               int at,
               uint_64 pivots[],
               int first_p,
               int last_p);
void Extract_Samples(uint_64* array_recv, uint_64* array_send, int thread_num, int end, int sample_dist, int parts);




/* declarations of quicksort's (auxiliary) funcs */
void quickSort(uint_64*, int, int);
int partition(uint_64*, int, int);
void swap(uint_64 *, uint_64 *);




/* general purpose auxiliary functions */
int BinarySearch(uint_64* array, uint_64 pivot, int left, int right);
void assert_other(int, char[]);





/* MIN and MAX macros. Used for their superficial destiny */
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


#endif
