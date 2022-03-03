#ifndef QUICKSORT_PARALLEL_H
#define QUICKSORT_PARALLEL_H



/* making life easier since I compile with C11, and it doesn't know uint64_t, 
 * yet "unsigned long" is quite a pain to maintain */ 
typedef unsigned long uint_64;




/* Main Mechanism's function declarations */
void parl_hyperQuickSort(uint_64*, int, int);
void Extract_Partition_Borders(uint_64 array[],
                            int start,			
                            int end, 
                            int partition_borders[],
                            int offset_in_arr, 
                            uint_64 pivots[],
                            int first_pivot,
                            int last_pivot);
               
void Extract_Samples(uint_64* array_recv,
					 uint_64* array_send,
					 int thread_num,
					 int end,
					 int sample_dist,
					 int parts);


void get_specs(int* start,
			   int* end,
			   int* size,
			   int* thread_num, 
			   int part_size,
			   int origin_length);


void Extract_Part(uint_64** local_part,
				  uint_64* original_arr, 
				  int start_copy,
				  int size_copy, 
				  uint_64** local_parts_array, 
				  int thread_num);



void Sort_Samples(uint_64* samples, int sample_size, uint_64* pivots, int parts);



void Accumulate_Partitions(uint_64** cummulative_partition,
						   uint_64* original_array, 
						   int thread_num,
						   int* partition_borders,
						   int start_copy,
						   uint_64** local_parts_array,
						   int parts);
						   
						   
void Partition_Length(int* length, 
					  int thread_num,
					  int* partition_borders,
					  int all_partition_borders,
					  int parts);




void Calc_Locations(int* partitions_locations, int* cummulative_partition_sizes, int parts);


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
