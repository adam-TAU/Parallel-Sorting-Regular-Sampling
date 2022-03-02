#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>


/* global variable declaration */
int processors;
int amount;
unsigned long *nums;
/******************************/


/* auxiliary funcs declaration */
void quickSort(int, int);
int partition(int, int);
void swap(unsigned long *, unsigned long *);
void parse_args(FILE*);
void printarr(unsigned long[], int);
void free_program();
void assert_other(int, char[]);
/******************************/





int main(int argv, char* args[]) {

	/* parse arguments */
	assert_other(argv != 3, "Incorrect amount of arguments!"); 
	int proc = atoi(args[1]);
	FILE* ifp = fopen(args[2], "r");
	parse_args(ifp);

	/* sort + measure time */
	struct timeval start, end;
	gettimeofday(&start, NULL);
	quickSort(0, amount - 1);
	gettimeofday(&end, NULL);
	double elapsed = (end.tv_sec - start.tv_sec)*1000000.0 + (end.tv_usec - start.tv_usec);

	/* print to stdout */
	printf("QuickSort: %.lf\n", elapsed);
    printarr(nums, amount);	


	/* coherently finish the program */
	free_program();
}






/****************  MAIN MECHANISM  *********************/






/* sorting the given file using multi-core parallelism */
void quickSort(int low, int high)
{

    if (low < high) {
        /* pi is partitioning index, arr[p] is now
         * at right place */
        int pi = partition(low, high);
 
        /* Separately sort elements before
         * partition and after partition */
        quickSort(low, pi);
        quickSort(pi + 1, high);
    }
}



/* hoare's partition */
int partition(int low, int high)
{
    unsigned long pivot = nums[low];
    int i = low - 1;
    int j = high + 1;
    while (1)
    {
        do {
            i++;
        } while (nums[i] < pivot);
 
        do {
            j--;
        } while (nums[j] > pivot);
 
        if (i >= j) {
            return j;
        }
 
        swap(&nums[i], &nums[j]);
    }
}




/* swap for quick sort */
void swap(unsigned long *xp, unsigned long *yp)
{
    unsigned long temp = *xp;
    *xp = *yp;
    *yp = temp;
}
/******************************************************/








/***************  GENERIC FUNCTIONS  *****************/



/* parsing crucial data of the input file */
void parse_args(FILE* ifp) {
	unsigned long i;
	int ind;	

	/* initialize the array of nums first */
	while (EOF != fscanf(ifp, "%lu", &i)) {
		amount++;
		if (EOF == fscanf(ifp, "\n")) break;
	}

	assert_other(amount == 0, "Empty input file");
	nums = (unsigned long*)calloc(amount, sizeof(unsigned long));
	
	
	/* populate the array of nums */
	rewind(ifp);
	while(EOF != fscanf(ifp, "%lu", &i)) {
		nums[ind++] = i;
		if (EOF == fscanf(ifp, "\n")) break;
	}

	/* closing the file */
	fclose(ifp);
}




/* print an array of unsigned longs separated by new lines */
void printarr(unsigned long nums[], int amount) {
	int i;
	for (i=0; i < amount - 1; i++) {
		printf("%lu\n", nums[i]);
	}
	printf("%lu", nums[i]);
}



/* freeing all resources used by the program (when done with sorting) */
void free_program() {
	if (nums != NULL) {
		free(nums);
	}
}




/* personally preferred assert func */
void assert_other(int cond, char stderr[]) {
	if (cond) {
		printf("%s\n", stderr);
		exit(1);
	}
}
/******************************************************/
