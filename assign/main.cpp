#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "PSRS.h"


/* global variable declarations */
uint_64 *nums;



/* auxiliary funcs declaration */
static void parse_args(FILE*, int*);
static void printarr(uint_64 *, int);
static void free_program();
/******************************/





int main(int argv, char* args[]) {
	int processors;
	int len = 0;

	/* parse arguments */
	char* str = (char*)"Incorrect amount of arguments!";
	assert_other(argv != 3, str); 
	processors = atoi(args[1]);
	FILE* ifp = fopen(args[2], "r");
	parse_args(ifp, &len);
	

	/* sort + measure time */
	struct timeval start, end;
	gettimeofday(&start, NULL);
	parl_PSRS(nums, processors, len);
	gettimeofday(&end, NULL);
	double elapsed = (end.tv_sec - start.tv_sec)*1000000.0 + (end.tv_usec - start.tv_usec);

	/* print to stdout */
	printf("PSRS (Parallel Sorting Regular Sampling): %.lf\n", elapsed);
      	printarr(nums, len);


	/* coherently finish the program */
	free_program();
}



/***************  GENERIC FUNCTIONS  *****************/



/* parsing crucial data of the input file */
static void parse_args(FILE* ifp, int* length) {
	uint_64 num;
	int ind = 0;	
	
	/* initialize the array of nums first */
	while (EOF != fscanf(ifp, "%lu", &num)) {
		(*length)++;
	}
	
	
	char* str = (char*)"Empty input file";
	assert_other((*length) == 0, str);
	nums = (uint_64*)calloc((*length), sizeof(uint_64));
	
	
	
	/* populate the array of nums */
	rewind(ifp);
	while(EOF != fscanf(ifp, "%lu", &num)) {
		nums[ind++] = num;
	}


	/* closing the file */
	fclose(ifp);
}




/* print an array of uint_64 integers separated by new lines */
static void printarr(uint_64 numbers[], int length) {
	int i;
	for (i=0; i < length; i++) {
		printf("%lu\n", numbers[i]);
	}
}



/* freeing all resources used by the program (when done with sorting) */
static void free_program() {
	if (nums != NULL) {
		free(nums);
	}
}




/* personally preferred assert func */
void assert_other(int cond, char strerr[]) {
	if (cond) {
		printf("%s\n", strerr);
		exit(1);
	}
}
/******************************************************/
