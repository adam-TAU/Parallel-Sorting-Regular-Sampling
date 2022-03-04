





// C++ program to merge k sorted arrays
// of size n each.
#include "mergeKarrays.h"
#include <bits/stdc++.h>
#include <stdio.h>
#include<iostream>
typedef unsigned long long uint_64;

 
 using namespace std;
// A pair of pairs, first element is going to
// store value, second element index of array
// and third element index in the array.
typedef pair<uint_64, pair<int, int> > ppi;
 
// This function takes an array of arrays as an
// argument and all arrays are assumed to be
// sorted. It merges them together and prints
// the final sorted output.
void mergeKArrays(int* partition_borders, uint_64** local_parts_array, int parts, uint_64* start_copy, int thread_num)
{

	/* build arr (vector<vector<uint_64> > arr) */
	vector<vector<uint_64> > arr;


	/* Extract the partitions in each thread, that are supposed to be handled by this thread (so all of the <thread_num> partitions, of all threads, shall be merged into one, here */
	for(int i = 0; i < parts; i++) {
		int bottom, top, thread_partition_size, offset;

		// get the i-th thread's, <thread_num> partition borders. While the 0-th partition is starts at index 0.
		offset = i * (parts + 1) + thread_num; 
		bottom = partition_borders[offset];
		top = partition_borders[offset+1];
		thread_partition_size = top - bottom;

		// copy the partition (using its borders - relative to the i-th thread's local array) into the final position of the cummulative <thread_num> partition (0th partition starts at index 0)
		if(thread_partition_size > 0) {
			vector<uint_64> tmp(local_parts_array[i] + bottom, local_parts_array[i] + top);
			arr.push_back(tmp);
		}
		
	}



    vector<uint_64> output;
 
    // Create a min heap with k heap nodes. Every
    // heap node has first element of an array
    priority_queue<ppi, vector<ppi>, greater<ppi> > pq;
 
    for (int i = 0; i < arr.size(); i++)
        pq.push({ arr[i][0], { i, 0 } });
 
    // Now one by one get the minimum element
    // from min heap and replace it with next
    // element of its array
    while (pq.empty() == false) {
        ppi curr = pq.top();
        pq.pop();
 
        // i ==> Array Number
        // j ==> Index in the array number
        int i = curr.second.first;
        int j = curr.second.second;
 
        output.push_back(curr.first);
 
        // The next element belongs to same array as
        // current.
        if (j + 1 < arr[i].size())
            pq.push({ arr[i][j + 1], { i, j + 1 } });
    }
    
    
    copy(output.begin(), output.end(), start_copy);
 
}
