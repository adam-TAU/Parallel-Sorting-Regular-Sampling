





// C++ program to merge k sorted arrays
// of size n each.
#include "mergeKarrays.h"
#include <bits/stdc++.h>
#include <stdio.h>
#include<iostream>
typedef unsigned long uint_64;

 
 using namespace std;
// A pair of pairs, first element is going to
// store value, second element index of array
// and third element index in the array.
typedef pair<uint_64, pair<int, int> > ppi;
 
// This function takes an array of arrays as an
// argument and all arrays are assumed to be
// sorted. It merges them together and prints
// the final sorted output.
uint_64* mergeKArrays(uint_64** partitions, int parts, int* sizes)
{


	/* build arr (vector<vector<uint_64> > arr) */
	vector<vector<uint_64> > arr;
	
	for (int i=0; i < parts; i++) {
		vector<uint_64> v(partitions[i], partitions[i] + sizeof partitions[i] / sizeof partitions[i][0]);
		arr.push_back(v);
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
 
	return output.data();
}
