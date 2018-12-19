#ifndef NEIGHBOR2D_H_INCLUDED
#define NEIGHBOR2D_H_INCLUDED

#include <iostream>
#include <stdio.h>
using namespace std;

template <class T, size_t R>
class Neighbor2D{ 
	public:
		vector<T> neighborVoxels2D;
		vector<vector<int>> index2D;
		vector<T> neighborVoxels3D;
		vector<vector<int>> index3D;

};

#endif