#pragma once

// Header file for coordinate class
// Coordinate class stores (x,y,z) Cartesian coords for vessel start and end

#include <iostream>
#include <cmath>

#define SMALL_NUM 0.000000001

class Coordinate
{
//private:
	// double x, y, z;	
public:
	double x, y, z;
	Coordinate() : x(0), y(0), z(0) {}
	Coordinate(double x, double y, double z)
	: x(x), y(y), z(z)
	{}
	void print_coordinate();
	void reset(double x1, double y1, double z1);
	double coordinate_x();
	double coordinate_y();
	double coordinate_z();
	Coordinate add(double x1, double y1, double z1);

	bool operator==(const Coordinate& rhs) const; 
};