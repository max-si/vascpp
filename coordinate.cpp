// This file contains the class methods for the Coordinate class
// Coordinate class stores (x,y,z) Cartesian coords for vessel start and end

#include <iostream>
#include "coordinate.h"

// Method: print the coordinate as x,y,z
void Coordinate::print_coordinate()
{
	std::cout << x << ", " << y << ", " << z << std::endl;
}

// Method: reset coordinate values 
void Coordinate::reset(double x1, double y1, double z1) 
{
	x = x1;
	y = y1;
	z = z1;
}

// Method: return x value of coordinate
double Coordinate::coordinate_x()
{
	return x;
}

// Method: return y value of coordinate
double Coordinate::coordinate_y()
{
	return y;
}

// Method: return z value of coordinate
double Coordinate::coordinate_z()
{
	return z;
}

// Method: returns new coordinate object with sum components of previous object
Coordinate Coordinate::add(double x1, double y1, double z1)
{
	double x2 = x + x1;
	double y2 = y + y1;
	double z2 = z + z1;
	return Coordinate(x2, y2, z2);
}