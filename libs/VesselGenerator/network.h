#pragma once

// Header file for Network class
// Includes constructors for class methods and constants

#include <cmath>
#include <vector>

#include "vessel.h"

typedef std::vector<Vessel> VesselVector;

class Network
{
public:
	Network();
	Network(int numLevels, int numDimensions);
	~Network();
	
	// a lot of get methods
	int getNumLevels()
	{
		return numberOfLevels;
	}
	int getNumDimensions()
	{
		return dims;
	}
	double getInitialBifurcationAngle()
	{
		return initialBifurcationAngle;
	}
	double getBifurcationReductionRatio()
	{
		return bifurcationReductionRatio;
	}
	double getLengthReductionRatio()
	{
		return lengthReductionRatio;
	}
	double getRadiusReductionRatio()
	{
		return radiusReductionRatio;
	}
	double getMaxVesselLength()
	{
		return maxVesselLength;
	}
	double getMaxVesselRadius()
	{
		return maxVesselRadius;
	}
	double getMinVesselLength()
	{
		return minVesselLength;
	}
	double getMinVesselRadius()
	{
		return minVesselRadius;
	}

	int getNumVessels()
	{
		return numVessels;
	}

	int getNumNodes()
	{
		return numNodes;
	}

	void generate_artery_root(VesselVector& vessels, double x_projection, double y_projection);
	void generate_artery_body(VesselVector& vessels, int level, double x_projection, double y_projection, double length, double radius);
	void generate_artery(VesselVector& vessels);
	void generate_vein_root(VesselVector& vessels, double x_projection, double y_projection);
	void generate_vein_body(VesselVector& vessels, int level, double x_projection, double y_projection, double length, double radius);
	void generate_vein(VesselVector& vessels);
	//! ****************
	void generate_parallel(std::string filename);
	void generate_root_parallel(VesselVector& arteries, VesselVector& veins, int rootNumLevels);
	void generate_sub_parallel(VesselVector& arteries, VesselVector& veins, int rootNumLevels, int subNumLevels);
	//! ****************
	// VesselVector generate(std::string filename);
	void generate(std::string filename);
	void print_vessels(VesselVector& vessels);
	double getOutputFlowRate();
	std::vector<double> getBoundingBox();


protected:
	int dims;
	int numberOfLevels;
	double initialBifurcationAngle;
	double bifurcationReductionRatio;
	double maxVesselLength;
	double maxVesselRadius;
	const int numNodes = 3 * pow(2, numberOfLevels - 1);
	const int numVessels = pow(2, numberOfLevels + 1) - 2;
	const double lengthReductionRatio = 0.8;
	const double radiusReductionRatio = 0.793700525984100;
	const double minVesselRadius = 0.0025;
	const double minVesselLength = 0.055;
	std::vector<double> boundingBox;
};