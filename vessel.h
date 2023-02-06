#pragma once

// Header file containing the Vessel class constructor and methods
//TODO: create constructor for PreprocessorVessel/nodes, change nodes to long long, add operators
//TODO: or create new vessel class that inherets from this vessel class?

#include <vector>
#include "coordinate.h"

class Vessel
{
public:
	Vessel() : idx(0), startingPoint(Coordinate()), endingPoint(Coordinate()), vesselLength(0), vesselRadius(0) {}
	Vessel(long long idx, Coordinate start, Coordinate end, double length, double radius);
	Vessel(Coordinate start, Coordinate end, double length, double radius);
	Vessel(Coordinate start, Coordinate end, double radius, long long n1, long long n2);
	~Vessel();

	void add_nodes(long long node1, long long node2);
	long long get_idx() const;
	void print_vessel();
	double get_conductance();
	Coordinate get_startingPoint();
	Coordinate get_endingPoint();
	double get_start_x();
	double get_start_y();
	double get_start_z();
	double get_end_x();
	double get_end_y();
	double get_end_z();
	long long getNode1() const;
	long long getNode2() const;
	std::vector<long long> const& get_nodes() const;
	double get_radius();
	double calculate_length();
	void add_connectedVessel(long long input);
	std::vector<long long> const& get_connectedVessels() const;
	std::vector<long long>::size_type get_numConnectedVessels() const;

protected:
	long long idx;
	long long node1;
	long long node2;
	Coordinate startingPoint;
	Coordinate endingPoint;
	double vesselLength;
	double vesselRadius;
	std::vector<long long> nodes;
	std::vector<long long> connectedVessels;
};

bool CompareIDX(const Vessel& left, const Vessel& right);

bool CompareEndNodes(const Vessel& left, const Vessel& right);

bool CompareEndNodes2(const Vessel& left, const Vessel& right);