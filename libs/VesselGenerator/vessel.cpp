
// This file contains the Vessel class and methods
// This class contains an instance of a vessel and its description
//		including index#, start & end coordinates, length, and radius
//		and functions to return these values, nodes, and the conductance
// TODO: Add comments from pyvascular

#include <cmath>
#include <vector>
#include <iostream>

#include "coordinate.h"
#include "vessel.h"

Vessel::Vessel(long long idx, Coordinate start, Coordinate end, double length, double radius)
	: idx(idx)
	, startingPoint(start)
	, endingPoint(end)
	, vesselLength(length)
	, vesselRadius(radius)
{}

Vessel::Vessel(Coordinate start, Coordinate end, double length, double radius)
	: startingPoint(start)
	, endingPoint(end)
	, vesselLength(length)
	, vesselRadius(radius)
{}

Vessel::Vessel(Coordinate start, Coordinate end, double radius, long long n1, long long n2)
	: startingPoint(start)
	, endingPoint(end)
	, vesselRadius(radius)
	, node1(n1)
	, node2(n2)
{}

Vessel::~Vessel() {}

void Vessel::add_nodes(long long node1, long long node2)
{
	nodes.push_back(node1);
	nodes.push_back(node2);
}

double Vessel::get_conductance()
{
	double diameter = 2 * vesselRadius * 1e3;
	double eta = 220 * exp(-1.3*diameter) + 3.2 - 2.44 * exp(-0.06 * pow(diameter, .645));
	double mu = 0.036 * eta;
	if (mu == 0) { return 0; }
	else {
		return M_PI * pow(vesselRadius, 4) / (8 * vesselLength * mu);
	}
}

long long Vessel::get_idx() const {
	return idx;
}

void Vessel::print_vessel()
{
	std::cout << "\nVessel " << std::endl;
	std::cout << "Radius = " << vesselRadius << std::endl;
	//std::cout << "Length = " << vesselLength << std::endl;
	std::cout << "Start Point: ";
	startingPoint.print_coordinate();
	std::cout << "End Point: ";
	endingPoint.print_coordinate();
	//std::cout << "Nodes: " << nodes[0] << ", " << nodes[1] << std::endl;
	std::cout << "Nodes: " << getNode1() << ", " << getNode2() << std::endl;
	//std::cout << "Conductance = " << get_conductance() << std::endl;
}

Coordinate Vessel::get_startingPoint() const {
	return startingPoint;
}

Coordinate Vessel::get_endingPoint() const {
	return endingPoint;
}

double Vessel::get_start_x() {
	return startingPoint.coordinate_x();
}

double Vessel::get_start_y() {
	return startingPoint.coordinate_y();
}

double Vessel::get_start_z() {
	return startingPoint.coordinate_z();
}

double Vessel::get_end_x() {
	return endingPoint.coordinate_x();
}

double Vessel::get_end_y() {
	return endingPoint.coordinate_y();
}

double Vessel::get_end_z() {
	return endingPoint.coordinate_z();
}

std::vector<long long> const& Vessel::get_nodes() const {
	return nodes;
}

long long Vessel::getNode1() const
{
    return node1;
}

long long Vessel::getNode2() const
{
    return node2;
}

double Vessel::get_radius() const {
	return vesselRadius;
}

double Vessel::calculate_length() {
	return sqrt(pow((get_end_x()-get_start_x()),2) + pow((get_end_y()-get_start_y()),2) + pow((get_end_z()-get_start_z()),2));
}

void Vessel::add_connectedVessel(long long input) {
	connectedVessels.push_back(input);
}

std::vector<long long> const& Vessel::get_connectedVessels() const {
	return connectedVessels;
}

std::vector<long long>::size_type Vessel::get_numConnectedVessels() const {
	return connectedVessels.size();
}

bool CompareIDX(const Vessel& left, const Vessel& right) {
	if (left.get_idx() < right.get_idx()) {
		return true;
	} else {
		return false;
	}
}

bool CompareEndNodes(const Vessel& left, const Vessel& right)
{
	long long leftNode = left.get_nodes()[1];
    long long rightNode = right.get_nodes()[1];

    if (leftNode < rightNode)
    {
        return true;
    }
    else if (leftNode > rightNode)
    {
        return false;
    }
	else if (left.get_nodes()[0] < right.get_nodes()[0])
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool CompareEndNodes2(const Vessel& left, const Vessel& right)
{
    long long leftNode = left.getNode2();
    long long rightNode = right.getNode2();

    if (leftNode < rightNode)
    {
        return true;
    }
    else if (leftNode > rightNode)
    {
        return false;
    }
    else if (left.getNode1() < right.getNode1())
    {
        return true;
    }
    else
    {
        return false;
    }
}