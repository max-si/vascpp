#pragma once

#include <iostream>
#include <vector>

#include "VesselGenerator/vessel.h"
#include "VesselGenerator/coordinate.h"

// Vessel object used in the network partitioning algorithm
// Inherits from the vessel class from generator algorithm

class PreprocessorVessel : public Vessel
{
public:
	PreprocessorVessel();
	PreprocessorVessel(Coordinate start, Coordinate end, double radius, long long n1, long long n2) 
		: Vessel(start, end, radius, n1, n2)
		, partitionNumber(-1)
		, initialArrayPosition(-1) {};

	int getPartitionNumber() const { return partitionNumber; }

	void setPartitionNumber(int val) 
	{
		partitionNumber = val;
	}

	void setConnectedVessels(std::vector<long long>& input);

	void addConnectedVessel(long long input)
	{
		connectedVessels.push_back(input);
	}

	void assignConnectedVessels(std::vector<long long>::iterator& start, std::vector<long long>::iterator& end);

	std::vector<long long> const& getConnectedVessels() const 
	{
		return connectedVessels;
	}

	std::vector<long long>& getConnectedVessels()
	{
		return connectedVessels;
	}

	std::vector<long long>::size_type getNumConnectedVessels() const
	{
		return connectedVessels.size();
	}

	long long getInitialArrayPosition() const
	{
		return initialArrayPosition;
	}

	void setInitialArrayPosition(size_t val)
	{
		initialArrayPosition = val;
	}
	bool operator==(PreprocessorVessel& rhs);
	bool operator!=(PreprocessorVessel& rhs);

protected:
	std::vector<long long> connectedVessels;
	long long initialArrayPosition;
	int partitionNumber;
};

bool CompareVesselPartition(const PreprocessorVessel& lhs, const PreprocessorVessel& rhs);

bool CompareVesselNodesPartition(const PreprocessorVessel& lhs, const PreprocessorVessel& rhs);

std::ostream& operator<<(std::ostream& os, const PreprocessorVessel& vessel);