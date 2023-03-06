#pragma once

#include <array>
#include <iostream>
#include <vector>

#include "EdgeNode.h"

/**
 * \class EdgeFinderVessel
 *
 * This class is a specialized vessel object for locating edge junctions.
 * It includes methods to help simplify the implementation of the finder algorithm
 * Could have probably simplified the interface by using a get node reference to allow access to the node props
 *
 */
class EdgeFinderVessel
{
public:
    EdgeFinderVessel();
    EdgeFinderVessel(long long startNodeIn, long long endNodeIn);
    EdgeFinderVessel(
        long long startNodeIn, long long endNodeIn, short partitionNumIn);
    ~EdgeFinderVessel();
    std::vector<EdgeNode> getArrayOfNodes() const;
    long long FindMatchingNode(const EdgeFinderVessel& in) const;
    long long FindMatchingNode(const std::array<long long, 2>& in) const;
    void setConnectedVessels(std::vector<long long>& input);
    const std::vector<long long>& getConnectedVessels() const
    {
        return connectedVessels;
    };
    std::vector<long long>& getConnectedVessels()
    {
        return connectedVessels;
    };
    long long getStartNode() const
    {
        return startNode;
    }
    void setStartNode(long long val)
    {
        startNode = val;
    }
    long long getEndNode() const
    {
        return endNode;
    }
    void setEndNode(long long val)
    {
        endNode = val;
    }
    bool getVesselLocal() const
    {
        return vesselLocal;
    }
    void setVesselLocal(bool val)
    {
        vesselLocal = val;
    }
    bool getStartNodeIsBcNode() const
    {
        return startNodeIsBcNode;
    }
    void setStartNodeIsBcNode(bool val)
    {
        startNodeIsBcNode = val;
    }
    bool getEndNodeisBcNode() const
    {
        return endNodeisBcNode;
    }
    void setEndNodeisBcNode(bool val)
    {
        endNodeisBcNode = val;
    }
    bool operator==(const EdgeFinderVessel& rhs) const;

protected:
    long long startNode;
    long long endNode;
    bool startNodeIsBcNode = false;
    bool endNodeisBcNode = false;
    bool vesselLocal = false;
    short partitionNum;
    std::vector<long long> connectedVessels;

private:
    friend std::ostream& operator<<(
        std::ostream& os, const EdgeFinderVessel& vessel);
};