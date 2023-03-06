#pragma once

#include <iostream>
#include <set>

/**
 * \class EdgeNode
 *
 * \breif Class to collect data on non-local nodes in a network.
 *
 * This class contains members and functions to help determine if a node is a partition edge or not.
 * The data stored in this class is partially serialized for use in the blood flow algorithm.
 * The nodes identified as an edge node are communication points in the blood flow algorithm.
 */
class EdgeNode
{
public:
    EdgeNode();
    EdgeNode(long long nodeIndex);
    ~EdgeNode();

    long long getNode() const
    {
        return node;
    }

    void setNode(long long val)
    {
        node = val;
    }

    const std::set<long long>& getLocalConnectedVessels() const
    {
        return localConnectedVessels;
    }

    std::set<long long>& getLocalConnectedVessels()
    {
        return localConnectedVessels;
    }

    const std::set<int>& getConnectedPartitions() const
    {
        return connectedPartitions;
    }

    int getNumLocalConnectedVessels() const
    {
        return localConnectedVessels.size();
    }

    std::set<int>& getConnectedPartitions()
    {
        return connectedPartitions;
    }

    std::set<long long>& getNonLocalConnectedVessels()
    {
        return nonLocalConnectedVessels;
    }

    const std::set<long long>& getNonLocalConnectedVessels() const
    {
        return nonLocalConnectedVessels;
    }

    void addConnectedPartition(int partNum)
    {
        connectedPartitions.insert(partNum);
    }

    void addNonLocalConnectedVessel(long long vesselIndex)
    {
        nonLocalConnectedVessels.insert(vesselIndex);
    }

    void addLocalConnectedVessel(long long vesselIndex)
    {
        localConnectedVessels.insert(vesselIndex);
    }

    int getHostPartition() const
    {
        return hostPartition;
    }
	
    void setHostPartition(int val)
    {
        hostPartition = val;
    }

    bool operator==(EdgeNode rhs) const;
    bool operator<(EdgeNode rhs) const;

protected:
    friend std::ostream& operator<<(std::ostream& os, const EdgeNode& nodeObj);
    friend bool EdgeNodesEquivilent(
        const EdgeNode& node1, const EdgeNode& node2);
    long long node = -1;
    int hostPartition = -1;
    std::set<int> connectedPartitions;
    std::set<long long> localConnectedVessels;
    std::set<long long> nonLocalConnectedVessels;
};