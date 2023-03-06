#pragma once

#include <iostream>
#include <set>

/**
 * \class BloodFlowEdgeNode
 *
 * \breif A class to describe edge nodes of a vascular network
 *
 * An edge node is any node that is present on more than one compute node after partitioing.
 * This class us used to thelp determine if a node is an edge node and specify the associated partitions
 *
 */
class BloodFlowEdgeNode
{
public:
    /// Default Constructor
    BloodFlowEdgeNode();
    /**
     * Constructor that initalizes the node index on creation
     * @param nodeIndex The index of the edge node
     */
    BloodFlowEdgeNode(long long nodeIndex);
    ///Destructor
    ~BloodFlowEdgeNode();

    friend std::ostream& operator<<(std::ostream& os, const BloodFlowEdgeNode& data);
    bool operator==(const BloodFlowEdgeNode& rhs) const;
    bool operator<(const BloodFlowEdgeNode& rhs) const;
    
	//Setters Getters
    long long getNodeIndex() const
    {
        return mNodeIndex;
    }
    void setNodeIndex(long long val)
    {
        mNodeIndex = val;
    }
    short getHostPartition() const
    {
        return mHostPartition;
    }
    void setHostPartition(short val)
    {
        mHostPartition = val;
    }
    /**
	 * Method that appends a partition on to the set
	 * @param partitionNum the number of the partition to add to set
	 */
    void addAssociatedPartition(short partitionNum);
    std::set<short>& getAssociatedPartitions()
    {
        return mAssociatedPartitions;
    }
    const std::set<short>& getAssociatedPartitions() const
    {
        return mAssociatedPartitions;
    }

private:
    long long mNodeIndex = -1;    ///< Index of the edge node
    short mHostPartition = -1;    ///< Assigned host Patition
    std::set<short> mAssociatedPartitions;    ///< Set of associated partitions for use identifying communication pattern
};