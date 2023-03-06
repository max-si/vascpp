

#include "EdgeNode.h"

EdgeNode::EdgeNode() {}

EdgeNode::EdgeNode(long long nodeIndex)
  : node(nodeIndex) {}

EdgeNode::~EdgeNode() {}

bool EdgeNode::operator==(EdgeNode rhs) const
{
    return (this->node == rhs.node);
}

bool EdgeNode::operator<(EdgeNode rhs) const
{
    return (this->node < rhs.node);
}

std::ostream& operator<<(std::ostream& os, const EdgeNode& nodeObj)
{
    os << nodeObj.node;
    if (nodeObj.hostPartition >= 0)
    {
        os << " " << nodeObj.hostPartition;
    }
    std::set<long long>::iterator itr;
    if (nodeObj.localConnectedVessels.size() > 0)
    {
        os << ":\n\tLocal Connections: ";
        for (itr = nodeObj.localConnectedVessels.begin();
             itr != nodeObj.localConnectedVessels.end(); itr++)
            os << *itr << " ";
    }
    if (nodeObj.nonLocalConnectedVessels.size() > 0)
    {
        os << "\n\tNon-Local Connections: ";
        for (itr = nodeObj.nonLocalConnectedVessels.begin();
             itr != nodeObj.nonLocalConnectedVessels.end(); itr++)
            os << *itr << " ";
    }
    if (nodeObj.connectedPartitions.size() > 0)
    {
        os << "\n\tConnected Partitions: ";
        for (std::set<int>::iterator itr2 = nodeObj.connectedPartitions.begin();
             itr2 != nodeObj.connectedPartitions.end(); itr2++)
            os << *itr2 << " ";
    }

    return os;
}

bool EdgeNodesEquivilent(const EdgeNode& node1, const EdgeNode& node2)
{
    bool equivilent = true;
    equivilent &= (node1.node == node2.node);
    equivilent &= (node1.connectedPartitions == node2.connectedPartitions);
    equivilent &= (node1.localConnectedVessels == node2.localConnectedVessels);
    equivilent &= (node1.nonLocalConnectedVessels == node2.nonLocalConnectedVessels);

    return equivilent;
}