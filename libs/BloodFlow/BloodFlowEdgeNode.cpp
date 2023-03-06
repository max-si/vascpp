
#include "BloodFlowEdgeNode.h"

BloodFlowEdgeNode::BloodFlowEdgeNode() {}

BloodFlowEdgeNode::BloodFlowEdgeNode(long long nodeIndex)
  : mNodeIndex(nodeIndex)
{}

BloodFlowEdgeNode::~BloodFlowEdgeNode() {}

std::ostream& operator<<(std::ostream& os, const BloodFlowEdgeNode& data)
{
    os << data.mNodeIndex << ": " << data.mHostPartition << " - ";
    for (auto itr = data.mAssociatedPartitions.begin();
         itr != data.mAssociatedPartitions.end(); itr++)
        os << *itr << " ";

    return os;
}

bool BloodFlowEdgeNode::operator<(const BloodFlowEdgeNode& rhs) const
{
    return mNodeIndex < rhs.getNodeIndex();
}

void BloodFlowEdgeNode::addAssociatedPartition(short partitionNum)
{
    mAssociatedPartitions.insert(partitionNum);
}

bool BloodFlowEdgeNode::operator==(const BloodFlowEdgeNode& rhs) const
{
    return mNodeIndex == rhs.getNodeIndex();
}