
#include "EdgeFinderVessel.h"
#include "EdgeNode.h"

EdgeFinderVessel::EdgeFinderVessel() {}

EdgeFinderVessel::EdgeFinderVessel(long long startNodeIn, long long endNodeIn)
  : startNode(startNodeIn)
  , endNode(endNodeIn)
{}

EdgeFinderVessel::EdgeFinderVessel(
    long long startNodeIn, long long endNodeIn, short partitionNumIn)
  : startNode(startNodeIn)
  , endNode(endNodeIn)
  , partitionNum(partitionNumIn)
{}

EdgeFinderVessel::~EdgeFinderVessel() {}

std::vector<EdgeNode> EdgeFinderVessel::getArrayOfNodes() const
{
    std::vector<EdgeNode> nodes;
    if (!getStartNodeIsBcNode())
    {
        nodes.emplace_back(getStartNode());
    }
    if (!getEndNodeisBcNode())
    {
        nodes.emplace_back(getEndNode());
    }
    return nodes;
}

long long EdgeFinderVessel::FindMatchingNode(const EdgeFinderVessel& in) const
{
    //std::cout << "FindMatchingNode This Vessel:\n" << *this << "\ntestVessel:\n " << in<< std::endl;
    if ((startNode == in.getStartNode()) || (startNode == in.getEndNode()))
    {
        return startNode;
    }
    else
        return endNode;

    throw "no matching node found";
}

long long EdgeFinderVessel::FindMatchingNode(
    const std::array<long long, 2>& in) const
{
    if ((startNode == in[0]) || (startNode == in[1]))
    {
        return startNode;
    }
    else
        return endNode;
}

void EdgeFinderVessel::setConnectedVessels(std::vector<long long>& input)
{
    connectedVessels.assign(input.begin(), input.end());
}

bool EdgeFinderVessel::operator==(const EdgeFinderVessel& rhs) const
{
    return (startNode == rhs.getStartNode()) && (endNode == rhs.getEndNode());
}

std::ostream& operator<<(std::ostream& os, const EdgeFinderVessel& vessel)
{
    os << vessel.getStartNode() << " " << vessel.getEndNode() << " "
       << "\t";
    std::vector<long long> temp = vessel.getConnectedVessels();
    for (int k = 0; k < temp.size(); k++)
    {
        os << temp[k] << " ";
    }
    return os;
}