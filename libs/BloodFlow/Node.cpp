#include "Node.h"

Node::Node()
  : index(-1)
{}

Node::Node(long long indexIn)
  : index(indexIn)
{}

bool Node::operator==(const Node& rhs) const
{
    return (this->index == rhs.GetIndex());
}