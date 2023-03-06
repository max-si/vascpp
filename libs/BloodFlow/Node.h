#pragma once

/**
 * Class that describes a network junction.
 *
 * Lets user probe whether anode is a local or a boundary condtion node.
 * Could probably be implemented with out the over head in the class and instead be handled by a set that is quickly
 * searchable.
 *
 */
class Node
{
public:
    Node();
    Node(long long index);
    //~Node();

    long long GetIndex() const
    {
        return index;
    }
    void SetIndex(long long val)
    {
        index = val;
    }

    bool GetIsBcNode() const
    {
        return isBcNode;
    }
    void SetIsBcNode(bool val)
    {
        isBcNode = val;
    }

    bool GetIsLocal() const
    {
        return isLocal;
    }
    void SetIsLocal(bool val)
    {
        isLocal = val;
    }

    bool operator==(const Node& rhs) const;

private:
    long long index = -1;
    bool isBcNode = false;
    bool isLocal = false;
};