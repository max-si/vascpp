#include <iostream>
#include <iterator>
#include <list>
#include <vector>

#include "MergeTriplets.h"
#include "Triplet.h"

using std::vector;

typedef vector<Triplet> TripletArray;

void MergeTriplets(std::vector<Triplet>& tripletArray)
{
    TripletArray::iterator itr1 = tripletArray.begin();
    TripletArray::iterator nextElementIter = tripletArray.begin();
    //TripletArray::size_type index = 0;
    while (itr1 != tripletArray.end())
    {
        TripletArray::iterator itr2 = itr1;
        itr1++;
        if (itr1 != tripletArray.end())
        {
            while (itr1 != tripletArray.end() && *itr1 == *itr2)
            {
                *itr2 += *itr1;
                itr1++;
            }
            std::swap(*nextElementIter, *itr2);
            ++nextElementIter;
        }
        else
        {
            *nextElementIter = *itr2;
            ++nextElementIter;
        }
    }
    tripletArray.erase(nextElementIter, tripletArray.end());
    //std::cout << "Merge Triplet Index: " << index << std::endl;

    return;
}

void MergeTripletsInPlace(std::vector<Triplet>& tripletArray)
{
    if (tripletArray.empty())
        return;
    std::list<Triplet> myList;
    myList.assign(tripletArray.begin(), tripletArray.end());
    std::list<Triplet>::iterator itr1, itr2;
    itr1 = myList.begin();
    while (itr1 != myList.end())
    {
        itr2 = itr1;
        itr1++;
        if (itr1 != myList.end())
        {
            while (itr1 != myList.end() && *itr1 == *itr2)
            {
                *itr2 += *itr1;
                itr1 = myList.erase(itr1);
            }
        }
    }
    tripletArray.assign(myList.begin(), myList.end());
    return;
}