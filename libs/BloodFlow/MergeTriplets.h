#pragma once

#include <vector>
#include "Triplet.h"

/**
 * function to merge triplets in a list.
 *
 * No Additional datastructre required, however does erase triplets
 *
 * @param[in,out] tripletArray
 */
void MergeTriplets(std::vector<Triplet>& tripletArray);

/**
 * function to merge triplets in a list.
 *
 * Uses a list of elements as intermediate structure, then copies data back to triplet array
 *
 * @param[in,out] tripletArray
 */
void MergeTripletsInPlace(std::vector<Triplet>& tripletArray);