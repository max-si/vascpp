#pragma once

#include <iostream>

// class to store triplet that described a matrix element in COO format
class Triplet
{
public:
	long long first;
	long long second;
	double third;
	Triplet();
	Triplet(long long, long long, double);
	~Triplet();

	bool operator<(const Triplet&) const;

	// adding triplets results in only third item being summed
	Triplet operator +=(const Triplet&);

	bool operator==(const Triplet& rhs) const {
		return (first == rhs.first) && (second == rhs.second);
	}

	bool operator !=(const Triplet& rhs) const {
		return !(*this == rhs);
	}
};

inline Triplet operator+(const Triplet& LHS, const Triplet& RHS) {
	Triplet r = LHS;
	r += RHS;
	return r;
};

bool CompareTriplet(const Triplet&x, const Triplet& y);

std::ostream& operator<<(std::ostream& os, const Triplet& val);