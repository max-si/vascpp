
#include <iostream>

#include "Triplet.h"

Triplet::Triplet()
	: first(-1)
	, second(-1)
	, third(0)
{}

Triplet::Triplet(long long i, long long j, double val)
	: first(i)
	, second(j)
	, third(val)
{}

Triplet::~Triplet() {}

bool Triplet::operator<(const Triplet& rhs) const {
	return CompareTriplet(*this, rhs);
}

Triplet Triplet::operator+=(const Triplet& RHS) {
	third += RHS.third;
	return *this;
}

bool CompareTriplet(const Triplet& x, const Triplet& y) {
	if (x.first < y.first) {
		return true;
	} else if (x.first > y.first) {
		return false;
	} else if (x.second < y.second) {
		return true;
	} else {
		return false;
	}
}

std::ostream& operator<<(std::ostream& os, const Triplet& val) {
	os << val.first << ", " << val.second << ", " << val.third;
	return os;
}