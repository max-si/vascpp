
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>

#include "BloodFlowVessel.h"
#include "VesselGenerator/vessel.h"
#include "VesselGenerator/coordinate.h"

#define NUM_VESSEL_SECTIONS 10
#define SMALL_NUM 0.000000001

using std::array;
typedef array<double, NUM_VESSEL_SECTIONS> sectionArray;

BloodFlowVessel::BloodFlowVessel()
  : node1(-1ll)
  , node2(-1ll)
  , damagedFlow(0)
  , healthyFlow(0)
  , vesselLength(0)
  , healthyRadius(0)
  , damagedRadius(0)
{}

BloodFlowVessel::BloodFlowVessel(double length, double radius, long long node1In, long long node2In)
  : vesselLength(length)
  , damagedFlow(0)
  , healthyFlow(0)
  , healthyRadius(radius)
  , node1(node1In)
  , node2(node2In)
  , damagedRadius(radius)
{}

BloodFlowVessel::BloodFlowVessel(double length, double radius,
    double damagedRadiusIn, long long node1In, long long node2In)
  : vesselLength(length)
  , damagedFlow(0)
  , healthyFlow(0)
  , healthyRadius(radius)
  , node1(node1In)
  , node2(node2In)
  , damagedRadius(damagedRadiusIn)
{}

BloodFlowVessel::~BloodFlowVessel() {}

double BloodFlowVessel::calculateHealthyConductance() const
{
    //std::cout << vesselLength << std::endl;
    return CalculateConductance(healthyRadius, vesselLength);
}

inline double CalculateDamagedConductanceUtility(const std::array<float, NUM_VESSEL_SECTIONS>& radiusArray, double totalLength)
{
    double coeff = 8 * totalLength / M_PI / static_cast<float>(NUM_VESSEL_SECTIONS);
    double temp = 0;

    for (int i = 0; i < NUM_VESSEL_SECTIONS; ++i)
    {
        temp += CalculateViscostity(radiusArray[i]) / pow(radiusArray[i], 4);
    }

    return 1 / (coeff * temp);
}

double BloodFlowVessel::calculateDamagedConductance() const
{
    if (damagedRadius < 1e-10)
        return 0;
    /*double temp = 0;
	for (int i = 0; i < NUM_VESSEL_SECTIONS; ++i) {
		temp += 1 / CalculateConductance(damagedRadiusArray[i], vesselLength/static_cast<float>(NUM_VESSEL_SECTIONS));
		
	}*/
    //return  1.0f / temp;
    //return static_cast<double> (damagedRadiusArray.size()) / temp;
    return CalculateConductance(damagedRadius, vesselLength);
}

void BloodFlowVessel::scaleDamagedRadii(double value)
{
    /*sectionArray::iterator itr;
	for (itr = damagedRadiusArray.begin(); itr != damagedRadiusArray.end(); itr++)
		*itr *= value;*/

    /*std::transform(damagedRadiusArray.begin(), damagedRadiusArray.end(), damagedRadiusArray.begin(),
		std::bind1st(std::multiplies<double>(), value));*/

    damagedRadius *= value;
}

/*void BloodFlowVessel::setDamagedRadiusArray(std::array<float, NUM_VESSEL_SECTIONS>& in)
{
	damagedRadiusArray = in;
}*/

/*
void BloodFlowVessel::damageSpecificSegment(short segmentIndex, double value)
{
	if (segmentIndex<0 || segmentIndex >= NUM_VESSEL_SECTIONS)
	{
		std::cout << "Vessel Segment Index outside of range" << std::endl;
		throw;
	}
	damagedRadiusArray[segmentIndex] *= value;
}
*/

std::ostream& operator<<(std::ostream& os, const BloodFlowVessel& vessel)
{
    os << vessel.getHealthyRadius() << "\t" << vessel.getNode1Index() << "\t"
       << vessel.getNode2Index() << "\n\t";
    //std::copy(vessel.damagedRadiusArray.begin(), vessel.damagedRadiusArray.end(), std::ostream_iterator<double>(os, " "));

    return os;
}

/* Helper functions */
//function to calculate viscosity
double CalculateViscostityInVivo(double radius)
{
    double diameter = 2 * radius * 1e3;
    double temp = diameter / (diameter - 1.1);
    temp *= temp;
    double eta = 6 * exp(-0.085 * diameter) + 3.2 -
        2.44 * exp(-.06 * pow(diameter, .645));
    return 0.036 * temp * (1.0f + temp * (eta - 1.0f));
}

double CalculateViscostity(double radius)
{
    double diameter = 2 * radius * 1e3;
    double eta = 220 * exp(-1.3 * diameter) + 3.2 - 2.44 * exp(-.06 * pow(diameter, .645));
    return 0.036 * eta;
}

double CalculateConductance(double radius, double length)
{
    double mu = CalculateViscostity(radius);
    //std::cout << mu << std::endl;
    if (mu == 0)
        return 0.0;
    double temp = M_PI * pow(radius, 4) / (8 * mu * length);
    //std::cout <<"in Conductance function " << temp << std::endl;
    return temp;
}

double CalculateLengthOfVessel(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double deltaX = x2 - x1;
    double deltaY = y2 - y1;
    double deltaZ = z2 - z1;
    return sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
}

double CalculateLengthOfVessel(Coordinate pt1, Coordinate pt2)
{
    return CalculateLengthOfVessel(pt1.x, pt1.y, pt1.z, pt2.x, pt2.y, pt2.z);
}