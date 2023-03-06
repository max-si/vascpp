#pragma once

#include <array>
#include <fstream>

#include "Node.h"
#include "VesselGenerator/coordinate.h"
//#include <vascular/core/Geometry.h>

#define SMALL_NUM 0.000000001

/**
 * \class BloodFlowVessel
 */
class BloodFlowVessel
{
    friend bool operator==(const BloodFlowVessel& lhs, const BloodFlowVessel& rhs);
    friend std::ostream& operator<<(std::ostream& os, const BloodFlowVessel& vessel);

public:
    BloodFlowVessel();
    BloodFlowVessel(double length, double radius, long long node1, long long node2);
    BloodFlowVessel(double length, double radius, double damagedRadii, long long node1, long long node2);
    ~BloodFlowVessel();
    /**
	 * Calculates the conductance of the vessel with the original radius
	 * @return conductance in mm^3/s
	 */
    double calculateHealthyConductance() const;
    /**
     * Calculates the conductance of the vessel with the damaged radius
     * @return conductance in mm^3/s
     */
    double calculateDamagedConductance() const;
    /**
	 * Scales the damaged radius by the supplied value
	 * @param value Scaling factor of radius
	 */
    void scaleDamagedRadii(double value);

    //Setters and Getters
    double getHealthyRadius() const
    {
        return healthyRadius;
    };
    double getDamagedRadius() const
    {
        return damagedRadius;
    };
    void setDamagedRadius(double val)
    {
        damagedRadius = val;
    };

    double GetHealthyFlow() const
    {
        return healthyFlow;
    }
    void SetHealthyFlow(double val)
    {
        healthyFlow = val;
    }
    double GetDamagedFlow() const
    {
        return damagedFlow;
    }
    void SetDamagedFlow(double val)
    {
        damagedFlow = val;
    }

    long long getNode1Index() const
    {
        return node1.GetIndex();
    }
    void setNode1Index(long long val)
    {
        node1.SetIndex(val);
    }
    long long getNode2Index() const
    {
        return node2.GetIndex();
    }
    void setNode2Index(long long val)
    {
        node2.SetIndex(val);
    }

    bool getNode1IsBcNode() const
    {
        return node1.GetIsBcNode();
    }
    bool getNode2IsBcNode() const
    {
        return node2.GetIsBcNode();
    }
    void setNode1IsBcNode(bool val)
    {
        node1.SetIsBcNode(val);
    }
    void setNode2IsBcNode(bool val)
    {
        node2.SetIsBcNode(val);
    }

    bool getNode1IsLocal() const
    {
        return node1.GetIsLocal();
    }
    bool getNode2IsLocal() const
    {
        return node2.GetIsLocal();
    }
    void setNode1IsLocal(bool val)
    {
        node1.SetIsLocal(val);
    }
    void setNode2IsLocal(bool val)
    {
        node2.SetIsLocal(val);
    }

    Node& getNode1()
    {
        return node1;
    }
    Node& getNode2()
    {
        return node2;
    }
    const Node& getNode1() const
    {
        return node1;
    }
    const Node& getNode2() const
    {
        return node2;
    }

    //void setDamagedRadiusArray(std::array<float, NUM_VESSEL_SECTIONS>& in);
    //void damageSpecificSegment(short segmentIndex, double value);

protected:
    //std::array<float, NUM_VESSEL_SECTIONS> damagedRadiusArray{};
    Node node1;
    Node node2;
    float healthyFlow;
    float damagedFlow;
    float vesselLength;
    float healthyRadius;
    float damagedRadius;
};

std::ostream& operator<<(std::ostream& os, const BloodFlowVessel& vessel);
double CalculateViscostity(double radius);
double CalculateViscostityInVivo(double radius);
double CalculateConductance(double radius, double length);

inline bool operator==(const BloodFlowVessel& lhs, const BloodFlowVessel& rhs)
{
    bool val1, val2, val3, val4, val5;
    val3 = fabs(lhs.getHealthyRadius() - rhs.getHealthyRadius()) <= SMALL_NUM;
    val4 = lhs.node1 == rhs.node1;
    val5 = lhs.node2 == rhs.node2;
    val1 = fabs(lhs.vesselLength - rhs.vesselLength) <= SMALL_NUM;
    val2 = fabs(lhs.getDamagedRadius() - rhs.getDamagedRadius()) <= SMALL_NUM;
    ;

    return val1 && val2 && val3 && val4 && val5;
}

double CalculateLengthOfVessel(Coordinate pt1, Coordinate pt2);
double CalculateLengthOfVessel(double x1, double y1, double z1, double x2, double y2, double z2);