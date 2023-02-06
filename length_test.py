import numpy as np

numLevels = 34

# Starting Parameters
bifAngle = 33.31											#Initial bifurcation angle
bifAngRatio = 0.88											#Bifurcation Reduction Ratio
lengthRatio = 0.8											#Length reduction ratio
minLength = 1											#Length of smallest capillary, 55 um
radiusRatio = 2**(-1/3)										#Radius reduction ratio
minRadius = .008 / 2											#Radius of smallest capillary, 2.5 um
numVessels = 2**(numLevels+1) - 2							#Number of vessels in the system
numNodes = 1.5*(2**numLevels)								#Number of nodes in the system

# Calculate Static Parameters
maxLength = minLength/(lengthRatio**(numLevels-1))
maxRadius = minRadius/(radiusRatio**(numLevels-1))

# Define x, y, and z extents	
yExtent = 0
zExtent = 0
xExtent = maxLength
vesselLength = maxLength

totalHalfLength = maxLength
volume = np.pi * (maxRadius*maxRadius) * maxLength
	
# Calculate x, y, and z Extents
for i in range(1, numLevels):
	vesselLength = maxLength*(lengthRatio**i)
	vesselRadius = maxRadius * (radiusRatio**i)
	xExtent += vesselLength*np.cos((np.pi*bifAngle)/180)
	yExtent += vesselLength*np.sin((np.pi*bifAngle)/180)
	bifAngle = bifAngle*bifAngRatio
	# print("length = ", vesselLength)
	# print("rad = ", vesselRadius)
	totalHalfLength = totalHalfLength + vesselLength * (2**i)
	vesselVolume = np.pi * vesselLength * (vesselRadius*vesselRadius)
	volume = volume + vesselVolume * (2**i)
	

print("total Length = ", 2*totalHalfLength)
print("total volume = ", 2*volume/1e6)