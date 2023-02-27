# Maxwell Cole
# This program takes geometry data from the h5 file and 
# plots it to visualize the partitions

import numpy as np
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


filename = '/work/maxwell/vascpp/build/test.h5'
# numDims = 2

# retrieve data from H5 file
k = h5py.File(filename, 'r')
geom_array = k['GEOM_GROUP/GEOM_ARRAY']
sz = len(geom_array)
node_array = k['GEOM_GROUP/NODE_ARRAY']
bbox = k['/boundingBox']

#conductance_array = k['FLOW_GROUP/HEALTHY_CONDUCTANCE_ARRAY']

#print(list(k.attrs.items()))

# create figure
fig = plt.figure()
fig.set_size_inches(6.5, 5)
ax = fig.add_subplot(111)

# plot data for distributed network (i.e., after step2)
#e = 'VESSELS_PER_PARTITION_ARRAY' in k
e = 0

if e:
	# retreive partition data
	partition_num = k['VESSELS_PER_PARTITION_ARRAY']
	p_sz = len(partition_num)
	cmap = plt.get_cmap('brg',p_sz)
 
 	# plot the geom array
	for i in range(0, sz):
		[x1, x2] = [geom_array[i][0], geom_array[i][3]]
		[y1, y2] = [geom_array[i][1], geom_array[i][4]]
		node1 = node_array[i][0]
		node2 = node_array[i][1]
		norm_rad = 1.5*(geom_array[i][6] / geom_array[0][6])
		
		#
		for n in range(0, p_sz):
			if i < partition_num[0]:
				col = 0
			elif i >= sum(partition_num[:(-1)]):
				col = p_sz
			elif sum(partition_num[:n]) <= i < sum(partition_num[:n+1]):
				col = n
				
		# plot vessels, nodes, and vessel numbers
		ax.plot([x1,x2], [y1,y2], c=cmap(col), linewidth=norm_rad)
		plt.text(x1, y1+.001, node1, fontsize=8)
		plt.text(x2, y2+.001, node2, fontsize=8)
		plt.text((x1+x2)/2, (y1+y2)/2 + 0.001, i, fontsize=8)
  
	# Plot properties/customizations
	# change properties depending on # of nodes
	if p_sz == 2:
		vm = 2
		tk = [0.5, 1.5]
		tklbl = ['1', '2']
	elif p_sz == 4:
		vm = 4
		tk = [0.5, 1.5, 2.5, 3.5] 
		tklbl = ['1', '2', '3', '4']
	norm = mpl.colors.Normalize(vmin=0,vmax=vm)
	sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
	sm.set_array([])
	cb = plt.colorbar(sm, ticks=tk, ax=ax)
	cb.set_label(label='Partition Containing Geometry Data', weight='bold')
	cb.ax.set_yticklabels(tklbl)
	#ttl = filename[11] + ' Layers on ' + str(p_sz) + ' Nodes, After NP'
	
  
#! For single node
else:
	for i in range(0, sz):
		[x1, x2] = [geom_array[i][0], geom_array[i][3]]
		[y1, y2] = [geom_array[i][1], geom_array[i][4]]
		node1 = node_array[i][0]
		node2 = node_array[i][1]
		#conductance = conductance_array[i][0]
		norm_rad = 3*(geom_array[i][6] / geom_array[0][6])
		ax.plot([x1,x2], [y1,y2], color='black', linewidth=norm_rad)
		plt.text(x1, y1, node1, fontsize=6, bbox = dict(boxstyle=f"circle", fc="white"))
		plt.text(x2, y2, node2, fontsize=6, bbox = dict(boxstyle=f"circle", fc="white"))
		# plt.text((x1+x2)/2, (y1+y2)/2 + 0.001, i, fontsize=7, weight='bold')
		#plt.text((x1+x2)/2, (y1+y2)/2 + 0.001, conductance, fontsize=6, weight='bold')
		plt.xticks(fontsize=6)
		plt.yticks(fontsize=6)
		plt.xlim([bbox[0]-.01, bbox[1]+.01])
		plt.ylim([bbox[2]-.01, bbox[3]+.01])

ax.set_xlabel('x [mm]', fontsize=8)
ax.set_ylabel('y [mm]', fontsize=8)
#plt.title("5 gen on 2 nodes - vascpp", weight='bold')
plt.savefig('5-2-geometry.png')