"""
Created on Mon Mar 14 13:48:09 2016

@author: emilyfairfax
with help from the LandLab tutorial for overland flow
and Greg Tucker
"""
#Imports
import numpy as np
from landlab import RasterModelGrid
import matplotlib.pylab as plt

"""Initialize"""
# Constants
num_rows = 20 # number of rows
num_cols = 100 # number of columns
dx = 10 # grid spacing
mx = 0.05 # slope of rock in x
my = 0.05# slope of rock in y
topomax = 100 # maximum elevation
dt = 0.1 # minutes
n = 0.05 # manning coefficient
R = 0.0003 # rainfall
I = 0 # infiltration

#Make the Grid
mg = RasterModelGrid(num_rows, num_cols, dx)
core_nodes = mg.core_nodes

# Create Variable Arrays
ground_elevation = mg.add_empty('node', 'land_surface_elevation')
water_elevation = mg.add_empty('node','water_surface_elevation')
H = mg.add_zeros('node','water_height')
Q = mg.add_zeros('link','flux')
dQdA = np.zeros(mg.number_of_links)

# Set Elevations
ground_elevation[:] = topomax - mx*mg.node_x
ground_elevation +=  my * np.abs(mg.node_y - mg.dx * ((num_rows - 1) / 2.0))
water_elevation[:] = ground_elevation+H

#Boundary Conditions - only open at downslope side
mg.set_closed_boundaries_at_grid_edges(False, True, True, True)


#Calculate Surface Slopes
all_gradients = mg.calculate_gradients_at_links(ground_elevation)
gradient = np.zeros(mg.number_of_links)
gradient[mg.active_links] = all_gradients[mg.active_links]
absgradient = np.abs(gradient) #taking sq rt later so needs to be abs value

#Make sure z and gradient look legit
print 'z', ground_elevation
print 'grad', gradient

"""Run"""
#Time Loop
for i in range(50000):
    # Calculate H at edges/links
    Hedge = mg.map_value_at_max_node_to_link('water_surface_elevation','water_height')
    #Calculate Flux
    Q = -np.sign(gradient)*1/n*((Hedge**(5/3)))*(absgradient**(1/2))
    #Calculate Flux Divergence
    dQdA = mg.calculate_flux_divergence_at_nodes(Q[mg.active_links])
    #Balance inputs and outputs
    dHdt = R-dQdA                 
    #Update elevations
    H[core_nodes] = np.maximum(H[core_nodes]+dHdt[core_nodes]*dt,0)
    water_elevation[core_nodes] = ground_elevation[core_nodes] + H[core_nodes]

"""Finalize"""
#plotting from LandLab Overland Flow Tutorial
# get a 2D array version of the water height
zwR = mg.node_vector_to_raster(ground_elevation)
# create raster image
plt.close()
image_extent = [0,num_cols*dx,0,num_rows*dx]
im = plt.imshow(zwR,cmap=plt.cm.RdBu,extent=image_extent)
plt.xlabel('Distance (m)', fontsize=12)
plt.ylabel('Distance (m)', fontsize=12)
# create contours
#cset = plt.contour(zwR,extent=image_extent)
#plt.clabel(cset,inline=True, fmt='%1.1f', fontsize=10)
# color bar on side
cb = plt.colorbar(im)
cb.set_label('Elevation (m)', fontsize=12)
# add title
plt.title('Channel Profile')
# Save plot
plt.savefig('Channel_Profile')
# Display plot
plt.show()


  
# get 2D array for contours
Hr = mg.node_vector_to_raster(H)

# Water height plot
image_extent = [0,num_cols*dx,0,num_rows*dx]
im = plt.imshow(Hr,cmap=plt.cm.RdBu,extent=image_extent)
plt.xlabel('Distance (m)', fontsize=12)
plt.ylabel('Distance (m)', fontsize=12)
# create contours
#cset = plt.contour(Hr,extent=image_extent)
#plt.clabel(cset,inline=True, fmt='%1.1f', fontsize=10)
# color bar on side with label 
cb = plt.colorbar(im)
cb.set_label('Water Height (m)', fontsize=12)
# plot title
plt.title('Water Height')
# save plot
plt.savefig('Water_Height')
# Display plot
plt.show()

plt.figure(2)