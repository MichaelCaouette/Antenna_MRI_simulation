# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 15:02:47 2022

@author: mcaoue2
"""
# Very usefull command to use for getting the last-not-printed error
# Just type _p() in the console 
import traceback
_p = traceback.print_last

# The design to test
from jefimenko_acsolver.model_concentric_loops import ConcentricLoops
# To study the rotating fied
from jefimenko_acsolver.rotating_field import RotatingVector
# To view the system
from jefimenko_acsolver.viewer_wire_field_3D import Show3DVectorField # To view the geometry
# To inspect the solution
from jefimenko_acsolver.imageviewer.images import  GUIimage # To view the 3D volume with personal script
# To calculate the field
from jefimenko_acsolver.field_lines import FieldLines

import affine_transfo as afifi


import nibabel as nib # To create a nifty file with the 3D volume

import numpy as np

# =============================================================================
# Get the reference Image to overal the loop on
# Especially to get the right FOV
# =============================================================================

def save_for_ITKsnap(volume, affine,
                     savename):
    """
    Explaination soon. 
    It is to save the volume in orther to meet the ITK-snap
    """
    # ITK-SNAP reverse the x and y from the numpy convention. 
    # So I need to swap them
    my3Dstuff = volume.transpose((1, 0, 2))
    
    # I should check the affine
    new_image = nib.Nifti1Image(my3Dstuff, 
                                affine=affine)
    nib.save(new_image, savename)

### Find the center of the MRI image

# ref_img = nib.load('affine_test_ori.nii')
ref_img = nib.load('NMT_v2.0_sym_fh_ORIGINAL.nii.gz')
# Get the translation to the center of FOV
img_data = ref_img .get_fdata()
vox_center = (np.array(img_data.shape) - 1) / 2.
ref_center = afifi.apply_affine(ref_img.affine, 
                          vox_center[0], # x and y are inverted in position from the shape. 
                          vox_center[1], 
                          vox_center[2])
# Get the position of the bounds of the image. 
print('WE ASSUM NO ROTATION !')
ref_Nx, ref_Ny, ref_Nz = np.array(img_data.shape)
ref_xmin, ref_ymin, ref_zmin = afifi.apply_affine(ref_img.affine, 
                                            0,  0, 0)
ref_xmax, ref_ymax, ref_zmax = afifi.apply_affine(ref_img.affine, 
                                            ref_Nx-1, ref_Ny-1, ref_Nz-1)
ref_dx = ref_xmax - ref_xmin
ref_dy = ref_ymax - ref_ymin
ref_dz = ref_zmax - ref_zmin
print('ref_dx = %.2f mm'%ref_dx)
print('ref_dy = %.2f mm'%ref_dy)
print('ref_dz = %.2f mm'%ref_dz)

# =============================================================================
# # Define the antennas, and the system
# =============================================================================

# Let's list a couple of axis on which we will check the B1 rotating component
# Axis of the bore hole, or the B0 field (for the rotating calculation)
axis_B0 = [0, 0, 1 ]

Rsphere = 0.045  # m, radius of the sphere on which lays the coil
R_loops = 0.020 # m, radius of each loop
frac_overlap = 0.3 # Rough estimate of the radius overlap between the layers
angle_layer = 2*np.arctan(R_loops*(1-frac_overlap)/Rsphere)

# Check if the clearance is good
# This is related to the distance of the layer loops with the top loop axis. 
x_constraint = Rsphere*np.sin(angle_layer) - R_loops*np.cos(angle_layer)
clearance = 2*x_constraint
print('')
print('Check-up the clearance')
print('Rsphere = %.2f mm'%(Rsphere*1e3))
print('R_loops = %.2f mm'%(R_loops*1e3))
print('angle_layer = %.1f degree'%(angle_layer*180/np.pi))
print('Clearane = %.2f mm '%(clearance*1e3))

self = ConcentricLoops(
                 Rsphere=Rsphere, # m, radius of the sphere on which lays the coil
                 list_Rcircle = [R_loops, R_loops],# m, radiues of each loop
                 list_Nloop = [6], # # We should verify the optimal number
                 list_angle_layer = [angle_layer], # This angle should be check
                 pos_center_sphere = [0, 0, 0], # Origin of the the sphere
                 n_vec = [0, -1, 0],  # Top loop direction
                 N_dscrt_pts_loop=85 )
my_loops = self.get_wires()


# Set the signal
# FOr now, the phase delays, amplitudes, frequencies

# Set the delays for optimizing the rotating component
zeros = np.zeros(len(my_loops)) # Just to have the right array dimension
my_freq = 297.2*1e6
my_delays = self.get_delays_optimal_rotation(axis_B0, my_freq)
my_freqs   = zeros + 297.2*1e6  # Same frequency in every loop
my_Is     =  zeros + 1 # The exact value doesn't matter for the overal shape. 

# =============================================================================
# # Define the spatial points to probe
# =============================================================================

# A volume
# Nx, Ny, Nz = 40, 40, 30 # Same ratio as reference image
# Nx, Ny, Nz = 80, 80, 60# Same ratio as reference image
Nx, Ny, Nz = 18, 19, 15# Same ratio as reference image

# This will be translated to match the center of the reference image
 # + self.n_vec[0]

xmin, xmax = -0.5*ref_dx*1e-3 , 0.5*ref_dx*1e-3
ymin, ymax = -0.5*ref_dy*1e-3 , 0.5*ref_dy*1e-3
zmin, zmax = -0.5*ref_dz*1e-3 , 0.5*ref_dz*1e-3



# Now everything should be done with the following variable.
x = np.linspace(xmin, xmax, Nx) 
y = np.linspace(ymin, ymax, Ny) 
z = np.linspace(zmin, zmax, Nz) 

x_probe, y_probe, z_probe = np.meshgrid(x, y, z)


# =============================================================================
# Take a look at look at the system
# =============================================================================
vv = Show3DVectorField()
# Show all the elements of the antennas
# Scale the color 
dd = np.array(my_delays)*my_freq
dd = dd - min(dd)
vv.add_list_wire(my_loops , list_color=dd, 
                 str_cmap='brg', label='My loop',
                 color_title='Time delay*freq (unitless))',
                 show_corner=False)
# Show the direction of the B0 field
vv.add_vector_field([0], [-0.03], [0],
                    [axis_B0[0]], [axis_B0[1]], [axis_B0[2]],
                    scale=0.02,
                    label='B0 direction', color='k')
# Check the probing boxe
vv.add_bonding_box(xmin, xmax, 
                   ymin, ymax, 
                   zmin, zmax)

# =============================================================================
# # Compute the field (only spatial component)
# =============================================================================
myB1 = FieldLines(my_loops,
                  my_delays,
                  my_Is,
                  my_freqs)

print('Probing the field into space...')
(list_vcos, 
 list_vsin) = myB1.get_field_spatial(x_probe, 
                                     y_probe,
                                     z_probe,
                                     print_progress=True) 
print('Done')

# =============================================================================
# # Compute the rotating component with respect to a couple of orientation
# =============================================================================


# define the Affine
# Get the original affine (in the current referential)
# IN mm ! nibabel works with millimeters
A_ori = afifi.get_affine(min(x)*1e3, max(x)*1e3, len(x),
                         min(y)*1e3, max(y)*1e3, len(y),
                         min(z)*1e3, max(z)*1e3, len(z))

# Get my center
vc = (np.array(x_probe.shape) - 1) / 2.
# ITK-Snap OR the numpy convention as x and y axis reversed in the shape
my_center = afifi.apply_affine(A_ori, vc[1], vc[0], vc[2])

vtrans = ref_center  - my_center

# Now we know the center of the image. So we will translate up to there !
A_trans = afifi.add_translation(A_ori, vtrans)


my_affine = A_trans

# Check where is our center with respect to the reference
my_xmin, my_ymin, my_zmin = afifi.apply_affine(my_affine, 0, 0, 0)
my_xmax, my_ymax, my_zmax = afifi.apply_affine(my_affine, Nx-1, Ny-1, Nz-1)
print()
print('Comparison of the FOV of the ref and field')
print('my_center = ', my_center)
print('spa_center = ', ref_center)
print()
print('ref_xmin, ref_xmax = ', ref_xmin, ref_xmax, ref_dx)
print('ref_xmin, ref_xmax = ', my_xmin, my_xmax, my_xmax-my_xmin)  
print()
print('ref_xmin, ref_ymax = ', ref_ymin, ref_ymax, ref_dy)
print('ref_xmin, ref_ymax = ', my_ymin, my_ymax, my_ymax-my_ymin)  
print()
print('ref_zmin, ref_zmax = ', ref_zmin, ref_zmax, ref_dz)
print('ref_zmin, ref_zmax = ', my_zmin, my_zmax, my_zmax-my_zmin)  
print()


r_calc = RotatingVector()
B_clock, _, B_antic, _ = r_calc.get_rotation_pure_osc(list_vcos,
                                                      list_vsin, 
                                                      axis_B0)
# Make sure that it is an array of float
B_clock = np.array(B_clock, dtype='float')
B_antic = np.array(B_antic, dtype='float')

# =============================================================================
# Create the filter map
# =============================================================================
# It is the volume that is 1 where the coil are and 0 otherwise. 
epsilon = 0.003 # m, radius of the wire
eps2 = epsilon**2 # Shortcut to reduce CPU time

print()
print('Computing the wire location...')
# For each point probed, we determined how close it is from the wire
wires_map = np.zeros(np.shape(x_probe)) # Start with zero everywhere
# Loop over each wire
for wire in my_loops:
    # Loop over each point defining the wire
    for pos_wire in wire:
        # Compute the distance between this point on the wire and the location
        # to probe
        xw, yw, zw = pos_wire
        # d2 is a meshgrid of the distance between this particile point on the 
        # wire and all the point probed on the grid. 
        d2 = (xw-x_probe)**2 + (yw-y_probe)**2 + (zw-z_probe)**2
        # Add the wire where they should be
        wires_map += np.where(d2<=eps2, 1, 0) # One is d2<=eps, zero otherwise
# Rescale to be one or zero
wires_map = np.where(wires_map>=1, 1, 0)
print('Computing the wire location...Done !')        


        

# =============================================================================
# Visualize that 
# =============================================================================
# Check the volume field
slicy = GUIimage()
slicy.show()
slicy.set_list_image(B_clock)
   
# =============================================================================
# Save as a nifty file 
# =============================================================================

# Save the nifty file
# Just a number to distinguish some results
index_version = 3


# List what we want to save
list_keyword = ['CLOCKY', 'ANTICLOCK', 
                'WIRE',
                'ClockyFiltered', 
                'AnticlockyFiltered']
list_volume = [B_clock, B_antic, wires_map, 
               B_clock*(1-wires_map),
               B_antic*(1-wires_map)]

for i in range(len(list_keyword)):
    
    vol = list_volume[i]
    shape = np.shape(vol)

    filename = 'B1_v%d_B0X%.3dY%.3dZ%.3d_%dx%dx%d_%s.nii'%(index_version,
                                                      axis_B0[0]*10,
                                                      axis_B0[1]*10,
                                                      axis_B0[2]*10,
                                                      shape[1],
                                                      shape[0],
                                                      shape[2],
                                                      list_keyword[i])   
    filename = 'sim results/'+filename  
    save_for_ITKsnap(list_volume[i], 
                     my_affine,
                     filename)    
























