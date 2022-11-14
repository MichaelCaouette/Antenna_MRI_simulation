# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 15:13:56 2022

Goal:
    Test the modules !

@author: mcaoue2
"""



# Very usefull command to use for getting the last-not-printed error
# Just type _p() in the console 
import traceback
_p = traceback.print_last


from model_5_loops import Flower5loop
from field_lines import FieldLines
from viewer_wire_field_3D import Show3DVectorField

import numpy as np

# =============================================================================
# Define the geometry of the loops
# =============================================================================

# myLoops = ConcentricLoops()    
Rsphere=0.06
R_loop = 0.0125
frac_overlap = 0.3 # Rough estimate of the radius overlap between the layers
angle_layer = 2*np.arctan(R_loop*(1-frac_overlap)/Rsphere)

myLoops = Flower5loop(Rsphere=Rsphere, 
                   R_loop = R_loop,
                   angle_layer = angle_layer,
                   pos_center_sphere = [0, 0, 0], 
                   n_vec = [0, 0, 1],  
                   N_dscrt_pts_loop=100)

my_loops = myLoops.get_loops()

# =============================================================================
# Define the delays and the amplitude and the frequencies
# One day it will be setting the exciation function, B1(t)
# Maybe an extra feature soon (require a couple of calculation)
# =============================================================================

# Set various delays
my_B0 = [0, 1, 0]
my_freq  = 8000 # The exact value doesn't really matter for the current test
my_delays = myLoops.get_delays_optimal_rotation(my_B0, my_freq)

# Set various amplitude
my_amplitudes = np.zeros(len(my_delays)) + 1 # Same amplitude everywhere

# Set the frequencies
my_amplitudes = np.zeros(len(my_delays)) + 1 # Same amplitude everywhere

freq = 297.2*1e6 # Hz
my_freqs = np.zeros(len(my_delays)) + freq # Same frequency everywhere

# =============================================================================
# Calculate the field !
# =============================================================================


myB1 = FieldLines(my_loops,
                  my_delays,
                  my_amplitudes,
                  my_freqs)
# We will probe along a line in space, fixed time. 
t_probe = 0.50/freq # Fix the time to probe

# A line cute into the cage
# Staight cut along the z axis
list_probe_x = np.linspace(-1.5*Rsphere, 1.5*Rsphere, 87)
list_probe_y = 0.1*Rsphere      + 0*list_probe_x
list_probe_z = 0*Rsphere + 0*list_probe_x

print('Probing the field into space...')
(list_Bx, 
 list_By, 
 list_Bz) = myB1.get_field(list_probe_x, 
                           list_probe_y, 
                           list_probe_z,
                           t_probe) 

print('Probing the field into space... Done !')

# =============================================================================
# View the masterpiece  (Loops and RF field !)  
# =============================================================================
vv = Show3DVectorField()
# Show all the elements of the antennas
dd = np.array(my_delays)*my_freq
dd = dd - min(dd)
vv.add_list_wire(my_loops, list_color=dd, 
                      str_cmap='brg', label='My loop',
                      color_title='Time delay*freq (unitless))',
                      show_corner=False)
# Show the field
vv.add_vector_field(list_probe_x, list_probe_y, list_probe_z,
                    list_Bx, list_By, list_Bz,
                    scale=0.02,
                    label='B1 at t*f = %.f'%(t_probe*freq), 
                    color='k')
# Show the line probed
vv.add_line(list_probe_x, list_probe_y, list_probe_z, 
             label='Line probed', color='orange', show_corner=False)










