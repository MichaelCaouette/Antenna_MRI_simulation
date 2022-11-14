# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 16:40:20 2022

Goal: 
    Implement a birdcage geomtery with the delay offset

@author: mcaoue2
"""


import numpy as np

class BirdCage():
    """
    
    Definition of the geometry for a set of antennas.
    
    This is the famous birdcage. Cylindrial. 
    
    """
    def __init__(self, 
                 R=1, 
                 height=2, 
                 N_rung=8,
                 N_smoothloop=100):
        """
        TODO: next update should include to set the center position and 
        orientiation.
        
        Input in SI unit !
        
        R:
            Meter
            Radius of the Birdcage
            
        height:
            Meter
            Height of the birdcage (along the cylindrical axis)
        
        N_rung:
            (int)
            Number of rung of the birdcage. 
            
            
        N_smoothloop:
            (int)
            Number of element to have a smooth end ring (number of segment used
            to define it). The more the better (because smoother circle).
        """
        
        self.R = R
        self.height = height
        self.N_rung = N_rung
        
        self.N_smoothloop = N_smoothloop
        
        # Define the loop
        # It will create the object "self.list_loop"
        # That is returned by the method "get_wires"     
        self._define_wires()
           

        
    def _define_wires(self):
        """
        Define each element of wires (list of list of points)
        
        In this case, I have a loop at the center and 5 loops that overlap with
        it. But the loops are not making a whole circle. 
        
        """        
        # =============================================================================
        #         # Prepare the geometry and delay. 
        # =============================================================================
        self.list_wire  = [] # Will hold the various wire geometry
        # Each element consist of a line path and has an associated delay.
        
        # Set the rung
        for i in range(self.N_rung):
            angle = (i/self.N_rung)*2*np.pi
            # 3D pts
            p1 = (self.R*np.cos(angle), 
                  self.R*np.sin(angle),
                  -self.height/2)
            p2 = (self.R*np.cos(angle), 
                  self.R*np.sin(angle),
                  +self.height/2)
            list_vec = [p1, p2]
            self.list_wire.append(list_vec)
            
        # Set the end loop
        list_pts_top_ring = [] # Will store the points defining the top end ring
        list_pts_bot_ring = [] # Will store the points defining the bottom end ring
        for i in range(self.N_smoothloop+1): # Include the last point !
            s = i/self.N_smoothloop
            x = self.R*np.cos(s*2*np.pi)
            y = self.R*np.sin(s*2*np.pi)
            z_top = +self.height/2
            z_bot = -self.height/2
            list_pts_top_ring.append( (x, y, z_top) )     
            list_pts_bot_ring.append( (x, y, z_bot) ) 
        # Add these wires     
        self.list_wire.append(list_pts_top_ring)     
        self.list_wire.append(list_pts_bot_ring)

        
    def get_current_optimal_rotation(self, I_rung):
        """
        Get the list of current in each rung that optimizes the rotating RF 
        field in the center of the birdcage. 
        
        To get the optimal current, I used chapter 11 of the book from 
        J. Thomas Vaughan; John R. Griffiths

        I_rung:
            (Amp)
            Maximum electrical current in a rung. 
            
            
        """
        self.I_rung = I_rung
        
        # Get other quantity from that
        self.I_ER = self.I_rung/(2*np.sin(np.pi/self.N_rung))
        
        self.list_current = []
        # The current in the rungs
        for i in range(self.N_rung):
            self.list_current.append(self.I_rung)
            
        # The current in each of the two end rings
        self.list_current.append(self.I_ER)
        # VERY IMPORTANT the current is minus the top one !
        self.list_current.append(-self.I_ER)
        
        return self.list_current
    
    def get_delays_optimal_rotation(self,  freq):
        """
        Get the list of delays (ie phase) of each line, such that they have 
        the RF field an the maximum rotating component along the bircage axis.

        Parameters
        ----------
        
        freq: float
            Hz, Frequency of the RF field. 

        Returns
        -------
        None.

        """
        self.list_delay = [] # Will hold the corresponding delay
        
        self.w = 2*np.pi*freq
        
        # The delays in each rung
        for i in range(self.N_rung):
            angle = (i/self.N_rung)*2*np.pi
            # Add the delay
            self.list_delay.append(angle/self.w)
            
        # Add the delay in the two end ring. 
        # (no delay in the end rings, 
        # see J. Thomas Vaughan; John R. Griffiths)
        self.list_delay.append(0)        
        self.list_delay.append(0)              
        
        return self.list_delay
        

    def get_wires(self):
        """
        Return the list of loops that define the geometry of the system. 

        It is basically a list of list of points. 

        """
        return self.list_wire      
            
            
        
        
        
        
if __name__ == '__main__':
    # Check that the geometry is right        
        
    # from viewer_3D import plot_3D_cool
    from viewer_wire_field_3D import Show3DVectorField    
    
    self = BirdCage(R=0.3, 
                    height=0.5, 
                    N_rung=16,
                    N_smoothloop=50)
    
    my_loops = self.get_wires()
    my_freq  = 8000 # The exact value doesn't really matter for the this test
    my_delays = self.get_delays_optimal_rotation(my_freq)
    
    # =============================================================================
    # View the masterpiece    
    # =============================================================================
    vv = Show3DVectorField()
    # Show all the elements of the antennas
    dd = np.array(my_delays)*my_freq
    dd = dd - min(dd)
    vv.add_list_wire(my_loops, list_color=dd, 
                          str_cmap='brg', label='My loop',
                          color_title='Time delay*freq (unitless))',
                          show_corner=False)
    
    
# # =============================================================================
# #     # Check if the field rotates at the center
# # =============================================================================
#     list_t_probe = np.linspace(0, 1, 27)/self.f 
#     (list_Bx, 
#      list_By, 
#      list_Bz) = self.get_field(0.2*self.R, 
#                                0, 
#                                0, 
#                                list_t_probe)     
#     list_probe_axis = list_t_probe
#     plt.figure(tight_layout=True)
#     plt.subplot(311)
#     # plt.title('x, y = %.1f, %.1f radius'%(x/R, y/R))
#     plt.ylabel('Bx (SI)')    
#     plt.plot(list_probe_axis, list_Bx, label='',  linewidth=4)
#     plt.subplot(312)
#     plt.ylabel('By (SI)')    
#     plt.plot(list_probe_axis, list_By, label='',  linewidth=4)
#     plt.subplot(313)
#     plt.ylabel('Bz (SI)')    
#     plt.plot(list_probe_axis, list_Bz, label='',  linewidth=4)
#     plt.xlabel('Time ')
#     plt.legend()      
        
        
        
# # =============================================================================
# #     # Probe the field at various position and fixed time
# # =============================================================================
#     t_probe = 0.50/self.f # Fix the time to probe
    
#     # A line cute into the cage
#     # Staight cut along the z axis
#     list_probe_x = np.linspace(-1.5*self.R, 1.5*self.R, 87)
#     list_probe_y = 0.1*self.R      + 0*list_probe_x
#     list_probe_z = 0*self.height + 0*list_probe_x
    
#     list_Bz = []
#     list_Bx = []
#     list_By = []
#     print('Probing the fiekld into space...')
#     (list_Bx, 
#      list_By, 
#      list_Bz) = self.get_field(list_probe_x, 
#                                list_probe_y, 
#                                list_probe_z,
#                                t_probe) 
#     print('Done')
    
    

    
    
    
    
    
        
        
        
        