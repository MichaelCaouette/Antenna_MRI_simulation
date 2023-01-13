# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 14:31:32 2022

Goal:
    Impalement the concentric array

@author: mcaoue2
"""

import numpy as np

#Debug
_debug_enabled           = False
def _debug(*a):
    if _debug_enabled: 
        s = []
        for x in a: s.append(str(x))
        print(', '.join(s))

def get_loop_pts(x0, y0, z0, n_vec, R, Npts=50):
    """
    

    Parameters
    ----------
    R : Float
        Radius of the loop
    Npts : int, optional
        Number of point to define the loop. The default is 50.

    Returns
    -------
    None.

    """
    
    # Determine the roating axis and angle
    n_hat = np.array(n_vec)/np.linalg.norm(n_vec)
    
    # # Some way to get the angle from complex number. 
    # cx = n_hat[2]
    # cy = (1 - cx**2)**0.5
    # angle_rot = np.angle(cx+cy*1j) #  # Angle 
    angle_rot = np.arccos(n_hat[2]) # Angle of rotation
    
    
    axis_rot = np.cross(n_hat, [0, 0, 1]) # The reference is the z axis
    
    # Define the loop
    list_pts = [] # Will store the points defining the loop
    for i in range(Npts+1): # Include the last point !
        s = i/Npts
        
        # pts in x-y plane
        x1 = R*np.cos(s*2*np.pi)
        y1 = R*np.sin(s*2*np.pi)
        z1 = 0
        
        # Rotate it only if the axis of rotation is not null (in case the n_hat
        # was the reference)
        if not (np.linalg.norm(axis_rot) ==0):
            x2, y2, z2 = rotate_pts([x1, y1, z1], 
                                    axis_rot, 
                                    angle_rot)
        else:
            x2, y2, z2 = x1, y1, z1
        # Translate it
        x3, y3, z3 = (x2+x0, y2+y0, z2+z0)
        
        list_pts.append( (x3, y3, z3) )
        
    return list_pts

def rotate_pts(pts, n_vec, angle):
    """
    
    Apply a rotation from an axis and a angle. 
    Thanks to wikipedia:
    https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle

    Parameters
    ----------
    pts : list (x,y,z)
        3D point to rotate.
    n_vec : list(nx, ny, nz)
        Vectore pointing in the direction of the axis to rotate
    angle : float
        Angle to rotate

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    # From
    # https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    # But note that we could also use the Rodriguess rotation formalism, which
    # is intuitive (and equivalent I suppose)
    # https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula   
    # Get the normalized vector
    ux, uy, uz = np.array(n_vec)/np.linalg.norm(n_vec)
    
    # n_norm = (n_vec[0]**2 + 
    #           n_vec[1]**2 + 
    #           n_vec[2]**2   )**0.5
    # ux = n_vec[0] / n_norm
    # uy = n_vec[1] / n_norm
    # uz = n_vec[2] / n_norm
    
    cosA = np.cos(angle)
    sinA = np.sin(angle)
    
    # Build the rotation matrix
    R11 = cosA + ux**2*(1-cosA)
    R12 = ux*uy*(1-cosA) - uz*sinA
    R13 = ux*uz*(1-cosA) + uy*sinA
    
    R21 = uy*ux*(1-cosA) + uz*sinA
    R22 = cosA + uy**2*(1-cosA)
    R23 = uy*uz*(1-cosA) - ux*sinA
    
    R31 = uz*ux*(1-cosA) - uy*sinA
    R32 = uz*uy*(1-cosA) + ux*sinA
    R33 = cosA + uz**2*(1-cosA)
    
    rot_matrix = np.matrix([[R11, R12, R13],
                            [R21, R22, R23],
                            [R31, R32, R33]])
    
    my_dot = np.dot(pts, rot_matrix)
    
    return np.array(my_dot)[0]    
    
class ConcentricLoops():
    """
    
    Definition of the geometry for a set of antennas. 
    
    The geometry is that we define the loops as circles lying on a sphere. 
    The first loop is at the top of the sphere (say the z axis). 
    The next set of loops are concentric around it, with a lower angle and 
    uniformly distributed around the z axis, still lying on a sphere. 
    There could be another set of loops at another angle. 
    
    """
    def __init__(self, 
                 Rsphere=1, 
                 list_Rcircle = [0.5, 0.45],
                 list_Nloop = [20],
                 list_angle_layer = [0.4*0.5*np.pi],
                 pos_center_sphere = [0, 0, 0], 
                 n_vec = [0, 0, 1], 
                 N_dscrt_pts_loop=100):
        """
        Description coming soon
        
        list_Rcircle: list of float
            meter
            List the radius of the loops for each layer. 
            Including the first loop at the top. 
            
        list_Nloop: list of int
            List the number of loop in each layer, starting by the second. 
            Therefore we exclude the top loop, because it is always 1.
            Therefore the lenght of the list is this lenght of
            list_Rcricle minus one. 
            
        list_angle_layer: list of float
            Radian
            List of angle between each layer. 
            The first angle of the list is the angle between the center of the 
            top loop and the center of all the other loops in the layer. 
            The second angle of the list is between the center of the loop of 
            the first layer and the center of the loop of the next layer. 
            Etc.
            Therefore the lenght of this list is the lenght of
            list_Rcricle minus one.            
        """
        _debug('ConcentricLoops.__init__')
        
        self.Rsphere = Rsphere
        self.list_Rcircle = list_Rcircle
        self.list_Nloop = list_Nloop
        self.list_angle_layer = list_angle_layer
        self.pos_center_sphere = pos_center_sphere
        self.n_vec = n_vec
        self.N_dscrt_pts_loop = N_dscrt_pts_loop

        # Number of layer, excluding the top loop. 
        self.N_layer = len(self.list_Nloop) 
        
        # Define my own basis, where the z points in the n_vec direction
        (self.z_hat, 
         self.x_hat, 
         self.y_hat) = self.get_basis(self.n_vec)
        
        # =============================================================================
        #         # Prepare the geometry and delay. 
        # =============================================================================
        self.list_wire  = [] # Will hold the various wire geometry
        # Each element consist of a line path and has an associated delay.
        
        self._define_concentric_loops()
        
    def cross_prod(self, v, u):
        """
        Compute the cross product between two vector.
        
        v and u:
            Both tuple (x, y z) defining a 3D vector.
        
        return:
            The cross product v x u
        """
        # Cross product between the lenght element along the loop and the 
        # relative position. 
        # Simpler to do it explicitely and element wise, because using 
        # np.cross() makes it hard to interprete multiple axis in x,y,z,t
        cross_x = v[1]*u[2] - v[2]*u[1]
        cross_y = v[2]*u[0] - v[0]*u[2]
        cross_z = v[0]*u[1] - v[1]*u[0]
        return  np.array([cross_x, cross_y, cross_z], dtype=np.ndarray)
    
    def magnitude(self, vec):
        """
        Compute the magnitude of the vector vec. 
        
        vec:
            Tuple (x, y, z) of the vector 
            
        return:
            The magnitude of vec. 
        """
        
        # Do it element wise, because this way it is easier if 'vec' contains
        # arrays.
        return (vec[0]*vec[0] + 
                vec[1]*vec[1] + 
                vec[2]*vec[2]   )**0.5

    def dot(self, v, u):
        # Do it element wise, because this way it is easier if the input 
        # contains arrays.
        return v[0]*u[0] + v[1]*u[1] + v[2]*u[2]
    
    def get_basis(self, a_axis):
        """
        Get 3 orthogonal and normalized vector, one which is the input axis. 
        The angle of the two other is a bit arbitrary.
        
        """
        a_axis = np.array(a_axis)
        
        a1, a2, a3 = a_axis
        
        # Find the second axis
        # The choose of the axis depends on which component of a_axis is 
        # non-null. 
        # At least one component of a_xis must be non zero
        if a1 != 0:
            b1 = -(a2+a3)/a1
            b2, b3 = b1*0+1, b1*0+1 #Not writting 1, in case the input as arrays. But anyway, I am already taking a condition on a1
        elif a2 != 0:
                b2 = -(a1+a3)/a2
                b1, b3 = b2*0+1, b2*0+1
        elif a3 != 0:
            b3 = -(a1+a2)/a3
            b1, b2 = b3*0+1, b3*0+1
        
        b_axis = np.array([b1, b2, b3], dtype=np.ndarray)
        
        # Get the third axis, which is perpendicular to the two first axis
        c_axis = self.cross_prod(a_axis, b_axis)
            
        # Normalize that
        a_hat = a_axis/self.magnitude(a_axis)
        b_hat = b_axis/self.magnitude(b_axis)
        c_hat = c_axis/self.magnitude(c_axis)
        
        return a_hat, b_hat, c_hat

    def get_wires(self):
        """
        Return the list of loops that define the geometry of the system. 

        It is basically a list of list of points. 

        """
        return self.list_wire   
    
    def _define_concentric_loops(self):
        """
        Set the wires to be the concentric loop (as defined in the description)
        of the class.
        
        """
        _debug('ConcentricLoops._define_bent_loops')
        
        
        # =============================================================================
        #         # First define the top loop
        # =============================================================================
        # This will list all the orientation of each wire. 
        # Usefull for the delay
        # Initiate the list with the top loop
        self.list_orientation = [self.n_vec]
        
        # The first loop lay on the sphere in the direction of the choosent vector
        x0 = self.pos_center_sphere[0] + self.Rsphere*self.n_vec[0]
        y0 = self.pos_center_sphere[1] + self.Rsphere*self.n_vec[1]
        z0 = self.pos_center_sphere[2] + self.Rsphere*self.n_vec[2]
        
        # Get the list of pts making the loop
        list_pts = get_loop_pts(x0, y0, z0, 
                                self.z_hat,
                                self.list_Rcircle[0], 
                                Npts=self.N_dscrt_pts_loop)        
        self.list_wire.append(list_pts )
        
        # =============================================================================
        # Add the other layers
        # =============================================================================
        
        # Orientation of the previous loop layer
        self.n_prime = self.z_hat 
        
        for i in range(self.N_layer):
            N_loop_layer = self.list_Nloop[i]
            
            # The layer is rotated with respect to the previus layer. 
            self.n_prime = rotate_pts(self.n_prime,
                                      self.x_hat, 
                                      self.list_angle_layer[i])
                
            for j in range(N_loop_layer):
                
                # Find the position of the center of the loop
                # We defined our own x,y,z axis with respect to a reference
                # So, within a layer, we rotate around the z axis.
                self.n_loop = rotate_pts(self.n_prime,
                                         self.z_hat , 
                                         j*2*np.pi/N_loop_layer)
                # Save this orientation
                self.list_orientation.append(self.n_loop)

                
                
                # The first loop lay on the sphere in the direction of the choosent vector
                x0 = self.pos_center_sphere[0] + self.Rsphere*self.n_loop[0]
                y0 = self.pos_center_sphere[1] + self.Rsphere*self.n_loop[1]
                z0 = self.pos_center_sphere[2] + self.Rsphere*self.n_loop[2]
               
                # Get the list of pts making the loop
                list_pts = get_loop_pts(x0, y0, z0, 
                                        self.n_loop,
                                        self.list_Rcircle[i+1], 
                                        Npts=self.N_dscrt_pts_loop)        
                self.list_wire.append(list_pts )  
                
        
        # Note the number of element that makes our field
        self.N_wire = len(self.list_wire)        

        return      
    
    def get_delays_optimal_rotation(self, n_B0, freq):
        """
        Get the list of delays (ie phase) of each line, such that they have 
        the RF field an the maximum rotating component aroung the B0 field. 

        Parameters
        ----------
        n_B0 : (nx, ny, nz)
            Define the direction of the B field. 
        
        freq: float
            Hz, Frequency of the RF field. 

        Returns
        -------
        None.

        """
        
        # Define a plane perpendicular to B0. 
        # n1 points in the B0 direrction
        # n2, n3 are perpend to B0 and orthogonal to each other
        n1, n2, n3 = self.get_basis(n_B0)
        # Reset the list of delay
        self.list_delay = []
        for n_loop in self.list_orientation:
            # Project this orientation on the plane perpendicular to B0
            # nx, ny are the orthogonal components in the plane perp to B0
            nx = self.dot(n_loop, n2)
            ny = self.dot(n_loop, n3)
            
            # Get the angle between this vector and a reference vector of 
            # the plane
            # Use the complex formalism to easily get the angle from 0 to 2*pi
            cplx_vec = nx + 1j*ny
            angle = np.angle(cplx_vec)
            
            # The delay should warp around a period (1/f) after an angle of 2*pi
            delay = angle/(2*np.pi*freq)
            
            self.list_delay.append(delay)  
            
        return self.list_delay
    
    
if __name__ == '__main__':
    _debug_enabled = True
       
    # from viewer_3D import plot_3D_cool
    from viewer_wire_field_3D import Show3DVectorField


    # self = ConcentricLoops()    
    Rsphere=0.06
    R_loops = 0.0125
    frac_overlap = 0.3 # Rough estimate of the radius overlap between the layers
    
    angle_layer = 2*np.arctan(R_loops*(1-frac_overlap)/Rsphere)
    
    # self = ConcentricLoops(
    #                  Rsphere=Rsphere, # m, radius of the sphere on which lays the coil
    #                  list_Rcircle = [R_loops, R_loops],# m, radiues of each loop
    #                  list_Nloop = [6], # # We should verify the optimal number
    #                  list_angle_layer = [angle_layer], # This angle should be check
    #                  pos_center_sphere = [0, 0, 0], # Origin of the the sphere
    #                  n_vec = [0, 0, -1],  # Top loop direction
    #                  I_loop=1,   # The current should not matter for the overlap shape
    #                  N_dscrt_pts_loop=85,
    #                  f = 297.2*1e6 )

    self = ConcentricLoops(
                     Rsphere=Rsphere, # m, radius of the sphere on which lays the coil
                     list_Rcircle = [R_loops, R_loops, R_loops],# m, radiues of each loop
                     list_Nloop = [6, 8], # # We should verify the optimal number
                     list_angle_layer = [angle_layer, angle_layer], # This angle should be check
                     pos_center_sphere = [0, 0, 0], # Origin of the the sphere
                     n_vec = [0, 0, -1],  # Top loop direction
                     N_dscrt_pts_loop=85)
    
    # Set various delays
    my_B0 = [0, 1, 0]
    my_freq  = 8000 # The exact value doesn't really matter for the current test
    my_delays = self.get_delays_optimal_rotation(my_B0, my_freq)
    my_loops = self.get_wires()
    
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
    vv.add_vector_field([0], [0], [-0.04],
                        [my_B0[0]], [my_B0[1]], [my_B0[2]],
                        scale=0.02,
                        label='B0 direction', color='k')

            
    # =============================================================================
    # Test some geomtry        
    # =============================================================================
    self = ConcentricLoops(
                     Rsphere=0.04, # m, radius of the sphere on which lays the coil
                     list_Rcircle = [0.020, 0.020],# m, radiues of each loop
                     list_Nloop = [6], # # We should verify the optimal number
                     list_angle_layer = [40*np.pi/180], # This angle should be check
                     pos_center_sphere = [0, 0, 0], # Origin of the the sphere
                     n_vec = [0, 0, -1],  # Top loop direction
                     N_dscrt_pts_loop=85)
    
    vv = Show3DVectorField()
    vv.add_list_wire(my_loops, list_color=dd, 
                          str_cmap='brg', label='My loop',
                          color_title='Time delay*freq (unitless))',
                          show_corner=False)      
        
    # Check the clearance at the middle
    # This is related to the distance of the layer loops with the top loop axis. 
    x_constraint = self.Rsphere*np.sin(self.list_angle_layer[0]) - self.list_Rcircle[0]*np.cos(self.list_angle_layer[0])
    clearance = 2*x_constraint
    print('')
    print('Check-up the clearance')
    print('Rsphere = %.2f mm'%(self.Rsphere*1e3))
    print('R_loops = %.2f mm'%(self.list_Rcircle[0]*1e3))
    print('angle_layer = %.1f degree'%(self.list_angle_layer[0]*180/np.pi))
    print('Clearane = %.2f mm '%(clearance*1e3))      
    # Check the height of the cap sphere
    h1 = self.Rsphere*np.cos(self.list_angle_layer[0]) - self.list_Rcircle[0]*np.sin(self.list_angle_layer[0])
    h2 = self.Rsphere-h1
    print('h1 = %.2f mm'%(h1*1e3))
    print('h2 = %.2f mm'%(h2*1e3))
        
        
        
        
        
        
        
        
        
    
    
    
    