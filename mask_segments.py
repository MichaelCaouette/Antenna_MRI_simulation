# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 14:39:41 2022

Goal:
    Define the function for making masks 

@author: mcaoue2
"""

import numpy as np

def vec_magnitude(vec):
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

def vec_dot(v, u):
    # Do it element wise, because this way it is easier if the input 
    # contains arrays.
    return v[0]*u[0] + v[1]*u[1] + v[2]*u[2]

def get_segment_pts(s, p1, p2):
    """
    Get the x, y, z coordinate on the segment parametrized as 
    r_seg(s) = p1 + s*(p2-p1)

    Parameters
    ----------
    s : TYPE
        DESCRIPTION.
    p1 : TYPE
        DESCRIPTION.
    p2 : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # Use numpy array for the computation on vector
    q1 = np.array(p1)
    q2 = np.array(p2)
    dp = q2 - q1
    
    x = s*dp[0] + q1[0]
    y = s*dp[1] + q1[1]
    z = s*dp[2] + q1[2]
    
    return x, y, z    
    
    
def get_s(p1, p2, r):
    """
    Get the parametrization parameter, s, that minimize the distance between 
    r and the segment joining p1 and p2. 
    If the minimum distance is outside the segment, s will either be 0 or 1, if
    the nearest point is p1 or p2, respectively.
    
    The segment is parametrized as r_seg(s) = p1 + s*(p2-p1)

    See my noteboke 2022-11-10 in the section 
    "Smallest distance between pts and line"
    
    More information soon

    Parameters
    ----------
    p1 : TYPE
        DESCRIPTION.
    p2 : TYPE
        DESCRIPTION.
    r : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # Use numpy array for the computation on vector
    q1 = np.array(p1)
    q2 = np.array(p2)
    x, y, z  = r
    
    # Find s that minimize the distance between R and r_seg
    # See notebook 2022-XX-XX for the derivation
    dq = q2-q1
    
    dq2 = vec_dot(dq, dq)
    
    # Do element wise in case the input are array
    bb = [x - q1[0],
          y - q1[1],
          z - q1[2]]
    
    s_optimal= vec_dot(dq,bb) / dq2
    
    # Clip the minimum to be at the end point if it exceed
    s_optimal = np.where(s_optimal<0, 0, s_optimal) 
    s_optimal = np.where(s_optimal>1, 1, s_optimal) 
        
    return s_optimal


def get_min_dist(p1, p2, r, get_closest_pts=False):
    """
    The segment is parametrized as r_seg(s) = p1 + s*(p2-p1)

    See my noteboke 2022-11-10 in the section 
    "Smallest distance between pts and line"    
    
    More information soon. 

    Parameters
    ----------
    p1 : TYPE
        DESCRIPTION.
    p2 : TYPE
        DESCRIPTION.
    r : TYPE
        DESCRIPTION.
        
    get_closest_pts: bool
        If True, the function returns also the point on the segment that has 
        the smallest distance. The default is False. 

    Returns
    -------
    None.

    """
    # Get the optimal parameter
    s0 = get_s(p1, p2, r)
    
    # Get the point that has the minimum distance
    r_min= get_segment_pts(s0, p1, p2)

    # Use numpy array for the computation on vector
    x, y, z  = r 
    # And the corresponding distance with R
    # Do element wise in case the input are array
    d_min = vec_magnitude([r_min[0]-x,
                           r_min[1]-y,
                           r_min[2]-z])
    
    if get_closest_pts:
        return r_min, d_min
    else:
        return d_min

if __name__ == '__main__':
    # Verify that the formula wrk
    
    # Test various cases
    my_p1 = np.array([0, 4, 7] )
    my_p2 = np.array([8, 7, -3])
    my_r  = np.array([9, -4, 0])
    # my_r  = np.array([9, -4, 20])
    
    my_rmin, my_dmin = get_min_dist(my_p1, my_p2, my_r,
                           get_closest_pts=True)
    s_opt = get_s(my_p1, my_p2, my_r) 
    
    
    print()
    print('s_opt = ', s_opt)
    print('my_dmin = ', my_dmin)
    # Test the dot product, to see the perpendicularity
    dot_prod = vec_dot(my_rmin- my_r, my_p1-my_p2)
    print('Dot product = ', dot_prod)
    print('Use 3D visualisation to check that !')
    
    from mpl_toolkits import mplot3d
    import matplotlib.pyplot as plt

    def set_axes_equal(ax):
        '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
        cubes as cubes, etc..  This is one possible solution to Matplotlib's
        ax.set_aspect('equal') and ax.axis('equal') not working for 3D.
        
        This function it taken from the internet
    
        Input
          ax: a matplotlib axis, e.g., as output from plt.gca().
        '''
    
        x_limits = ax.get_xlim3d()
        y_limits = ax.get_ylim3d()
        z_limits = ax.get_zlim3d()
    
        x_range = abs(x_limits[1] - x_limits[0])
        x_middle = np.mean(x_limits)
        y_range = abs(y_limits[1] - y_limits[0])
        y_middle = np.mean(y_limits)
        z_range = abs(z_limits[1] - z_limits[0])
        z_middle = np.mean(z_limits)
    
        # The plot bounding box is a sphere in the sense of the infinity
        # norm, hence I call half the max range the plot radius.
        plot_radius = 0.5*max([x_range, y_range, z_range])
    
        ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    
    # Create the segment
    s_parms = np.linspace(0, 1, 100)
    xline, yline, zline = get_segment_pts(s_parms, my_p1, my_p2)
    

    ax.plot3D(xline, yline, zline, 'gray', label='segment')
    ax.plot3D(my_r[0], my_r[1], my_r[2], 'or', label='Point probed')
    ax.plot3D(my_rmin[0], my_rmin[1], my_rmin[2], 'ob', label='Smallest dist pts')
    # The line on the smallest distance
    xd, yd, zd = get_segment_pts(s_parms, my_rmin, my_r)
    ax.plot3D(xd, yd, zd, '--b', label='Smallest dist = %.1f\nDot produc =%.1f'%(my_dmin, dot_prod))
    
    ax.legend()
    set_axes_equal(ax)
    
    
    
    
    
    
    
