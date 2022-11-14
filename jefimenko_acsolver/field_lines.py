# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 11:33:08 2022

Get the RF field from a set of lines

@author: mcaoue2
"""


from field_line import FieldLine

import numpy as np
import traceback
# Very usefull command to use for getting the last-not-printed error
# Just type _p() in the console 
_p = traceback.print_last

#Debug
_debug_enabled           = False
def _debug(*a):
    if _debug_enabled: 
        s = []
        for x in a: s.append(str(x))
        print(', '.join(s))
        
class FieldLines():
    """
    
    Get the RF field from a set of lines. 
    
    The user input the set of lines, and the class calculate the RF field. 
    
    """
    def __init__(self, 
                 list_lines,
                 list_delays,
                 list_I,
                 list_freq,
                 u0 = 4*np.pi*1e-7,
                 c = 299792458 ):
        """
        More Description coming soon. 
        The units must be SI !
        
        list_lines: list
            Each element in the list is a list of points defining the line. 
        list_delay: list
            Each element in the list is the delays of the RF signal in the 
            line.
        list_I: list
            Each element in the list is the amplitude of the electric current 
            in the line. Ampere. (SI unit!)
        
        list_freq: list of freq
            Each frequency f is described as follow. 
            f:
                (Hz)
                Frequency of the current oscillation in the lines. For now it is 
                the same current in each line. 
                Note that t is not w !  
                For example,  2*pi*f = w, such that gamma*B=w=2*pi*f 
                (For n magnetization dynamic)

        u0 = 4*np.pi*1e-7:
            H/m Magnetic constant
            
        c =  299792458:
             m/s Speed of light in vaccuum
        """
        _debug('FieldLines.__init__')
        
        #TODO verify that the lenght match
        
        self.list_lines  = list_lines
        self.list_delays = list_delays
        self.list_I      = list_I
        self.list_freq   = list_freq
        
        # A line is a wire
        self.N_wire = len(self.list_lines)
        
        # We are not super to use this variable other then sending them to 
        # the field line calculator
        self.u0 = u0
        self.c  = c 
        
        # Initiate the field calulcation
        self.list_field_line = []
        for i in range(self.N_wire):
            field_calc = FieldLine(list_vec=self.list_lines[i],
                                   I  =self.list_I[i], 
                                   f  =self.list_freq[i],
                                   u0 = self.u0,
                                   c  = self.c)
            self.list_field_line.append(field_calc)
            
    
    def get_field_from_spatial(self, vcos, vsin, t):
        """
        Convenient method to get the field when we already have the spatial 
        components. 
        """
        _debug('FieldLines.get_field_from_spatial')
        
        for i in range(self.N_wire):
            field_calc = self.list_field_line[i]
            out = field_calc.get_field_from_spatial(vcos, vsin, t)
            if i == 0:
                # Initiate the field
                self.Bx, self.By, self.Bz = out
            else:
                # Add the field
                self.Bx += out[0]
                self.By += out[1]
                self.Bz += out[2]            
        # Done !
        return self.Bx, self.By, self.Bz               
        
    def get_field_spatial(self, x, y, z, print_progress=False):
        """
        Get only the spatial contribution of the field (ie the two vector:
            the cos and sin term)
         Useful for studying the rotating components   
         
         
         return self.vcos, self.vsin
        """
        _debug('FieldLines.get_field_spatial')
        
        for i in range(self.N_wire):
            if print_progress:
                print('\nCalculating wire # %d / %d'%(i, self.N_wire))
            # Calculation of the field without the delay
            field_calc = self.list_field_line[i]  
            
            out =  field_calc.get_field_spatial(x, y, z, 
                                                print_progress=print_progress)
            # Add the effect of the phase delay
            w = self.list_freq[i]*2*np.pi
            phase = w*self.list_delays[i]
            A, D = out
            # The effect is a slight shift in the cos and sin component
            # (See my notebook 2022-07-13 for the trigonometry trick)
            vcos = A*np.cos(phase) - D*np.sin(phase)
            vsin = A*np.sin(phase) + D*np.cos(phase)
            
            # Add them
            if i == 0:
                # Initiate the answer
                self.vcos, self.vsin = vcos, vsin
            else:
                self.vcos += vcos
                self.vsin += vsin
        
        return self.vcos, self.vsin
        
    def get_field(self, x, y, z, t, print_progress=False):
        """
        Get the magnetic field a time t and position x,y,z. 
        
        """
        _debug('FieldLines.get_field')
        
        # We will add the field from each element plus the delays
        
        for i in range(self.N_wire):
            if print_progress:
                print('\nCalculating wire # %d / %d'%(i, self.N_wire))
            # Important to add the delay
            field_calc = self.list_field_line[i]  
            out =  field_calc.get_field(x, y, z, 
                                        t-self.list_delays[i],
                                        print_progress=print_progress)
            if i == 0:
                # Initiate the field
                self.Bx, self.By, self.Bz = out
            else:
                # Add the field
                self.Bx += out[0]
                self.By += out[1]
                self.Bz += out[2]
        # Done !
        return self.Bx, self.By, self.Bz    
    

if __name__ == '__main__':
    # I should test some field calculation (require to )    
    print('Please ')
    
    