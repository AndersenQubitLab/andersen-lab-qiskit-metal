# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# This class was created by Figen YILMAZ
"""Grounded Tmon"""

from operator import length_hint
import numpy as np
from qiskit_metal import draw, Dict
from math import *
from qiskit_metal.draw.basic import buffer
from qiskit_metal.qlibrary.core import BaseQubit


class Tmon(BaseQubit):
    """The base `Tmon` class.

    Inherits `BaseQubit` class.

    Description:
        Create a grounded tunable tmon qubit with two josephson junction
        (see drawing below).

    Has only one connector line and it is added by using the `flux-bias line` dictionaries.
    The connector pad has a name and a list of default properties.

    Sketch:
        Below is a sketch of the qubit
        ::
                             0
                            
                                 |  Flux bias line
    fake_flux_bias_line  |      | |
                         _      | |
                        | |     | |
                        | |_____| |  
                        |_________| 
                         _________
                 +1      |  |  |  |       +1
                         |  x  x  |
            -1           |  |__|  |              +1    ^ Y
                         |  |  |  |                    |           
                         |  |  |  |                    |           
                         |  |  |  |                    |----->  X
                         |  |  |  |
                _________|  |  |  |________
               /   _________|  |_________   \           
              /  /                       \   \
             |  |                         |  |
              \  \_______________________/  /
           -1  \___________________________/     +1
             
                 -1                       -1


    .. image::
        Tmon.png
            
    .. meta::
        Tmon

    Default Options:
        * inductor_width: '10um' -- Width of the pseudo junction on the x-axis. Really just for simulating in HFSS / other EM software
        * jj_gap: '30um' -- Width of the pseudo junction on the y-axis. Really just for simulating in HFSS / other EM software
        * pad_head_width: '400um' -- The 'head' capacitance length along the x-axis
        * pad_head_length: '800um' -- The 'head' capacitance length along the y-axis
        * pad_arm_width: '100um' -- The 'equator' capacitance length along the y-axis
        * pad_arm_length: 1000um' -- The 'equator' capacitance length along the x-axis
        * palm_radius: '50um' -- Radius of the circle at the end of the pads (equator)
        * pad_gap: '80um' -- The distance between capacitance to ground plane, AKA etch gap
        * pad_armgap: '10um' -- The distance around the arms to the ground, arms will be in FX pocket.
        * orientation: '0' -- Degree of qubit rotation
        * flux_bias_line_options=Dict
            * make_fbl = True -- Boolean to make the flux bias line 
            * fbl_sep: '100um' -- The separation between the flux bias line and the jj along the y-axis
            * fbl_height: '50um' -- The height of the flux bias line along the y-axis
            * fbl_length: '30um' --  The width of the flux bias line along the x-axis
            * cpw_width: 'cpw_width' -- The width of the flux bias line
            * cpw_gap: 'cpw_gap' -- The dielectric gap width of the flux bias line
    """

    component_metadata = Dict(short_name='Tmon',
                              _qgeometry_table_path='True',
                              _qgeometry_table_poly='True',
                              _qgeometry_table_junction='True',
                             )
    """Component metadata"""

    # Default drawing options
    default_options = Dict(
        chip='main',
        inductor_width='10um',
        jj_gap='30um',
        pad_head_width='15um',
        pad_head_length='400um',
        pad_arm_width='15um',
        pad_arm_length='1400um',
        palm_radius='25um',
        pad_gap='100um',
        pad_armgap='5um',
        # 90 has dipole aligned along the +X axis,
        # while 0 has dipole aligned along the +Y axis
        orientation='0',
        flux_bias_line_options=Dict(
            make_fbl = False,
            fbl_sep='60um',
            fbl_height ='80um',
            fbl_length = '30um',
            cpw_width ='cpw_width',
            cpw_gap = 'cpw_gap',
        ))
    """Default drawing options"""

    TOOLTIP = """The base `Base Qubit` class."""

    def make(self):
        """Define the way the options are turned into QGeometry.

        The make function implements the logic that creates the geometry
        (poly, path, etc.) from the qcomponent.options dictionary of
        parameters, and the adds them to the design, using
        qcomponent.add_qgeometry(...), adding in extra needed
        information, such as layer, subtract, etc.
        """
        self.make_pocket()
       
        if self.p.flux_bias_line_options.make_fbl == True:
            self.make_flux_bias_line()
        
    def make_pocket(self):
        """Makes a Tmon shape grounded qubit with 3 arm, half-cross."""

        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p

        # since we will reuse these options, parse them once and define them as variables
        jj_gap=p.jj_gap
        pad_head_width = p.pad_head_width
        pad_head_length = p.pad_head_length
        pad_arm_width = p.pad_arm_width
        pad_arm_length = p.pad_arm_length
        palm_radius = p.palm_radius
        pad_gap = p.pad_gap
        pad_armgap = p.pad_armgap
        
        # Draw 'the arms', the capacitance. 
        pad_north = draw.rectangle(pad_head_width, pad_head_length, 0, pad_head_length/2) # pad_head drawned and define as pad_north
        pad_equator=draw.rectangle(pad_arm_length, pad_arm_width, 0, 0) # pad_arms will lean on the x-axis and defined as pad_equator
        # since arms defined before then bulp shape on the corners will be palm.
        pad_palm_left = draw.Point(-pad_arm_length/2, 0).buffer(palm_radius/2) 
        pad_palm_right = draw.Point(pad_arm_length/2, 0).buffer(palm_radius/2)
        # we union all the capacitance as 't' shape
        tmon = draw.union(pad_north, pad_equator, pad_palm_left, pad_palm_right)

        # The gap around it
        pad_gap_north = draw.rectangle(pad_head_width+2*pad_gap, pad_head_length+2*jj_gap, 0, pad_head_length/2)
        pad_gap_equator = draw.rectangle(pad_arm_length+2*pad_armgap, pad_arm_width+pad_gap*2, 0, 0)
        pad_palm_gap_left = draw.Point(-pad_arm_length/2, 0).buffer(palm_radius*1.5)
        pad_palm_gap_right = draw.Point(pad_arm_length/2, 0).buffer(palm_radius*1.5)
        # again we union all the gap part as one and called pad_etch. Will be etched away during fab
        pad_etch = draw.union(pad_gap_north, pad_gap_equator, pad_palm_gap_left, pad_palm_gap_right)

        # Draw the junction
        rect_jj = draw.LineString([(0, pad_head_length+jj_gap), 
                                    (0, pad_head_length)])
        # the draw.rectangle representing the josephson junction
        # rect_jj = draw.rectangle(p.inductor_width, pad_gap)

        # Rotate and translate all qgeometry as needed.
        polys = [rect_jj, tmon, pad_etch]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [rect_jj, tmon, pad_etch] = polys

        # Use the geometry to create Metal qgeometry
        self.add_qgeometry('poly', dict(tmon=tmon))
        self.add_qgeometry('poly',
                           dict(pad_etch=pad_etch),
                           subtract=True)
        # self.add_qgeometry('poly', dict(
        #     rect_jj=rect_jj), helper=True)
        self.add_qgeometry('junction',
                           dict(rect_jj=rect_jj),
                           width=p.inductor_width,
                           )

    def make_flux_bias_line(self):
        """ Adds flux bias line to the tmon."""
        # self.p allows us to directly access parsed values (string -> numbers) form the user option
 
        p = self.p
        pfb = self.p.flux_bias_line_options # parser on connector options

        # define commonly used variables once
        fbl_sep = pfb.fbl_sep
        fbl_length = pfb.fbl_length
        fbl_height = pfb.fbl_height
        cpw_width = pfb.cpw_width
        cpw_gap = pfb.cpw_gap

        # Define the geometry
        # Flux Bias Line
        
        # Draw the flux-bias line
        d = p.pad_head_length+p.jj_gap
        flux_bias_line = draw.Polygon([
             (-fbl_length/2-cpw_width, d+fbl_height+fbl_sep),   # point a
             (-fbl_length/2, d+fbl_height+fbl_sep),    # point b
             (-fbl_length/2, d+fbl_sep+cpw_width),   # point c
             (fbl_length/2, d+fbl_sep+cpw_width),   # point d
             (fbl_length/2, d+2*(fbl_height+fbl_sep)),   # point e
             (fbl_length/2+cpw_width, d+2*(fbl_height+fbl_sep)),   # point f
             (fbl_length/2+cpw_width, d+fbl_sep),   # point g
             (-fbl_length/2-cpw_width, d+fbl_sep),   # point i
        ])
         
        # Flux Bias line's gap part
        flux_bias_line_gap = draw.Polygon([
             (-fbl_length/2-cpw_width-cpw_gap, d+fbl_height+fbl_sep),   # point a
             (-fbl_length/2+cpw_gap, d+fbl_height+fbl_sep),    # point b
             (-fbl_length/2+cpw_gap, d+fbl_sep+cpw_width+cpw_gap),   # point c
             (fbl_length/2-cpw_gap, d+fbl_sep+cpw_width+cpw_gap),   # point d
             (fbl_length/2-cpw_gap, d+2*(fbl_height+fbl_sep)),   # point e
             (fbl_length/2+cpw_width+cpw_gap, d+2*(fbl_height+fbl_sep)),   # point f
             (fbl_length/2+cpw_width+cpw_gap, d+fbl_sep-cpw_gap),   # point g
             (-fbl_length/2-cpw_width-cpw_gap, d+fbl_sep-cpw_gap),   # point i
        ])
        
        # Flux-Bias Line CPW wire
        port_line = draw.LineString([(-(fbl_length/15-cpw_width), d+2*(fbl_height+fbl_sep)), 
                                    (fbl_length+cpw_gap, d+2*(fbl_height+fbl_sep))]) #point e

        # Fake Flux-Bias Line CPW wire
        fake_port_line = draw.LineString([(-(fbl_length+cpw_width), d+(fbl_height+fbl_sep)), 
                                    (0, d+(fbl_height+fbl_sep))]) #point e, kinda

        objects = [flux_bias_line, flux_bias_line_gap, port_line]
        objects = draw.rotate(objects, p.orientation, origin=(0, 0))
        objects = draw.translate(objects, p.pos_x, p.pos_y)
        [flux_bias_line, flux_bias_line_gap, port_line] = objects

        self.add_qgeometry('poly', {'flux_bias_line': flux_bias_line})
        self.add_qgeometry('poly', {'flux_bias_line_gap': flux_bias_line_gap}, subtract=True)        

        ####################################################################

        # add pins
        port_line_cords = list(draw.shapely.geometry.shape(port_line).coords)
        self.add_pin('flux_bias_line', 
                    port_line_cords, cpw_width)

        fake_port_line_cords = list(draw.shapely.geometry.shape(fake_port_line).coords)
        self.add_pin('fake_flux_bias_line', 
                    fake_port_line_cords, cpw_width)
