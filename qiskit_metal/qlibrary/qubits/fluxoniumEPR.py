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

# This class was created by Roald van den Boogaart, Christian Kraglund Andersen, Figen YILMAZ

from re import sub
import numpy as np
from qiskit_metal import draw, Dict
from math import *
from qiskit_metal.draw.basic import buffer
from qiskit_metal.qlibrary.core import BaseQubit, QRoutePoint

from collections import OrderedDict
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander


class FluxoniumPocket(BaseQubit):
    """The base `FluxoniumPocket` class.

    Inherits `BaseQubit` class.

    Create a standard pocket fluxonium qubit for a ground plane,
    with two pads connected by a junction (see drawing below).

    Connector lines can be added using the `connection_pads`
    dictionary. Each connector pad has a name and a list of default
    properties.

    Sketch:
        Below is a sketch of the qubit
        ::
                               0
        
                             | | |  Charge line
                 +1          |___|           +1
                _______________________________
            -1 |              ___              | +1        Y
               |             /   \             |           ^   
               |             \   /             |           |
               |              | |___           |           |----->  X
               |              |_|   |    ______|
               |               |    |   |  ____|
               |               x    |   | |____|__
               |               |    |   |_________-- Flux bias line
               |              | |___|          |
               |              | |              |
               |             /   \             |
               |             \___/             |
            -1 |_______________________________|  +1
                 -1        ___________      -1
                          /  _______  \
                         /  /       \  \            
                         \  \_______/  /
                          \___  |  ___/
                              | | | 
                          Read out line

    .. image::
        FluxoniumPocket.png

    BaseQubit Default Options:
        * pos_x: '0um'
        * pos_y: '0um'
        * connection_pads: Empty Dict -- The dictionary which contains all active connection lines for the qubit.
        * _default_connection_pads: Empty Dict -- The default values for the (if any) connection lines of the qubit.

    Default Options:
        * pos_x: '0um' -- Where the center of the pocket should be located on chip, on x axis
        * pos_y: '0um' -- Where the center of the pocket should be located on chip, on y axis
        * pad_gap: '30um' -- The distance between the two charge islands, which is also the resulting 'length' of the pseudo junction
        * inductor_width: '20um' -- Width of the pseudo junction between the two charge islands (if in doubt, make the same as pad_gap). Really just for simulating in HFSS / other EM software
        * pad_width: '30um' -- The width (x-axis) of the charge island pads
        * pad_height: '300um' -- The size (y-axis) of the charge island pads
        * pad_radius: '80um' -- Radius of the circle at the end of the pads
        * l_length: '1um' -- Length of the kinetic inductor along the y-axis
        * l_arm_length: '200um' -- Length of the arm of the kinetic inductor along x-axis
        * l_inductance:
        * l_ind_per_square:
        * L_j:
        * l_fillet:
        * pocket_width: '800um' -- Size of the pocket (cut out in ground) along x-axis
        * pocket_height: '800um' -- Size of the pocket (cut out in ground) along y-axis
        * orientation: '0' -- Degree of qubit rotation
        * flux_bias_line_options=Dict
            * make_fbl = True -- Boolean to make the flux bias line 
            * fbl_sep='50um' -- The separation between the flux bias line and the inductor along the x-axis
            * fbl_height ='50um' -- The height of the flux bias line along the y-axis
            * cpw_width ='cpw_width' -- The width of the flux bias line
            * cpw_gap = 'cpw_gap' -- The dielectric gap width of the flux bias line    
        * charge_line_options=Dict
            * make_cl = True -- Boolean to make the charge line
            * cl_length
            * cl_sep ='15um' -- The separation between the connection pade and the pocket the y-axis
            * cpw_width='10um' -- The width of the charge line
            * cpw_gap= '10um' -- The dielectric gap width of the charge line
            * loc_W=
            * loc_H=
        * readout_line_options=Dict(
            * make_rol = True -- Boolean to make the readout line
            * pad_sep='50um' -- The separation between the connection pad and the capacitor pad the y-axis
            * pad_width='150um' -- Width of the connection pad along the x-axis  
            * pad_height='50um', -- Height of the connection pad along the y-axis
            * pad_gap='10um' -- Dielectric gap width of the connection pad
            * cpw_width='cpw_width', -- The width of the charge line
            * cpw_gap='cpw_gap' -- The dielectric gap width of the readout line
            * loc_W=
            * loc_H=
    """

    component_metadata = Dict(short_name='FluxoniumPocket',
                              _qgeometry_table_path='True',
                              _qgeometry_table_poly='True',
                              _qgeometry_table_junction='True',
                             # _qgeometry_table_inductor='True'
                             )
    """Component metadata"""

    # Default drawing options
    default_options = Dict(
        chip='main',
        pos_x='0um',
        pos_y='0um',
        pad_gap='30um',
        inductor_width='10um',
        pad_width='15um',
        pad_height='200um',
        pad_radius='50um',
        l_length='1um',
        l_arm_length='50um',
        l_inductance='200nH',
        l_ind_per_square='2nH',
        L_j = '16.35nH',
        l_fillet = '5um', 
        pocket_width='800um',
        pocket_height='800um',
        # 90 has dipole aligned along the +X axis,
        # while 0 has dipole aligned along the +Y axis
        orientation='0',
        flux_bias_line_options=Dict(
            make_fbl = False,
            fbl_sep='100um',
            fbl_height ='50um',
            cpw_width ='cpw_width',
            cpw_gap = 'cpw_gap',
        ),
        charge_line_options=Dict(
            make_cl = False,
            cl_length = '100um',
            cl_sep ='15um',
            cpw_width='cpw_width',
            cpw_gap= 'cpw_gap',
            loc_W='0',  # width location only 0,
            loc_H='+1',  # height location only +1 or -1
        ),
        readout_line_options=Dict(
            make_rol = False,
            pad_sep='100um',
            pad_width = '150um',
            pad_height = '50um',
            cpw_width='cpw_width',
            cpw_gap='cpw_gap',
            loc_W='0',  # width location only 0,
            loc_H='-1',  # height location only +1 or -1
        ))
    """Default drawing options"""

    TOOLTIP = """The base `FluxoniumPocket` class."""

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
        if self.p.charge_line_options.make_cl == True:
            self.make_charge_line()
        if self.p.readout_line_options.make_rol == True:
            self.make_readout_line()
        
    def make_pocket(self):
        """Makes standard fluxonium in a pocket."""
        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p

        # Draw the junction
        rect_jj = draw.LineString([(0, -p.pad_gap / 2), (0, +p.pad_gap / 2)])

        # Draw the pads (shapely polygons)
        pad_rect = draw.rectangle(p.pad_width, p.pad_height - p.pad_radius)
        pad_circle = draw.Point(0,(p.pad_height - p.pad_radius) / 2.).buffer(p.pad_radius)
        pad = draw.union(pad_rect, pad_circle)
        pad_top = draw.translate(pad, 0, +(p.pad_height + p.pad_gap - p.pad_radius) / 2.)
        pad_bot = draw.rotate(pad_top, 180, origin=(0, 0))

        # Calculating total length of the inductor
        ind_per_square = float(p.l_ind_per_square.replace('nH',''))
        l_inductance = float(p.l_inductance.replace('nH',''))
        L_tot = l_inductance*p.l_length/ind_per_square

        # Draw fillets for inductor meander
        outer_circle = draw.Point(0,0).buffer(p.l_fillet + p.l_length/2)
        inner_circle = draw.Point(0,0).buffer(p.l_fillet - p.l_length/2)
        circle = draw.subtract(outer_circle, inner_circle)
        semicircle = draw.subtract(circle, draw.rectangle(p.l_fillet + p.l_length / 2., 2. * p.l_fillet + p.l_length, - (p.l_fillet + p.l_length / 2.) / 2.,0))
        poly_seg = draw.subtract(semicircle, draw.rectangle(2. * p.l_fillet + p.l_length, p.l_fillet + p.l_length / 2., 0, - (p.l_fillet + p.l_length / 2.) / 2.))
        
        # Place and orientate copies of the poly_seg 
        poly_wire_top = []
        poly_wire_bot = []
        mount_top = draw.translate(draw.rectangle(p.l_arm_length, p.l_length), (p.pad_width + p.l_arm_length)/ 2., +L_tot / 2.+ p.l_fillet)
        mount_bot = draw.translate(draw.rectangle(p.l_arm_length, p.l_length), (p.pad_width + p.l_arm_length)/ 2., -L_tot / 2.- p.l_fillet)
        seg_top = draw.translate(poly_seg, p.l_arm_length + p.pad_width / 2. , L_tot / 2.)
        seg_bot = draw.translate(draw.rotate(poly_seg, -90, origin = (0,0)), p.l_arm_length + p.pad_width / 2. , -L_tot / 2.)

        poly_wire_top.append(mount_top)
        poly_wire_top.append(seg_top)
        poly_wire_top = draw.union(poly_wire_top)
        poly_wire_bot.append(seg_bot)
        poly_wire_bot.append(mount_bot)
        poly_wire_bot = draw.union(poly_wire_bot)

        # Create linestring inductor segments
        c = (p.pad_width / 2. + p.l_arm_length + p.l_fillet, +L_tot / 2. )
        d = (p.pad_width / 2. + p.l_arm_length + p.l_fillet, -L_tot / 2. )  
        inductor = draw.LineString([c,d])
        
        # Draw the pocket
        rect_pk = draw.rectangle(p.pocket_width, p.pocket_height)

        # Rotate and translate all qgeometry as needed.
        polys = [rect_jj, pad_top, pad_bot, rect_pk, poly_wire_top, poly_wire_bot, inductor]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [rect_jj, pad_top, pad_bot, rect_pk, poly_wire_top, poly_wire_bot, inductor] = polys

        # Use the geometry to create Metal qgeometry
        
        self.add_qgeometry('poly', dict(pad_top=pad_top, pad_bot=pad_bot, poly_wire_top = poly_wire_top, poly_wire_bot = poly_wire_bot))
        self.add_qgeometry('poly', dict(rect_pk=rect_pk), subtract=True)
        
        self.add_qgeometry('path', 
                           dict(inductor=inductor),
                           width = p.l_length,
                           hfss_inductance = str(l_inductance)+'nH')
        self.add_qgeometry('junction',
                           dict(rect_jj=rect_jj),
                           width=p.inductor_width,
                           hfss_inductance = p.L_j)

    def make_flux_bias_line(self):
        """ Adds flux bias line to fluxonium pocket."""
        # Grab option values
        p = self.p
        pfb = self.p.flux_bias_line_options

        fbl_sep = pfb.fbl_sep
        fbl_height = pfb.fbl_height
        cpw_width = pfb.cpw_width
        cpw_gap = pfb.cpw_gap

        fbl_length = p.pocket_width/2 - fbl_sep
        # Define the geometry
        # Flux Bias Line

        d = p.pocket_width/2 # The position of flux bias line on the x-axis, starting point inside the pocket

        flux_bias_lineup = draw.Polygon([
             (d, fbl_height/2),   # point a
             (d, fbl_height/2+cpw_width),    # point b
             (fbl_sep, fbl_height/2+cpw_width),   # point c
             (fbl_sep, fbl_height/2),   # point d
        ])
        flux_bias_linemid = draw.Polygon([
             (fbl_sep-cpw_width, fbl_height/2),   # point e
             (fbl_sep, fbl_height/2),    # point f
             (fbl_sep, -fbl_height/2),   # point g
             (fbl_sep-cpw_width, -fbl_height/2),   # point h
        ])

        # Here we make flux-bias line curvy for the top side and also bottom side and the union all of them
        circle_top = draw.shapely.geometry.Point(fbl_sep, fbl_height/2).buffer(cpw_width)
        cut_ply = draw.Polygon([
             (fbl_sep*2, fbl_height+cpw_width),   # point o
             (fbl_sep, fbl_height+cpw_width ),    # point p
             (fbl_sep, fbl_height/2),   # point r same with point d
             (fbl_sep/2, fbl_height/2),   # point s
             (fbl_sep/2, -fbl_height+cpw_width),  # point t
             (fbl_sep*2, -fbl_height+cpw_width),  # point u
        ])
        circle_top = draw.subtract(circle_top, cut_ply)
        # same goes for bottom edge
        circle_bot = draw.shapely.geometry.Point(fbl_sep, -fbl_height/2).buffer(cpw_width)
        cut_ply2 = draw.Polygon([
             (fbl_sep-cpw_width, -fbl_height/2),   # point v same with h or i
             (fbl_sep-cpw_width, fbl_height*2),    # point y
             (fbl_sep, fbl_height),   # point z 
             (fbl_sep, fbl_height),   # point w
             (fbl_sep*2, fbl_height),  # point x
             (fbl_sep*2, -fbl_height/2),  # point k
        ])
        circle_bot = draw.subtract(circle_bot, cut_ply2)
 
        flux_bias_linebot = draw.Polygon([
             (d+fbl_sep/2, -fbl_height/2),   # point i
             (d+fbl_sep/2, -fbl_height/2-cpw_width),    # point k
             (fbl_sep, -fbl_height/2-cpw_width),   # point l
             (fbl_sep, -fbl_height/2),   # point m
        ])

        flux_bias_line = draw.union(flux_bias_lineup, flux_bias_linemid,  flux_bias_linebot, circle_top, circle_bot)


        # Flux Bias line's gap part, outside the GND
        flux_bias_line_gap = draw.rectangle(fbl_sep/2, cpw_width+cpw_gap*2, d+fbl_sep/4 ,-fbl_height/2-cpw_width/2)

        # Flux-Bias Line CPW wire
        port_line = draw.LineString([((d+fbl_sep/2), 0), 
                                    ((d+fbl_sep/2), -(fbl_height+cpw_width))])

        objects = [flux_bias_line, flux_bias_line_gap, port_line]
        objects = draw.rotate(objects, p.orientation, origin=(0, 0))
        objects = draw.translate(objects, p.pos_x, p.pos_y)
        [flux_bias_line, flux_bias_line_gap, port_line] = objects


        self.add_qgeometry('poly', {'flux_bias_line': flux_bias_line})
        self.add_qgeometry('poly', {'flux_bias_line_gap': flux_bias_line_gap}, subtract=True)        

        ############################################################

        # add pins
        port_line_cords = list(draw.shapely.geometry.shape(port_line).coords)
        self.add_pin('flux_bias_line', 
                    port_line_cords, cpw_width)

    def make_charge_line(self):
        """ Adds charge line to fluxonium pocket."""
        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p
        pc = self.p.charge_line_options # parser on charge line options
 
        # define commonly used variables once
        cl_length = pc.cl_length
        cl_sep = pc.cl_sep
        cpw_width = pc.cpw_width
        cpw_gap = pc.cpw_gap

        # For this design the loc_W has to be in 0 but loc_H can be -1 or +1, 
        # For all other direction One can change with the oritation of qubit.
        loc_W = float(pc.loc_W)
        loc_W, loc_H = float(pc.loc_W), float(pc.loc_H)
        if float(loc_W) not in [0] or float(loc_H) not in [-1., +1.]:
            self.logger.info(
                'Warning: Did you mean to define a transmon qubit with loc_W is not 0 and'
                ' loc_H is not +1 or -1 ? Are you sure you want to do this?'
            )

        # Define the geometry
        # Charge Line
    #    h = (p.pocket_height)/2 # The position of charge line
        charge_line = draw.rectangle(cpw_width, cl_length, 0, 0)
        charge_line_gap = draw.rectangle(cpw_width+2*cpw_gap, cl_length+cpw_gap, 0, -cpw_gap/2)

        # Charge Line CPW wire
        port_line = draw.LineString([(-cpw_width/2, cl_length/2),
                                     (cpw_width/2, cl_length/2)])
        
        # Position the charge_line, rotate and translate
        objects = [charge_line, charge_line_gap, port_line]
        objects = draw.scale(objects, 1, loc_H, origin=(0, 0))
        objects = draw.translate(
            objects,
            0,
            loc_H * (p.pocket_height/2 + cl_sep + cl_length))
        objects = draw.rotate_position(objects, p.orientation,
                                       [p.pos_x, p.pos_y])
        [charge_line, charge_line_gap, port_line] = objects

        self.add_qgeometry('poly', {'charge_line': charge_line})
        self.add_qgeometry('poly', {'charge_line_gap': charge_line_gap}, subtract=True)

        ############################################################

        # add pins
        port_line_cords = list(draw.shapely.geometry.shape(port_line).coords)
        port_line_cords = port_line_cords if loc_H==1 else port_line_cords[::-1]
        points = list(port_line_cords)
        self.add_pin('charge_line', 
                    points, cpw_width)

    
    def make_readout_line(self):
        """ Adds readout line to fluxonium pocket."""
        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p
        pr = self.p.readout_line_options # parser on readout line options

        # define commonly used variables once
        pad_sep = pr.pad_sep
        pad_width = pr.pad_width
        pad_height = pr.pad_height
        cpw_width = pr.cpw_width
        cpw_gap = pr.cpw_gap

        # For this design the loc_W has to be in 0 but loc_H can be -1 or +1, 
        # for all other direction One can change with the oritation of qubit.
        loc_W = float(pr.loc_W)
        loc_W, loc_H = float(pr.loc_W), float(pr.loc_H)
        if float(loc_W) not in [0] or float(loc_H) not in [-1., +1.]:
            self.logger.info(
                'Warning: Did you mean to define a transmon qubit with loc_W is not 0 and'
                ' loc_H is not +1 or -1 ? Are you sure you want to do this?'
            )

        # Define the geometry
        # Readout pad
#        h = (p.pocket_height)/2+pad_sep # The position of readout pad or line
        readout_pad = draw.rectangle(pad_width, pad_height, 0, 0)
        readout_pad_circle_left = draw.Point(-pad_width/2, 0).buffer(pad_height/2) # making the pad circle for left side
        readout_pad_circle_right = draw.Point(pad_width/2, 0).buffer(pad_height/2) # making the pad circle for right side

        # Readout pad's gap
        readout_pad_gap = draw.rectangle(pad_width+2*cpw_gap, pad_height+2*cpw_gap, 0, 0)
        readout_pad_gap_circle_left = draw.Point(-(pad_width+2*cpw_gap)/2, 0).buffer((pad_height+2*cpw_gap)/2) # making the pad's gap circle for left side
        readout_pad_gap_circle_right = draw.Point((pad_width+2*cpw_gap)/2, 0).buffer((pad_height+2*cpw_gap)/2) # making the pad's gap circle for right side
        
        # Defining the geometry for readout pad line and it's gap
        readout_line = draw.rectangle(cpw_width, pad_height, 0, pad_height)
        readout_line_gap = draw.rectangle(cpw_width+2*cpw_gap, pad_height, 0, pad_height)

        # Here, we union the readout pad and reaout line and second line excatly the same for the gap
        readout_padNline = draw.union(readout_pad, readout_pad_circle_left, readout_pad_circle_right, readout_line)
        readout_padNline_gap = draw.union(readout_pad_gap, readout_pad_gap_circle_left, readout_pad_gap_circle_right, readout_line_gap)
    
        # Readout Line CPW wire
        port_line = draw.LineString([(cpw_width/2, pad_height*1.5),
                                     (-cpw_width/2, pad_height*1.5)])
       
        # Position the readout, rotate and translate
        objects = [readout_padNline, readout_padNline_gap, port_line]
        objects = draw.scale(objects, 1, loc_H, origin=(0, 0))
        objects = draw.translate(
            objects,
            0,
            loc_H * (p.pocket_height/2 + pad_sep))
        objects = draw.rotate_position(objects, p.orientation,
                                       [p.pos_x, p.pos_y])
        [readout_padNline, readout_padNline_gap, port_line] = objects

        self.add_qgeometry('poly', {'readout_padNline': readout_padNline})
        self.add_qgeometry('poly', {'readout_padNline_gap': readout_padNline_gap}, subtract=True)

        ############################################################

        # add pins
        port_line_cords = list(draw.shapely.geometry.shape(port_line).coords)
        port_line_cords = port_line_cords if loc_H==-1 else port_line_cords[::-1]
        points = list(port_line_cords)
        self.add_pin('readout_line', 
                    points, cpw_width)
