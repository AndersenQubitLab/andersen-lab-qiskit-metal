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

# This class was created by Roald van den Boogaart, Christian Kraglund Andersen

import numpy as np
from qiskit_metal import draw, Dict
from math import *
from qiskit_metal.qlibrary.core import BaseQubit


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
    .. image::
                  | | |  Charge line
                  |___|  
     _______________________________
    |              ___              |
    |             /   \             |   
    |             \   /             |
    |              | |______        |
    |              |_|      |     __|
    |               |       |    |  |
    |               x       |    |  |__
    |               |       |    |_____ Flux bias line
    |              | |______|        __
    |              | |              |
    |             /   \             |
    |             \___/             |
    |_______________________________|
               ___________
              |  _______  |
              | |       | |            
              | |_______| |
              |___  |  ___|
                  | | | 
                    Read out line

    BaseQubit Default Options:
        * pos_x: '0um'
        * pos_y: '0um'
        * connection_pads: Empty Dict -- The dictionary which contains all active connection lines for the qubit.
        * _default_connection_pads: Empty Dict -- The default values for the (if any) connection lines of the qubit.

    Default Options:
        * pos_x: '0um' -- Where the center of the pocket should be located on chip
        * pos_y: '0um' -- Where the center of the pocket should be located on chip
        * width: '1200um' -- The width of the fluxonium qubit along x-axis
        * height: '1200um' -- The height of the fluxonium qubit along y-axis
        * pad_gap: '30um' -- The distance between the two charge islands, which is also the resulting 'length' of the pseudo junction
        * inductor_width: '20um' -- Width of the pseudo junction between the two charge islands (if in doubt, make the same as pad_gap). Really just for simulating in HFSS / other EM software
        * pad_width: '30um' -- The width (x-axis) of the charge island pads
        * pad_height: '300um' -- The size (y-axis) of the charge island pads
        * pad_radius: '80um' -- Radius of the circle at the end of the pads
        * l_width: '5um' -- Width of the kinectic inductor 
        * l_arm_length: '200um' -- Length of the arm of the kinetic inductor along x-axis
        * l_meander_height ='150um' -- Height of the meander of the kinetic inductor along y-axis
        * l_meander_length ='350um' -- Length of the meander of the kinetic inductor 
        * pocket_width: '800um' -- Size of the pocket (cut out in ground) along x-axis
        * pocket_height: '800um' -- Size of the pocket (cut out in ground) along y-axis
        * orientation: '0' -- Degree of qubit rotation
        * flux_bias_line_options=Dict
        *    make_fbl = True -- Boolean to make the flux bias line 
        *    fbl_sep='20um' -- The separation between the flux bias line and the inductor along the x-axis
        *    fbl_height ='50um' -- The height of the flux bias line along the y-axis
        *    cpw_width ='10um' -- The width of the flux bias line
        *    cpw_gap = '10um' -- The dielectric gap width of the flux bias line    
        * charge_line_options=Dict
        *    make_cl = True -- Boolean to make the charge line
        *    cl_sep ='15um' -- The separation between the flux bias line and the pocket the y-axis
        *    cpw_width='10um' -- The width of the charge line
        *    cpw_gap= '10um' -- The dielectric gap width of the charge line
        * readout_line_options=Dict(
        *    make_rol = True -- Boolean to make the readout line
        *    pad_sep='20um' -- The separation between the connection pad and the capacitor pad the y-axis
        *    pad_width='150um' -- Width of the connection pad along the x-axis  
        *    pad_height='50um', -- Height of the connection pad along the y-axis
        *    pad_gap='10um' -- Dielectric gap width of the connection pad
        *    cpw_width='10um', -- The width of the charge line
        *    cpw_gap='10um' -- The dielectric gap width of the readout line
    """

    component_metadata = Dict(short_name='FluxoniumPocket',
                              _qgeometry_table_path='True',
                              _qgeometry_table_poly='True',
                              _qgeometry_table_junction='True')
    """Component metadata"""

    # Default drawing options
    default_options = Dict(
        chip='main',
        pos_x='0um',
        pos_y='0um',
        width='1200um',
        height='1200um',
        pad_gap='30um',
        inductor_width='10um',
        pad_width='15um',
        pad_height='200um',
        pad_radius='50um',
        l_width='1um',
        l_arm_length='50um',
        l_inductance='200nH',
        l_ind_per_square='2nH',
        l_fillet = '5um', 
        L_j = '16.35nH',
        pocket_width='800um',
        pocket_height='800um',
        # 90 has dipole aligned along the +X axis,
        # while 0 has dipole aligned along the +Y axis
        orientation='0',
        flux_bias_line_options=Dict(
            make_fbl = False,
            fbl_sep='50um',
            fbl_height ='50um', 
            cpw_width ='10um',
            cpw_gap = '10um'
        ),
        charge_line_options=Dict(
            make_cl = False,
            cl_sep ='15um',
            cpw_width='10um',
            cpw_gap= '10um'
        ),
        readout_line_options=Dict(
            make_rol = False,
            pad_sep='20um',
            pad_width='150um',
            pad_height='50um',
            pad_gap='10um',
            cpw_width='10um',
            cpw_gap='10um'
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
       
        if self.p.flux_bias_line_options.make_fbl:
            self.make_flux_bias_line()
        if self.p.charge_line_options.make_cl:
            self.make_charge_line()
        if self.p.readout_line_options.make_rol:
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
        L_tot = l_inductance*p.l_width/ind_per_square

        # Draw fillets for inductor meander
        outer_circle = draw.Point(0,0).buffer(p.l_fillet + p.l_width/2)
        inner_circle = draw.Point(0,0).buffer(p.l_fillet - p.l_width/2)
        circle = draw.subtract(outer_circle, inner_circle)
        semicircle = draw.subtract(circle, draw.rectangle(p.l_fillet + p.l_width / 2., 2. * p.l_fillet + p.l_width, - (p.l_fillet + p.l_width / 2.) / 2.,0))
        poly_seg = draw.subtract(semicircle, draw.rectangle(2. * p.l_fillet + p.l_width, p.l_fillet + p.l_width / 2., 0, - (p.l_fillet + p.l_width / 2.) / 2.))
        
        # Place and orientate copies of the poly_seg 
        poly_wire_top = []
        poly_wire_bot = []
        mount_top = draw.translate(draw.rectangle(p.l_arm_length, p.l_width), (p.pad_width + p.l_arm_length)/ 2., +L_tot / 2.+ p.l_fillet)
        mount_bot = draw.translate(draw.rectangle(p.l_arm_length, p.l_width), (p.pad_width + p.l_arm_length)/ 2., -L_tot / 2.- p.l_fillet)
        seg_top = draw.translate(poly_seg, p.l_arm_length + p.pad_width / 2. , L_tot / 2.)
        seg_bot = draw.translate(draw.rotate(poly_seg, -90, origin=(0, 0)), p.l_arm_length + p.pad_width / 2. , -L_tot / 2.)

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
        polys = [rect_jj, pad_top, pad_bot, rect_pk, poly_wire_top, poly_wire_top, inductor]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [rect_jj, pad_top, pad_bot, rect_pk, poly_wire_top, poly_wire_top, inductor] = polys

        # Use the geometry to create Metal qgeometry
        
        self.add_qgeometry('poly', dict(pad_top=pad_top, pad_bot=pad_bot, poly_wire_top = poly_wire_top, poly_wire_bot = poly_wire_bot))
        self.add_qgeometry('poly', dict(rect_pk=rect_pk), subtract=True)
        
        self.add_qgeometry('junction', 
                           dict(inductor=inductor),
                           width=p.l_width,
                           hfss_inductance=str(l_inductance)+'nH')
        self.add_qgeometry('junction',
                           dict(rect_jj=rect_jj),
                           width=p.inductor_width,
                           hfss_inductance=p.L_j)

    def make_flux_bias_line(self):
        """ Adds flux bias line to fluxonium pocket."""
        # Grab option values
        pf = self.p.flux_bias_line_options
        p = self.p

        # draw the flux bias line
        a = (p.pocket_width / 2., pf.fbl_height / 2.)
        b = (p.pad_width / 2. + p.l_arm_length + pf.fbl_sep + (p.l_width + pf.cpw_width) / 2., 0.5 * pf.fbl_height)
        c = (p.pad_width / 2. + p.l_arm_length + pf.fbl_sep + (p.l_width + pf.cpw_width) / 2., 0.5 * -pf.fbl_height)
        d = (p.pocket_width / 2., -pf.fbl_height / 2.)
        e = (p.width / 2., -pf.fbl_height / 2.)

        fbl = draw.LineString([a, b, c, d, e])
        fbl_gap = draw.LineString([d, e])

        # Translate and rotate all elements
        segments = [fbl, fbl_gap]
        segments = draw.rotate(segments, p.orientation, origin=(0, 0))
        segments = draw.translate(segments, p.pos_x, p.pos_y)
        [fbl, fbl_gap] = segments

        self.add_qgeometry('path', {'path':fbl}, layer=1, subtract=False, width=pf.cpw_width)
        self.add_qgeometry('path', {'path':fbl_gap}, layer=1, subtract=True, width=pf.cpw_width + 2. * pf.cpw_gap)

        # add pins
        fbl_pin1 = self.qpin_rotate_translate(d)
        fbl_pin2 = self.qpin_rotate_translate(e)

        self.add_pin('Flux bias line',
                     points=[fbl_pin1, fbl_pin2],
                     width=pf.cpw_width,
                     input_as_norm=True)

    def make_charge_line(self):
        """ Adds charge line to fluxonium pocket."""
        # Grab option values
        pc = self.p.charge_line_options
        p = self.p

        # draw the charge line
        a = p.height / 2.
        b = p.pocket_height / 2. + pc.cl_sep 
        
        cl = draw.box(-pc.cpw_width/2, a, pc.cpw_width/2, b)
        cl_gap = draw.box(-pc.cpw_gap/2-pc.cpw_width/2, a, pc.cpw_gap/2+pc.cpw_width/2, b-pc.cpw_gap)

        # Translate and rotate all elements
        segments = [cl, cl_gap]
        segments = draw.rotate(segments, p.orientation, origin=(0, 0))
        segments = draw.translate(segments, p.pos_x, p.pos_y)
        [cl, cl_gap] = segments

        # # draw the charge line
        self.add_qgeometry('poly', dict(cl=cl))
        self.add_qgeometry('poly', dict(cl_gap=cl_gap), subtract=True)
        
       # add pins
        cl_pin1 = self.qpin_rotate_translate((0,b))
        cl_pin2 = self.qpin_rotate_translate((0,a))

        self.add_pin('Charge line',
                     points=[cl_pin1, cl_pin2],
                     width=pc.cpw_width,
                     input_as_norm=True)
    
    def make_readout_line(self):
        """ Adds readout line to fluxonium pocket."""
        # Grab option values
        pr = self.p.readout_line_options
        p = self.p
        
        # Draw the readout pad
        pad = draw.rectangle(pr.pad_width, pr.pad_height)
        cutout_pad = draw.rectangle(pr.pad_width + 2. * pr.pad_gap, pr.pad_height + 2. * pr.pad_gap)
        pad = draw.translate(pad, 0, -(pr.pad_height + p.pocket_height) / 2. - pr.pad_sep - pr.pad_gap)
        cutout_pad = draw.translate(cutout_pad, 0, -(pr.pad_height + p.pocket_height) / 2. - pr.pad_sep - pr.pad_gap)

        # Draw the coplanar waveguide
        a = (0, -p.height / 2.)
        b = (0, - pr.pad_height - p.pocket_height / 2. - pr.pad_sep - pr.pad_gap )
        
        cpw_readout = draw.LineString([a,b])

        # Translate and rotate all elements
        segments = [pad, cutout_pad, cpw_readout]
        segments = draw.rotate(segments, p.orientation, origin=(0, 0))
        segments = draw.translate(segments, p.pos_x, p.pos_y)
        [pad, cutout_pad, cpw_readout] = segments

        self.add_qgeometry('poly', dict(pad=pad))
        self.add_qgeometry('poly', dict(cutout_pad=cutout_pad), subtract=True)
        self.add_qgeometry('path',{'path':cpw_readout}, layer=1, subtract=False, width=pr.cpw_width)
        self.add_qgeometry('path',{'path':cpw_readout}, layer=1, subtract=True, width=pr.cpw_width + 2. * pr.cpw_gap)

        # add pins
        ro_pin1 = self.qpin_rotate_translate(b)
        ro_pin2 = self.qpin_rotate_translate(a)

        self.add_pin('Readout line',
                     points=[ro_pin1, ro_pin2],
                     width=pr.cpw_width,
                     input_as_norm=True)

    def qpin_rotate_translate(self, x):
        p = self.p
        y = list(x)
        z = [0.0, 0.0]
        z[0] = y[0] * cos(p.orientation * 3.14159 / 180) - y[1] * sin(
                p.orientation * 3.14159 / 180)
        z[1] = y[0] * sin(p.orientation * 3.14159 / 180) + y[1] * cos(
                p.orientation * 3.14159 / 180)
        z[0] = z[0] + p.pos_x
        z[1] = z[1] + p.pos_y
        x = (z[0], z[1])
        return x