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

#  This a launch structure used on BlueJayV2, used for wire bonding
#  There is no CPW tee attached to this p#

# This class was created by Figen YILMAZ

# Imports required for drawing

# import numpy as np # (currently not used, may be needed later for component customization)
from turtle import width
from qiskit_metal import draw
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent

# Define class and options for the launch geometry


class Frame(QComponent):
    """Launch pad to feed/read signals to/from the chip.

    Inherits 'QComponent' class.

    .. image:
        LaunchpadWirebondCoupled.png

    .. meta::
        Launchpad Wirebond Coupled

    Creates a 50 ohm launch pad with a ground pocket cutout.
    Limited but expandable parameters to control the launchpad polygons.
    The (0,0) point is the center of the necking of the launch tip.
    The pin attaches directly to the built in lead length at its midpoint
    This launch has an inductive coupler section.

    Pocket and pad:
        Pocket and launch pad geometries are currently fixed.
        (0,0) point is the midpoint of the necking of the launch tip.
        Pocket is a negative shape that is cut out of the ground plane

    Values (unless noted) are strings with units included, (e.g., '30um')

    Sketch:
        Below is a sketch of one marker
        ::

             __________                  y                                 |
            |          |                  ^                            2   |   1
            |          |                  |                         _______|_______
            |          |                  |                                | 
            |__________|                  |------> x                  3    |   4
      

    .. image::
        Markers.png

    Default Options:
        * marker: '20um' -- center trace width of the marker
    """

    default_options = Dict(
        chip='main',
        pos_x = '0.0 mm',
        pos_y = '0.0 mm',
        frame_w = '9mm',
        frame_h = '9mm',
        f_width = '70um',
        )
    """Default options"""

    TOOLTIP = """Markers for the chip."""

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""

        p = self.p
        frame_w = p.frame_w
        frame_h = p.frame_h
        f_width = p.f_width

        # Creates a frame around the device, size of 9*9mm
        line_top = draw.rectangle(frame_w, f_width, 0, 4.5)
        line_left = draw.rectangle(f_width, frame_h, 4.5, 0)
        line_right = draw.rectangle(f_width, frame_h, -4.5, 0)
        line_bottom = draw.rectangle(frame_w, f_width, 0, -4.5)

        frame = draw.union(line_top, line_left, line_right, line_bottom)

        # Create polygon object list
        polys = [frame]

        # Rotates and translates all the objects as requested. Uses package functions
        # in 'draw_utility' for easy rotation/translation
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, xoff=p.pos_x, yoff=p.pos_y)
        [frame] = polys

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(frame=frame),
                           layer=p.layer)