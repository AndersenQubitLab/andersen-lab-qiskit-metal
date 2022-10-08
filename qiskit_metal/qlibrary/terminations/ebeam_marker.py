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


class Markers(QComponent):
    """ * to/from the chip.

    Inherits 'QComponent' class.

    .. image:
        Frame.png

    .. meta::
        Frame

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
        marker_sep = '30um',
        marker_w='20um',
        marker_h='20um',
         )
    """Default options"""

    TOOLTIP = """Markers for the chip."""

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""

        p = self.p
        marker_w = p.marker_w
        marker_h = p.marker_h
        pos_x = p.pos_x
        pos_y = p.pos_y
        marker_sep = p.marker_sep

        # Square shape markers for E-Beam marker search, check ...
        marker_1 = draw.rectangle(marker_w, marker_h, pos_x, pos_y)
        marker_2 = draw.rectangle(marker_w, marker_h, pos_x - marker_sep, pos_y)
        marker_3 = draw.rectangle(marker_w, marker_h, pos_x - marker_sep, pos_y - marker_sep)
        marker_4 = draw.rectangle(marker_w, marker_h, pos_x, pos_y - marker_sep)

        marker = draw.union(marker_1, marker_2, marker_3, marker_4)

        # Create polygon object list
        polys = [marker]

        # Rotates and translates all the objects as requested. Uses package functions
        # in 'draw_utility' for easy rotation/translation
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, xoff=p.pos_x, yoff=p.pos_y)
        [marker] = polys

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(marker=marker),
                           layer=p.layer)
