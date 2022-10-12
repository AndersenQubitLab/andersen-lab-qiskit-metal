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

# Imports required for drawing

# import numpy as np # (currently not used, may be needed later for component customization)
from qiskit_metal import draw
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent

# Define class and options for the markers geometry

class Markers(QComponent):
    """Markers for the EBPG marker search; positive markers.

    Inherits 'QComponent' class.

    .. image:
        Markers.png

    .. meta::
        Markers

    Creates 4 sets of markers for EBPG marker search. 
    The (0,0) point is the centre of these 4 markers.
    For how their named one can check the drawing below. 
    2 Dimensional Cartesian Coordinate System is the way how they are named.

    Sketch:
        Below is a sketch of one the marker set.
        ::

             __    __          y                            |
            |__|  |__|          ^                       2   |   1
                                |                    _______|_______
             __    __           |                           | 
            |__|  |__|          |------> x              3   |   4
      

    .. image::
        Markers.png

    Default Options:
        * pos_x = '0.0mm' -- the x position of the marker
        * pos_y = '0.0mm' -- the y position of the marker
        * marker_sep = '20um' -- the distance between each marker
        * marker_w = '20um' -- marker width
        * marker_h = '20um' -- marker height
    """

    default_options = Dict(
        pos_x = '0.0mm',
        pos_y = '0.0mm',
        marker_sep = '20um',
        marker_w = '20um',
        marker_h = '20um',
         )
    """Default options"""

    TOOLTIP = """Markers for the chip"""

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""

        p = self.p
        pos_x = p.pos_x
        pos_y = p.pos_y
        marker_sep = p.marker_sep
        marker_w = p.marker_w
        marker_h = p.marker_h
        
        # Square shape markers for E-Beam marker search, see the drawing above for their positions
        # on the cartesian coordinate system
        marker_1 = draw.rectangle(marker_w, marker_h, pos_x/2 + marker_sep, pos_y/2 + marker_sep)
        marker_2 = draw.rectangle(marker_w, marker_h, pos_x/2 - marker_sep, pos_y/2 + marker_sep)
        marker_3 = draw.rectangle(marker_w, marker_h, pos_x/2 + marker_sep, pos_y/2 - marker_sep)
        marker_4 = draw.rectangle(marker_w, marker_h, pos_x/2 - marker_sep, pos_y/2 - marker_sep)

        markers = draw.union(marker_1, marker_2, marker_3, marker_4)

        # Create polygon object list
        polys = [markers]
        # Rotates and translates all the objects as requested. Uses package functions
        # in 'draw_utility' for easy rotation/translation
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, xoff=pos_x/2, yoff=pos_y/2)
        [markers] = polys

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(markers=markers),
                           layer=p.layer)
