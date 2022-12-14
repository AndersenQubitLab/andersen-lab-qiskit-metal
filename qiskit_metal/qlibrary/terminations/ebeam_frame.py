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
from turtle import width
from qiskit_metal import draw
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent

# Define class and options for the launch geometry

class Frame(QComponent):
    """The frame will be on the device as 'positive'.

    Inherits 'QComponent' class.

    .. image:
        Frame.png

    .. meta::
        Frame

    Sketch:
        Below is a sketch of the frame
        ::
             _______________
            |  ___________  |                 y                             |
            | |           | |                 ^                         2   |   1
            | |           | |                 |                             |
            | |           | |                 |                      _______|_______
            | |           | |                 |                             | 
            | |___________| |                 |------> x               3    |   4
            |_______________|                                               |

    .. image::
        Frame.png

    Default Options:
        * frame_w: '9mm' -- the width of the frame. That's our device size; we decided to be 9*9mm square
        * frame_h: '9mm' -- the width of the frame. That's our device size; we decided to be 9*9mm square
        * f_width: '70um' -- the thickness of the frame line can be arranged accordingly dicing blade thickness 
    """

    default_options = Dict(
        frame_w = '9mm',
        frame_h = '9mm',
        f_width = '70um',
        )
    """Default options"""

    TOOLTIP = """Frame on the chip."""

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
