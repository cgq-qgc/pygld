# -*- coding: utf-8 -*-

"""
PyGLD License Agreement (MIT License)
--------------------------------------

Copyright (c) Louis Lamarche
Copyright (c) Jean-Sébastien Gosselin
https://github.com/jnsebgosselin/pygld

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
"""

from pygld.properties import HeatCarrierFluid
# from geothermie.thermal_R import *
# from geothermie.design_glhe import (calcul_Tp_ils, calcul_Tp_fls,
#                                     design_borefield)
# from geothermie.g_function import (G_function_ICS, g_function_cooper,
#                                    g_function_bernier, G_function_ILS,
#                                    G_function_FLS, G_function_lamarche)
# from geothermie.geometrie import (TLB, calcul_sdr, plot_borehole,
#                                   FigBoreholeGeo)

# from geothermie.glhe import OpenGLDesignerBase

import os
import sys

version_info = (0, 1, 0, "dev")
__version__ = '.'.join(map(str, version_info))
__appname__ = 'PyGLD'
__namever__ = __appname__ + " " + __version__
__date__ = '15/12/2017'
__project_url__ = "https://github.com/jnsebgosselin/pygld"
__releases_url__ = __project_url__ + "/releases"
__releases_api__ = "https://api.github.com/repos/jnsebgosselin/pygld/releases"


def is_frozen():
    """
    Return whether the application is running from a frozen exe or if it
    is running from the Python source files.

    See: https://stackoverflow.com/a/42615559/4481445
    """
    return getattr(sys, 'frozen', False)


if is_frozen():
    __rootdir__ = sys._MEIPASS
else:
    __rootdir__ = os.path.dirname(os.path.realpath(__file__))
