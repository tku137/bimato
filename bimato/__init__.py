# Copyright (C) 2022 Tony Fischer
#
# This file is part of Bio Matrix Topology (BiMaTo).
#
# Bio Matrix Topology (BiMaTo) is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Bio Matrix Topology (BiMaTo) is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Bio Matrix Topology (BiMaTo).  If not, see <http://www.gnu.org/licenses/>.

'''The module :mod:`bimato` is the top-level namespace. It contains the following sub-modules:

- :mod:`bimato.core`: fundamental functions to the bimato project, mainly the custom binarization
- :mod:`bimato.poresize`: pore-size algorithm published at https://www.nature.com/articles/s41598-019-44764-5
- :mod:`bimato.network`: algorithms to analyze the network structure
- :mod:`bimato.utils`: various statistical and other utility functions

See documentation for a general user guide.
'''


__author__ = "Tony Fischer (tku137)"
__copyright__ = "Copyright 2022, Tony Fischer (tku137)"
__license__ = "GPLv3"
__email__ = "tonyfischer@mailbox.org"
__status__ = "Development"
__version__ = "2022.1.2"
__credits__ = ["Tony Fischer (tku137)", "Alexander Hayn"]


from bimato.api import *  # noqa: F401,F403
