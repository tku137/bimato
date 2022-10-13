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

'''The module :mod:`bimato.core` contains some fundamental functions for the `bimato` project. Mainly it contains the
custom binarization technique tailored towards 3D biompolymer network images and others.

See documentation for a general user guide.
'''


__author__ = "Tony Fischer (tku137)"
__copyright__ = "Copyright 2022, Tony Fischer (tku137)"
__license__ = "GPLv3"
__email__ = "tonyfischer@mailbox.org"
__status__ = "Development"
__version__ = "2022.1.2"
__credits__ = ["Tony Fischer (tku137)", "Alexander Hayn"]


from skimage.filters import threshold_local  # pylint: disable=no-name-in-module
from skimage.morphology import binary_closing, disk, ball, remove_small_objects
from scipy.ndimage.morphology import distance_transform_edt as edt

import numpy as np


def get_binary(data, bw_closing_frac=1024, min_obj_frac=256, win_frac=8, thresh_shift_percentage=0.01, morph_3d=True):
    '''Function to get a segmentation of input biopolymer network image data. Parameters are optional, but should be
    considered. Optimal values can vary for drastically differing image data.

    Notes
    -----
    All Parameters are determined automatically depending on input image properties. However, they can be influenced by
    specifying a factor.

    Parameters
    ----------
    data : numpy.ndarray
        Input image data
    bw_closing_frac : int
        Factor for calculation of the ball size during binary closing morphological operations
    min_obj_frac : int
        Factor for determining the minimum object size
    win_frac : int
        Factor for determining the local adaptive threshold window size
    thresh_shift_percentage : float
        Percentage factor in [0,1] specifying the amount the threshold should be shifted
    morph_3d : bool
        If morphological operations should be applied in 2D or 3D, default is 3D

    Returns
    -------
    numpy.ndarray
        Binary segmentation output

    '''

    # calc the threshold shift
    thresh_shift = thresh_shift_percentage * data.max()

    # determine optimal parameters
    small_edge_size = np.min(data.shape[:2])  # size of the smaller dim of an x-y image
    threshold_win = 2 * int(small_edge_size / (win_frac * 2)) + 1  # threshold_local window
    bw_closing_size = int(np.ceil(small_edge_size / bw_closing_frac))  # ball size for morphological closing

    min_obj_exp = 3 if morph_3d else 2  # wanna morph in 3D or 2D?
    min_obj = int(np.ceil(small_edge_size / min_obj_frac) ** min_obj_exp)  # objects must be this large

    # per plane adaptive threshold
    binary = np.zeros_like(data, dtype=np.bool)
    for z in range(data.shape[2]):
        # get the local threshold map
        t = threshold_local(data[..., z], threshold_win)
        # threshold the image with it
        binary[..., z] = data[..., z] > (t + thresh_shift)
        if not morph_3d:
            # do morphological operation in 2D, because 3D would cut of the bw_closing_size size in z direction!
            binary[..., z] = binary_closing(binary[..., z], footprint=disk(bw_closing_size))

    # morphological closing in 3D
    if morph_3d:
        binary = binary_closing(binary, footprint=ball(bw_closing_size))

    # remove small objects and holes
    binary = remove_small_objects(binary, min_size=min_obj)

    #
    return binary


def get_edm(binary, sampling):
    '''It takes a binary image and a dictionary of sampling rates and returns the Euclidean distance map

    Parameters
    ----------
    binary : numpy.ndarray
        a 3D binary array
    sampling : dict
        a dictionary with keys 'x', 'y', and 'z' that contain the voxel size in each dimension

    Returns
    -------
    numpy.ndarray
        The Euclidean distance map (EDM) of the binary image.

    '''
    edm = edt(
        np.logical_not(binary),
        sampling=[sampling['x'], sampling['y'], sampling['z']]
    )
    return edm
