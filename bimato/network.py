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

"""
docstring of module
"""


__author__ = "Tony Fischer (tku137)"
__copyright__ = "Copyright 2022, Tony Fischer (tku137)"
__license__ = "GPLv3"
__email__ = "tonyfischer@mailbox.org"
__status__ = "Development"
__version__ = "2022.1.2"
__credits__ = ["Tony Fischer (tku137)", "Alexander Hayn"]

import numpy as np
from skimage.morphology import ball, remove_small_objects, skeletonize_3d
from skimage.measure import label


def get_network_skeleton(binary, min_size=10):
    """
    Calculates the network skeletonization of an input network binary.

    Parameters
    ----------
    binary : numpy.ndarray
        Input binary segmentation of typed `bool`
    min_size : int
        Minimum size of objects in `px`

    Returns
    -------
    numpy.ndarray
        Skeleton output

    """
    sk = skeletonize_3d(binary).astype(np.bool)
    sk_refined = remove_small_objects(sk, min_size=min_size, connectivity=3)
    return sk_refined


def classify_network(skeleton):
    """
    Classifies each binary fibril pixel of input skeleton depending on number of neighbouring pixels.

    Class 0 is fluid phase. Class 2 is a fibril with class 1 being a fibril ending. Classes higher than 2 mean that this
    is a node, while denoting the number of connections outgoing.

    Warning: the center pixel is always 1, so the sum of each cut out fibril environment needs to be reduced by 1

    Parameters
    ----------
    skeleton : numpy.ndarray
        Input skeleton

    Returns
    -------
    numpy.ndarray
        Classification matrix

    """
    fibril_coords = np.argwhere(skeleton > 0)

    # filter border coordinates
    x_size, y_size, z_size = skeleton.shape
    if not (x_size == y_size == z_size):
        raise ValueError("image dimensions are not equal")

    fibril_coords = fibril_coords[np.all(fibril_coords > 0, axis=1), :]
    fibril_coords = fibril_coords[np.all(fibril_coords < (x_size-1), axis=1), :]

    fibril_class = np.zeros_like(skeleton, dtype=np.uint8)
    for fib in fibril_coords:
        x, y, z = fib
        fib_cut = skeleton[x - 1:x + 2, y - 1:y + 2, z - 1:z + 2]
        fibril_class[x, y, z] = np.sum(fib_cut) - 1

    return fibril_class


def get_nodes(fibril_class, se_radius=1):
    # node analysis
    node_coords = np.argwhere(fibril_class > 2)  # coordinates of 'nodes' aka fibril_class > 2

    node_stack = np.zeros_like(fibril_class)  # empty stack
    for row in node_coords:  # for each node coordinate, draw ball
        x, y, z = row.astype('int')
        node_indices = np.argwhere(ball(se_radius))
        node_indices[:, 0] += x - se_radius
        node_indices[:, 1] += y - se_radius
        node_indices[:, 2] += z - se_radius

        # filter edge coords
        node_indices = node_indices[node_indices[:, 0] < node_stack.shape[0]]
        node_indices = node_indices[node_indices[:, 1] < node_stack.shape[1]]
        node_indices = node_indices[node_indices[:, 2] < node_stack.shape[2]]

        # binarize node stack
        node_stack[node_indices[:, 0], node_indices[:, 1], node_indices[:, 2]] = 1

    return node_stack


def get_single_fibrils_labeled(skeleton, node_stack):
    """
    Returns a label matrix with single fibrils. Each integer number in the label matrix represents pixels of a single,
    individual fibril.

    Parameters
    ----------
    fibril_class : numpy.ndarray
        Matrix with classified skeleton

    Returns
    -------
    numpy.ndarray
        Labelled matrix with single fibrils

    """
    # skeleton * inverse of node stack yields single fibril pixels
    single_fibrils = label(skeleton * np.logical_not(node_stack))

    return single_fibrils


def get_edge_fibril_ids(fibril_labels):
    edge_fibril_ids = np.unique(np.concatenate((
        np.unique(fibril_labels[:2, :, :]),  # front
        np.unique(fibril_labels[-2:, :, :]),  # back

        np.unique(fibril_labels[:, :2, :]),  # left
        np.unique(fibril_labels[:, -2:, :]),  # right

        np.unique(fibril_labels[:, :, :2]),  # bottom
        np.unique(fibril_labels[:, :, -2:]),  # top
    )))
    return edge_fibril_ids
