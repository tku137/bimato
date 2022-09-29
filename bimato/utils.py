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
# TODO: write helpers module docstring
docstring of module
"""


__author__ = "Tony Fischer (tku137)"
__copyright__ = "Copyright 2022, Tony Fischer (tku137)"
__license__ = "GPLv3"
__email__ =  "tonyfischer@mailbox.org"
__status__ = "Development"
__version__ = "2022.1"
__credits__ = ["Tony Fischer (tku137)", "Alexander Hayn"]


import numpy as np
import pandas as pd
from skimage.restoration import denoise_tv_chambolle


def read_data(lif_stack):
    """
    Read image data from Leica lif file.

    Parameters
    ----------
    lif_stack : LifImage
        Instance of LeicaBioFormats.LifImage

    Returns
    -------
    numpy.ndarray
        Image data in numerical form

    """
    # TODO: adapt readlif
    z_size = int(lif_stack.info['SizeZ'])
    data = lif_stack.read_image(z=(0, z_size)).squeeze()
    return data


def normalize_image(image, scale):
    """
    Standard function to normalize image data values to [0,1].

    Parameters
    ----------
    image : numpy.ndarray
        Input image data
    scale : int or float
        Scaling factor, e.g. 255 would scale to 8bit integer matrices

    Returns
    -------
    numpy.ndarray
        Normalized image data

    """
    im_min = image.min()
    im_max = image.max()
    return scale * ((image - im_min) / (im_max - im_min))


def get_weight_factor(size):
    """
    Calculates the weight factor for TV-denoising using scikit-image. Values have been determined empirically.

    Notes
    -----
    Formulas of fitting empirical data:

    - mean values: 3E-08x^2 + 8E-06x + 0,0037

    - more denoising: 3E-08x^2 + 6E-06x + 0,0127

    - linear fit was: 6E-05x - 0,0215

    Parameters
    ----------
    size: int
        Relevant size of image data
    Returns
    -------
    float
        TV-denoising weight factor
    """
    weight = 3e-8 * (size ** 2) + 8e-6 * size + 0.0037
    return weight


# TODO: inplace possible?
def denoise_image(image, denoise_weight=None):
    """
    Apply TV-denoising on input image data. Denoising weight factor is determined automatically, except explicitly
    specified

    Parameters
    ----------
    image : numpy.ndarray
        Input raw image data
    denoise_weight : int or float, optional
        Weight factor for scikit image TV-denoising

    Returns
    -------
    numpy.ndarray
        Denoised image data

    """
    # get weight factor
    small_edge_size = np.min(image.shape[:2])  # size of the smaller dim of an x-y image
    if denoise_weight:
        weight = denoise_weight
    else:
        weight = get_weight_factor(small_edge_size)

    # preallocate image array
    image_denoised = np.zeros_like(image, dtype=np.dtype('uint8'))

    # denoise per-plane for better results
    for z in range(image_denoised.shape[2]):
        # denoise_tv_chambolle always gives float64 in range [0,1]
        # rescale to unit8 [0,255] range
        image_denoised[:, :, z] = denoise_tv_chambolle(image[:, :, z], weight=weight) * 255

    #
    return image_denoised


def calc_weighted_median(df, data_col, weight_col):
    """
    Calculates the weighted median of specified DataFrame column with respect to another specified column-

    Parameters
    ----------
    df: pandas.DataFrame
        Input DataFrame that came from initial pore detection function
    data_col: str
        Name of column from which the median is calculated
    weight_col: str
        Name of column that is used as weight

    Returns
    -------
    pandas.Series
        Pandas Series with calculated weighted median

    """
    tmp_df = df.sort_values(data_col)
    cumsum = tmp_df[weight_col].values.cumsum()
    cutoff = tmp_df[weight_col].values.sum() / 2.0
    return tmp_df[data_col][cumsum >= cutoff].iloc[0]


def get_voxel_volume(lif_stack):
    vx = float(lif_stack.info['PhysicalSizeX'])
    vy = float(lif_stack.info['PhysicalSizeY'])
    vz = float(lif_stack.info['PhysicalSizeZ'])
    return vx * vy * vz


# construct list of intervals as numpy array
def get_interval_list(part_size, number_parts):
    return [(np.array([0, 1]) + n) * part_size for n in range(number_parts)]


def get_intervals(lif_stack, part_size):
    # only depends on the part_size, the part numbers are calculated based on micron sizes
    number_parts_x, number_parts_y, number_parts_z = get_possible_part_numbers(lif_stack, part_size)
    part_size_x_px, part_size_y_px, part_size_z_px = get_part_size_px(lif_stack, number_parts_x, number_parts_y,
                                                                      number_parts_z)

    intervals_x_px = get_interval_list(part_size_x_px, number_parts_x)
    intervals_y_px = get_interval_list(part_size_y_px, number_parts_y)
    intervals_z_px = get_interval_list(part_size_z_px, number_parts_z)

    return intervals_x_px, intervals_y_px, intervals_z_px


def get_cube_size_in_micron(lif_stack):
    sx = int(lif_stack.info['SizeX']) * float(lif_stack.info['PhysicalSizeX'])
    sy = int(lif_stack.info['SizeY']) * float(lif_stack.info['PhysicalSizeY'])
    sz = int(lif_stack.info['SizeZ']) * float(lif_stack.info['PhysicalSizeZ'])
    return sx, sy, sz


def get_possible_part_numbers(lif_stack, part_size):
    sx, sy, sz = get_cube_size_in_micron(lif_stack)
    nx = np.floor(sx / part_size).astype('int')
    ny = np.floor(sy / part_size).astype('int')
    nz = np.floor(sz / part_size).astype('int')
    return nx, ny, nz


def get_part_size_px(lif_stack, nx, ny, nz):
    px = np.floor(int(lif_stack.info['SizeX']) / nx).astype('int')
    py = np.floor(int(lif_stack.info['SizeY']) / ny).astype('int')
    pz = np.floor(int(lif_stack.info['SizeZ']) / nz).astype('int')
    return px, py, pz
