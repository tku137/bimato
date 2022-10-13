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

'''The module :mod:`bimato.poresize` contains the complete pore-size algorithm published at
https://www.nature.com/articles/s41598-019-44764-5. The most important is :func:`bimato.poresize.get_pore_sizes` and
:func:`bimato.poresize.get_fiber_thickness`.

See documentation for a general user guide.
'''


__author__ = "Tony Fischer (tku137)"
__copyright__ = "Copyright 2022, Tony Fischer (tku137)"
__license__ = "GPLv3"
__email__ = "tonyfischer@mailbox.org"
__status__ = "Development"
__version__ = "2022.1.2"
__credits__ = ["Tony Fischer (tku137)", "Alexander Hayn"]


from itertools import product
from skimage.feature import peak_local_max
from skimage.filters import gaussian  # pylint: disable=no-name-in-module
from skimage.morphology import ball
import numpy as np
import pandas as pd
from .core import get_edm
from .utils import get_voxel_volume, get_intervals


def calc_free_pore_space(df):
    '''Calculates the fraction of pore volume to sample volume.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame that came from initial pore detection function

    Returns
    -------
    pandas.Series
        Pandas Series with calculated free pore space

    '''
    sum_pore_vol = df.loc[:, "Pore Volume [µm³]"].sum()
    cube_vol = df.loc[:, "Cube Volume [µm³]"].iloc[0]
    return sum_pore_vol / cube_vol


def get_pore_sizes(binary, sampling, sigma_frac=128, residual_pore_detection=False):
    '''This function takes a binary image of a 3D biopolymer network image and calculates the pore sizes given the
    pixel-micron conversion.

    Parameters
    ----------
    binary
        The binary image containing the segmentation of fluid (=0) and polymer phase (=1)
    sampling
        The physical size of each voxel in the image
    sigma_frac, optional
        The sigma parameter for the gaussian filter. The smaller the value, the more pores are detected
    residual_pore_detection, optional
        If True, the algorithm will try to detect even smaller and obscure pores

    Returns
    -------
        A dataframe with the following columns, each row represents data for a single pore:
        x [px] : x coordinate in pixel
        y [px] : y coordinate in pixel
        z [px] : z coordinate in pixel
        Diameter [µm] : pore diameter
        Residual Degree : number of iterative pore detections, 0 for first step, 1 for second
        Size x [px] : total x size of 3D image
        Size y [px] : total y size of 3D image
        Size z [px] : total z size of 3D image
        PhysicalSize x : pixel-micron conversion factor used
        PhysicalSize y : pixel-micron conversion factor used
        PhysicalSize z : pixel-micron conversion factor used
        Number Of Pores : total number of detected pores in sample
        Cube Volume [µm³] : physical volume of 3D image
        Single Pore Volume [µm³] : volume of detected pore
        Real Pore Volume [µm³] : physical volume occupied by detected pores, without overlapping of individual pores
        Collagen Volume [µm³] : volume of the polymer phase
        Fluid Volume [µm³] : volume of the fluid phase
        Single Pore Volume Fraction : space occupied in relation to the sample volume
        Real Pore Volume Fraction : space occupied in relation to the sample volume
        Collagen Volume Fraction : space occupied in relation to the sample volume
        Fluid Volume Fraction : space occupied in relation to the sample volume
        Zeta Single Pores : pore diameter scaled with its corresponding column fraction
        Zeta Real Pores : pore diameter scaled with its corresponding column fraction
        Zeta Collagen : pore diameter scaled with its corresponding column fraction
        Zeta Fluid : pore diameter scaled with its corresponding column fraction
        Pseudo Pore Diameter [µm] : theoretical pore diameter calculated by assuming that all detected pores would occupy the fluid-phase fully

    '''
    # calc fluid phase volume
    voxel_volume = get_voxel_volume(sampling)
    voxel_radius = np.cbrt((3*voxel_volume)/(4*np.pi))

    # calculate the EDM
    small_edge_size = np.min(binary.shape[:2])  # size of the smaller dim of an x-y image
    sigma = small_edge_size / sigma_frac
    edm = get_edm(binary, sampling)  # edm has microns as values

    # get pore coordinates
    pore_coords = peak_local_max(gaussian(edm, sigma=sigma))

    # construct DataFrame with diameter and other infos
    df_pores = pd.DataFrame(pore_coords, columns=['x [px]', 'y [px]', 'z [px]'])  # coords are px coords!
    # we take values from edm, so the result is microns again!
    df_pores['Diameter [µm]'] = edm[pore_coords[:, 0], pore_coords[:, 1], pore_coords[:, 2]] * 2  # Radius! So x2!
    # TODO: multi-degree
    df_pores['Residual Degree'] = 0
    # calc real pore volume
    pore_stack = np.zeros_like(binary)

    for _, pore in df_pores.iterrows():
        x, y, z, d, _ = pore.astype(np.int)  # need integers to index the array
        r = int((d / 2) / float(voxel_radius))  # convert to pixel-radius

        # construct coordinate arrays of a sphere at (x,y,z) with radius r
        pore_indices = np.argwhere(ball(r))
        pore_indices[:, 0] += x - r
        pore_indices[:, 1] += y - r
        pore_indices[:, 2] += z - r

        # filter out of range pore voxel
        pore_indices = pore_indices[pore_indices[:, 0] < pore_stack.shape[0]]
        pore_indices = pore_indices[pore_indices[:, 1] < pore_stack.shape[1]]
        pore_indices = pore_indices[pore_indices[:, 2] < pore_stack.shape[2]]

        # draw sphere
        pore_stack[pore_indices[:, 0], pore_indices[:, 1], pore_indices[:, 2]] = 1

    # residual poresize detection
    if residual_pore_detection:

        # this is collagen plus already detected pores
        residual_fluid = binary + pore_stack

        # calculate the EDM
        edm_residual = get_edm(residual_fluid, sampling)  # edm has microns as values

        # get pore coordinates
        pore_coords_residual = peak_local_max(gaussian(edm_residual, sigma=sigma))

        # construct DataFrame with diameter and other infos
        df_pores_residual = pd.DataFrame(pore_coords_residual, columns=['x [px]', 'y [px]', 'z [px]'])  # coords are px!
        # we take values from edm, so the result is microns again!
        df_pores_residual['Diameter [µm]'] = edm_residual[
                                                 pore_coords_residual[:, 0],
                                                 pore_coords_residual[:, 1],
                                                 pore_coords_residual[:, 2]] * 2  # Radius! So x2!
        df_pores_residual['Residual Degree'] = 1

        # calc real pore volume
        for _, pore in df_pores_residual.iterrows():  # TODO: refactor to function
            x, y, z, d, _ = pore.astype(np.int)  # need integers to index the array
            r = int((d / 2) / float(voxel_radius))  # convert to pixel-radius
            pore_indices = np.argwhere(ball(r))
            pore_indices[:, 0] += x - r
            pore_indices[:, 1] += y - r
            pore_indices[:, 2] += z - r
            pore_indices = pore_indices[pore_indices[:, 0] < pore_stack.shape[0]]
            pore_indices = pore_indices[pore_indices[:, 1] < pore_stack.shape[1]]
            pore_indices = pore_indices[pore_indices[:, 2] < pore_stack.shape[2]]
            pore_stack[pore_indices[:, 0], pore_indices[:, 1], pore_indices[:, 2]] = 1

        df_pores = pd.concat([df_pores, df_pores_residual], axis=0)

    # parsed meta-data
    # TODO: meta-data in comment lines above DataFrame
    size_x, size_y, size_z = binary.shape
    df_pores['Size x [px]'] = size_x
    df_pores['Size y [px]'] = size_y
    df_pores['Size z [px]'] = size_z
    df_pores['PhysicalSize x'] = sampling['x']
    df_pores['PhysicalSize y'] = sampling['y']
    df_pores['PhysicalSize z'] = sampling['z']

    # number of pores
    df_pores["Number Of Pores"] = df_pores.shape[0]

    # sample volume (the image cube)
    df_pores["Cube Volume [µm³]"] = (df_pores['Size x [px]'] * df_pores['Size y [px]'] * df_pores[
        'Size z [px]']) * voxel_volume

    # volume of each individual pore, the volumes might overlap!
    df_pores['Single Pore Volume [µm³]'] = (4 / 3) * np.pi * ((df_pores["Diameter [µm]"] / 2) ** 3)

    # calculate the real volume that all pores occupy
    df_pores['Real Pore Volume [µm³]'] = pore_stack.sum() * voxel_volume

    # calculate volume of collagen and fluid phase
    df_pores['Collagen Volume [µm³]'] = binary.sum() * voxel_volume
    df_pores['Fluid Volume [µm³]'] = df_pores["Cube Volume [µm³]"] - df_pores['Collagen Volume [µm³]']

    # how much space each component volume occupies in relation to the sample (cube) volume
    df_pores["Single Pore Volume Fraction"] = df_pores["Single Pore Volume [µm³]"] / df_pores["Cube Volume [µm³]"]
    df_pores["Real Pore Volume Fraction"] = df_pores["Real Pore Volume [µm³]"] / df_pores["Cube Volume [µm³]"]
    df_pores["Collagen Volume Fraction"] = df_pores["Collagen Volume [µm³]"] / df_pores["Cube Volume [µm³]"]
    df_pores["Fluid Volume Fraction"] = df_pores["Fluid Volume [µm³]"] / df_pores["Cube Volume [µm³]"]

    # the actual scaled pore-sizes Zeta. this is the pore diameter scaled with its corresponding column fraction.
    df_pores["Zeta Single Pores"] = df_pores["Diameter [µm]"] * df_pores["Single Pore Volume Fraction"]
    df_pores["Zeta Real Pores"] = df_pores["Diameter [µm]"] * df_pores["Real Pore Volume Fraction"]
    df_pores["Zeta Collagen"] = df_pores["Diameter [µm]"] * df_pores["Collagen Volume Fraction"]
    df_pores["Zeta Fluid"] = df_pores["Diameter [µm]"] * df_pores["Fluid Volume Fraction"]

    # pseudo pore diameter
    df_pores["Pseudo Pore Diameter [µm]"] = np.cbrt(
        (6 * (df_pores["Fluid Volume [µm³]"] / df_pores["Number Of Pores"])) / np.pi)

    #
    return df_pores


def get_fiber_thickness(binary, sampling, sigma_frac=512):
    '''This function takes a binary image of a 3D biopolymer network image and returns the fiber thickness

    Parameters
    ----------
    binary
        The binary image containing the segmentation of fluid (=0) and polymer phase (=1)
    sampling
        The physical size of each voxel in the image
    sigma_frac, optional
        The sigma parameter for the gaussian filter. The smaller the value, the more fibers are detected

    '''

    # calculate the EDM
    small_edge_size = np.min(binary.shape[:2])  # size of the smaller dim of an x-y image
    sigma = small_edge_size / sigma_frac
    edm = get_edm(np.logical_not(binary), sampling)  # edm has microns as values

    # get fiber coordinates
    fiber_coords = peak_local_max(gaussian(edm, sigma=sigma))

    # construct DataFrame with diameter and other infos
    df_fibers = pd.DataFrame(fiber_coords, columns=['x [px]', 'y [px]', 'z [px]'])  # coords are px!
    # we take values from edm, so the result is microns again!
    df_fibers['Diameter [µm]'] = edm[fiber_coords[:, 0], fiber_coords[:, 1], fiber_coords[:, 2]] * 2  # Radius! So x2!

    return df_fibers


def get_fragmented_poresize(binary, sampling, part_size_micron, residual_pore_detection=False, sigma_frac=128):
    data_shape = binary.shape
    intervals_x, intervals_y, intervals_z = get_intervals(data_shape, sampling, part_size_micron)

    df = list()

    for ix, iy, iz in product(intervals_x, intervals_y, intervals_z):
        # each interval (ix...) contains start and stop indices of the respective slice interval
        # cube[ interval[start] : interval[end], ... ) is the correct thing
        cube_slice = binary[ix[0]:ix[1], iy[0]:iy[1], iz[0]:iz[1]]

        # construct DataFrame with measures
        tmp_df = get_pore_sizes(cube_slice, sampling, sigma_frac, residual_pore_detection)

        tmp_df['ix_start'] = ix[0]
        tmp_df['ix_end'] = ix[1]
        tmp_df['iy_start'] = iy[0]
        tmp_df['iy_end'] = iy[1]
        tmp_df['iz_start'] = iz[0]
        tmp_df['iz_end'] = iz[1]
        tmp_df['Slice ID'] = f"{ix[0]}-{ix[1]}_{iy[0]}-{iy[1]}_{iz[0]}-{iz[1]}"

        # save data
        df.append(tmp_df)

    df = pd.concat(df, axis=0, ignore_index=True)

    return df
