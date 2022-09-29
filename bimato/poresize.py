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
# TODO: write module docstring
docstring of module
"""


__author__ = "Tony Fischer (tku137)"
__copyright__ = "Copyright 2022, Tony Fischer (tku137)"
__license__ = "GPLv3"
__email__ =  "tonyfischer@mailbox.org"
__status__ = "Development"
__version__ = "2022.1"
__credits__ = ["Tony Fischer (tku137)", "Alexander Hayn"]


from itertools import product
from skimage.feature import peak_local_max
from skimage.filters import gaussian
from skimage.morphology import ball
import numpy as np
import pandas as pd
from .core import get_edm


def calc_free_pore_space(df):
    """

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame that came from initial pore detection function

    Returns
    -------
    pandas.Series
        Pandas Series with calculated free pore space

    """
    sum_pore_vol = df.loc[:, "Pore Volume [µm³]"].sum()
    cube_vol = df.loc[:, "Cube Volume [µm³]"].iloc[0]
    return sum_pore_vol / cube_vol


def get_pore_sizes(lif_stack, binary, sigma_frac=128, residual_pore_detection=False):

    # calc fluid phase volume
    voxel_volume = float(lif_stack.info['PhysicalSizeX']) * float(lif_stack.info['PhysicalSizeY']) * float(
        lif_stack.info['PhysicalSizeZ'])
    # print(f"    voxel_volume={voxel_volume}")

    # calculate the EDM
    small_edge_size = np.min(binary.shape[:2])  # size of the smaller dim of an x-y image
    sigma = small_edge_size / sigma_frac
    edm = get_edm(lif_stack, binary)

    # get pore coordinates
    pore_coords = peak_local_max(gaussian(edm, sigma=sigma))

    # construct DataFrame with diameter and other infos
    df_pores = pd.DataFrame(pore_coords, columns=['x [px]', 'y [px]', 'z [px]'])
    df_pores['Diameter [µm]'] = edm[pore_coords[:, 0], pore_coords[:, 1], pore_coords[:, 2]] * 2  # Radius! So x2!
    # TODO: multi-degree
    df_pores['Residual Degree'] = 0
    # calc real pore volume
    pore_stack = np.zeros_like(binary)

    for _, pore in df_pores.iterrows():
        x, y, z, d, _ = pore.astype(np.int)  # need integers to index the array
        r = int((d / 2) / float(lif_stack.info['PhysicalSizeX']))  # convert to pixel-radius

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
        edm_residual = get_edm(lif_stack, residual_fluid)

        # get pore coordinates
        pore_coords_residual = peak_local_max(gaussian(edm_residual, sigma=sigma))

        # construct DataFrame with diameter and other infos
        df_pores_residual = pd.DataFrame(pore_coords_residual, columns=['x [px]', 'y [px]', 'z [px]'])
        df_pores_residual['Diameter [µm]'] = edm_residual[
                                                 pore_coords_residual[:, 0],
                                                 pore_coords_residual[:, 1],
                                                 pore_coords_residual[:, 2]] * 2  # Radius! So x2!
        df_pores_residual['Residual Degree'] = 1

        # calc real pore volume
        for _, pore in df_pores_residual.iterrows():  # TODO: refactor to function
            x, y, z, d, _ = pore.astype(np.int)  # need integers to index the array
            r = int((d / 2) / float(lif_stack.info['PhysicalSizeX']))  # convert to pixel-radius
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
    df_pores['Sample Name'] = lif_stack.image_name
    size_x, size_y, size_z = binary.shape
    df_pores['Size x [px]'] = size_x
    df_pores['Size y [px]'] = size_y
    df_pores['Size z [px]'] = size_z
    df_pores['PhysicalSize x'] = float(lif_stack.info['PhysicalSizeX'])
    df_pores['PhysicalSize y'] = float(lif_stack.info['PhysicalSizeY'])
    df_pores['PhysicalSize z'] = float(lif_stack.info['PhysicalSizeZ'])
    df_pores['PhysicalSizeUnit (x,y,z)'] = str(
        (lif_stack.info['PhysicalSizeXUnit'], lif_stack.info['PhysicalSizeYUnit'], lif_stack.info['PhysicalSizeZUnit']))

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


def get_fiber_thickness(lif_stack, binary, sigma_frac=512):

    # calculate the EDM
    small_edge_size = np.min(binary.shape[:2])  # size of the smaller dim of an x-y image
    sigma = small_edge_size / sigma_frac
    edm = get_edm(lif_stack, np.logical_not(binary))

    # get fiber coordinates
    fiber_coords = peak_local_max(gaussian(edm, sigma=sigma))

    # construct DataFrame with diameter and other infos
    df_fibers = pd.DataFrame(fiber_coords, columns=['x [px]', 'y [px]', 'z [px]'])
    df_fibers['Diameter [µm]'] = edm[fiber_coords[:, 0], fiber_coords[:, 1], fiber_coords[:, 2]] * 2  # Radius! So x2!

    return df_fibers


def get_fragmented_poresize(lif_stack, binary, part_size_micron):
    # voxel_volume = get_voxel_volume(lif_stack)
    intervals_x, intervals_y, intervals_z = get_intervals(lif_stack, part_size_micron)

    df = list()

    for ix, iy, iz in product(intervals_x, intervals_y, intervals_z):
        # each interval (ix...) contains start and stop indices of the respective slice interval
        # cube[ interval[start] : interval[end], ... ) is the correct thing
        cube_slice = binary[ix[0]:ix[1], iy[0]:iy[1], iz[0]:iz[1]]

        # construct DataFrame with measures
        tmp_df = get_pore_sizes(lif_stack, cube_slice)

        tmp_df['ix_start'] = ix[0]
        tmp_df['ix_end'] = ix[1]
        tmp_df['iy_start'] = iy[0]
        tmp_df['iy_end'] = iy[1]
        tmp_df['iz_start'] = iz[0]
        tmp_df['iz_end'] = iz[1]
        tmp_df['Slice ID'] = f"{ix[0]}-{ix[1]}_{iy[0]}-{iy[1]}_{iz[0]}-{iz[1]}"

        tmp_df['Sample name'] = lif_stack.image_name

        # save data
        df.append(tmp_df)

    df = pd.concat(df, axis=0, ignore_index=True)

    return df
