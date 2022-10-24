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
__version__ = "2022.2"
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


def get_fragmented_poresizes(binary, sampling, part_size_micron, residual_pore_detection=False):
    """Calculates fragmented pore-sizes of a sample, given the desired part size.

    First, based on the given part size, the possible number of parts and their respective coordinate intervals are
    determined. Using the previous coordinate intervals, the binary segmentation of the original sample is split into
    multiple smaller samples and pore size calculations and statistics are calculated for each individual part. This
    is the basis for further inhomogeneity analyses, as described at
    https://www.frontiersin.org/articles/10.3389/fcell.2020.593879/.

    The resulting :class:`pandas.DataFrame` contains mainly the same columns as desccribed in :func:`bimato.poresize.get_pore_sizes`.
    However, it is extended by several columns containing the start and end coordinates of each part inside the original
    data, a unique slice-ID which can be used for example in :func:`pandas.DataFrame.groupby`, and several other data
    columns that are neccessary to calculate the inhomogeneity using :func:`bimato.poresize.calc_inhomogeneity`.

    Parameters
    ----------
    binary
        The binary image containing the segmentation of fluid (=0) and polymer phase (=1)
    sampling
        The physical size of each voxel in the image
    part_size_micron
        Desired size of a single part in microns
    residual_pore_detection
        If True, the algorithm will try to detect even smaller and obscure pores

    Returns
    -------
        DataFrame with dedicated pore sizes and statistics for individual parts of the original sample
    """
    data_shape = binary.shape
    intervals_x, intervals_y, intervals_z = get_intervals(data_shape, sampling, part_size_micron)

    df = list()
    for ix, iy, iz in product(intervals_x, intervals_y, intervals_z):
        # each interval (ix...) contains start and stop indices of the respective slice interval

        # query data with coordinates in specified intervals, represents slice of image cube
        cube_slice = binary[ix[0]:ix[1], iy[0]:iy[1], iz[0]:iz[1]]

        # construct DataFrame with measures
        tmp_df = get_pore_sizes(cube_slice, sampling, sigma_frac=128, residual_pore_detection=residual_pore_detection)

        # store metadata of slice
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

    df = get_part_amount(df)

    return df


def get_part_amount(df):
    """Integral part of :func:`bimato.poresize.get_fragmented_poresizes`. Calculates specific parameters about the
    collagen content and pore sizes of each fragment.


    Parameters
    ----------
    df
        DataFrame with pore sizes

    Returns
    -------
        Extended DataFrame with specific parameters fundamental to the inhomogeneity calculation
    """
    tmp_df = df.copy()
    tmp_df.loc[:, "Relative Number Of Pores"] = tmp_df.loc[:, "Number Of Pores"] / tmp_df["Number Of Pores"].sum()
    tmp_df.loc[:, "Relative Real Pore Volume"] = tmp_df.loc[:, "Real Pore Volume [µm³]"] / tmp_df["Real Pore Volume [µm³]"].sum()
    tmp_df.loc[:, "Relative Collagen Volume"] = tmp_df.loc[:, "Collagen Volume [µm³]"] / tmp_df["Collagen Volume [µm³]"].sum()
    tmp_df.loc[:, "Relative Diameter"] = tmp_df.loc[:, "Diameter [µm]"] / tmp_df["Diameter [µm]"].sum()
    return tmp_df


def calc_inhomogeneity(df):
    """Calculates the inhomogeneity parameter based on sophisticated statistics and parameters, as described at
    https://www.frontiersin.org/articles/10.3389/fcell.2020.593879/. Takes a :class:`pandas.DataFrame` that has been
    extended by :func:`bimato.poresize.get_fragmented_poresizes` and returns the inhomogeneity parameter of that
    DataFrame.

    Note
    ----
    This function is especially useful when combined with :func:`pandas.DataFrame.groupby` and
    :func:`pandas.core.groupby.GroupBy.apply` operations. See :doc:`user-guide` for more information.

    Parameters
    ----------
    df
        pandas.DataFrame with specific parameters calculated by :func:`bimato.poresize.get_fragmented_poresizes`

    Returns
    -------
        Inhomogeneity parameter
    """
    std_df = df.loc[:, ['Relative Number Of Pores', 'Relative Collagen Volume', 'Relative Diameter']].std()
    inhomogeneity = np.sqrt(
        std_df['Relative Number Of Pores'] ** 2 +
        std_df['Relative Collagen Volume'] ** 2 +
        std_df['Relative Diameter'] ** 2
    )
    return inhomogeneity
