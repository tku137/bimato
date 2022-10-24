BiMaTo
======

|DOI1| |DOI2| |GPLv3 license| |PyPI version shields.io| |Documentation Status|

.. |DOI1| image:: https://zenodo.org/badge/DOI/10.1038/s41598-019-44764-5.svg
   :target: https://doi.org/10.1038/s41598-019-44764-5

.. |DOI2| image:: https://zenodo.org/badge/DOI/10.3389/fcell.2020.593879.svg
   :target: https://doi.org/10.3389/fcell.2020.593879

.. |PyPI version shields.io| image:: https://img.shields.io/pypi/v/bimato.svg
   :target: https://pypi.python.org/pypi/bimato/

.. |GPLv3 license| image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: http://perso.crans.org/besson/LICENSE.html

.. |Documentation Status| image:: https://readthedocs.org/projects/bimato/badge/?version=latest
   :target: http://bimato.readthedocs.io/?badge=latest

Bio Matrix Topology (BiMaTo) is a library containing all the biopolymer matrix topology analyses published by the Biological Physics Group, (BIP), Peter Debye Institute, University Leipzig, Germany.

Documentation can be found `here <https://bimato.readthedocs.io/>`__.

How to install
--------------

**bimato** uses Python3.8 and up. Installation is trivial::

    pip install bimato

Exemplary analysis workflow
---------------------------

This is an exemplary workflow to analyze pore sizes of two different collagen scaffolds. The matrices have been fluorescently stained and 3D images were recorded using an LSM.

Usually, we have for example different collagen scaffolds and want to compare structural parameters. For this, we would load several images, calculate their structural parameters and plot them. Below is an exemplary workflow for this:

- load each image in the LIF file
- analyze it
- extract meta-data such as collagen concentration from image name
- concatenate this data to global DataFrame
- plot comparison boxplot

Pore-Size
^^^^^^^^^

We load a lif file with multiple samples per collagen concentration and analyze these in a loop:

..  code-block:: python

    import pandas as pd
    from readlif.reader import LifFile
    import seaborn as sns
    import bimato

    lif_file = LifFile("/path/to/sample.lif")

    df_poresize = list()
    for lif_image in lif_file.get_iter_image():

        data = bimato.utils.read_lif_image(lif_image)
        bw = bimato.get_binary(data)

        sampling = {
            'x': 1/lif_image.info["scale"][0],
            'y': 1/lif_image.info["scale"][1],
            'z': 1/lif_image.info["scale"][2]
        }

        df_tmp = bimato.get_pore_sizes(binary=bw, sampling=sampling)

        df_tmp['Concentration [g/l]'] = lif_image.name
        df_poresize.append(df_tmp)

    df_poresize = pd.concat(df_poresize)

    g = sns.catplot(
        data=df_poresize,
        kind='box',
        x='Concentration [g/l]',
        y='Diameter [µm],
    )
    g.set_ylabels("Pore-size [µm]")

Resulting in the following plot:

.. image:: https://github.com/tku137/bimato/raw/main/docs/source/poresize_m.jpeg
  :width: 200
  :align: center
  :alt: boxplot of poresize between two differently concentrated collagen matrices

Inhomogeneity
^^^^^^^^^^^^^

We load a lif file with multiple samples per collagen concentration and analyze these in a loop:

..  code-block:: python

    import pandas as pd
    from readlif.reader import LifFile
    import seaborn as sns
    import bimato

    lif_file = LifFile("/path/to/sample.lif")

    df_inhomogeneity = list()
    for lif_image in lif_file.get_iter_image():

        data = bimato.utils.read_lif_image(lif_image)
        bw = bimato.get_binary(data)

        sampling = {
            'x': 1/lif_image.info["scale"][0],
            'y': 1/lif_image.info["scale"][1],
            'z': 1/lif_image.info["scale"][2]
        }

        df_tmp = bimato.poresize.get_fragmented_poresizes(binary=bw, sampling=sampling, part_size_micron=30)
        df_tmp['Inhomogeneity'] = bimato.poresize.calc_inhomogeneity(df_tmp)

        df_tmp['Concentration [g/l]'] = lif_image.name
        df_inhomogeneity.append(df_tmp)

    df_inhomogeneity = pd.concat(df_inhomogeneity)

    g = sns.catplot(
        data=df_poresize,
        kind='box',
        x='Concentration [g/l]',
        y='Inhomogeneity,
    )

Resulting in the following plot:

.. image:: https://github.com/tku137/bimato/raw/main/docs/source/inhomogeneity_m.jpeg
  :width: 200
  :align: center
  :alt: boxplot of inhomogeneity between two differently concentrated collagen matrices

How to cite
-----------

Fischer T, Hayn A, Mierke CT (2019) Fast and reliable advanced two-step pore-size analysis of biomimetic 3D extracellular matrix scaffolds. Scientific Reports 9:8352. https://doi.org/10.1038/s41598-019-44764-5
Hayn A, Fischer T, Mierke CT (2020) Inhomogeneities in 3D Collagen Matrices Impact Matrix Mechanics and Cancer Cell Migration. Front Cell Dev Biol 8:593879. https://doi.org/10.3389/fcell.2020.593879

