BiMaTo
======

|DOI|

.. |DOI| image:: https://zenodo.org/badge/DOI/10.1038/s41598-019-44764-5.svg
   :target: https://doi.org/10.1038/s41598-019-44764-5

|PyPI version shields.io|

.. |PyPI version shields.io| image:: https://img.shields.io/pypi/v/bimato.svg
   :target: https://pypi.python.org/pypi/bimato/

Bio Matrix Topology (BiMaTo) is a library containing all the biopolymer matrix topology analyses published by the Biological Physics Group, (BIP), Peter Debye Institute, University Leipzig, Germany.

Documentation can be found `here <https://bimato.readthedocs.io/>`__.

Exemplary analysis workflow
---------------------------

This is an exemplary workflow to analyze pore sizes of two different collagen scaffolds. The matrices have been fluorescently stained and 3D images were recorded using an LSM.

Usually, we have for example different collagen scaffolds and want to compare their pore-sizes. For this, we would load several images, calculate their pore-sizes and plot them. Below is an exemplary workflow for this:

- load each image in the LIF file
- analyze it
- extract meta-data such as collagen concentration from image name
- concatenate this data to global DataFrame
- plot comparison boxplot

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

.. image:: docs/source/poresize_m.jpeg
  :width: 200
  :align: center
  :alt: boxplot of poresize between two differently concengtrated collagen matrices


.. _article: https://www.nature.com/articles/s41598-019-44764-5