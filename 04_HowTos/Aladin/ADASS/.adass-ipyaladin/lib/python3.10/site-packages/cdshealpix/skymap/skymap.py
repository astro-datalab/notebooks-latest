"""Manipulation of skymaps.

SkyMaps are described in <Data formats for gamma-ray astronomy
 https://gamma-astro-data-formats.readthedocs.io/en/latest/skymaps/healpix/index.html>_
This sub-module supports skymaps in the nested scheme, and in the implicit format.
The coordinates system should be 'CEL'.
"""
from .. import cdshealpix

from pathlib import Path
from typing import Union

try:
    import matplotlib.pyplot as plt
except ImportError:
    _matplotlib_missing = True
import numpy as np


class Skymap:
    """A Skymap, containing values to associate to healpix cells."""

    def __init__(self, values):
        self.values = values

    @property
    def depth(self):
        """The depth of the skymap.

        Avoids the costly log calculation.

        Returns
        -------
        `int`
            The depth of the skymap.

        Examples
        --------
        >>> from cdshealpix.skymap import Skymap
        >>> map = Skymap.from_array([0]*12)
        >>> map.depth
        0
        """
        return cdshealpix.depth_skymap(self.values)

    @classmethod
    def from_fits(cls, path: Union[str, Path]):
        """Read a skymap in the nested schema from a FITS file.

        This reader supports files which are:

        - all sky maps
        - in the nested scheme
        - and the implicit format

        Parameters
        ----------
        path : str, `pathlib.Path`
            The file's path.

        Returns
        -------
        `Skymap`
            A skymap. Its values are in a numpy array which data type in inferred from
            the FITS header.
        """
        return cls(cdshealpix.read_skymap(str(path)))

    @classmethod
    def from_array(cls, values):
        """Create a skymap from an array.

        Parameters
        ----------
        values : `numpy.array`
            An array-like object. It should be one-dimensional, and its length should be
            the number of cells in a HEALPix order.
            It should be in the nested ordering (not tested).

        Returns
        -------
        `Skymap`
            A skymap object.

        Examples
        --------
        >>> from cdshealpix.skymap import Skymap
        >>> import numpy as np
        >>> skymap =Skymap.from_array(np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], dtype=np.uint8))
        """
        # only makes a copy if it was not C-contiguous in the first place
        values = np.ascontiguousarray(values)
        if values.ndim != 1:
            raise ValueError(
                "Skymap values should be one-dimensional. Got an array of "
                f"shape {values.shape}."
            )
        n = int(len(values) / 12)
        # test if it is a power of two (1000 & 0111 = 0000)
        if n & (n - 1) != 0 or n == 0:
            raise ValueError(
                "The length of values should be a valid number of cells in "
                "a given HEALPix order, i.e something like 12, 48, 192... "
                f"Got '{len(values)}'."
            )

        if values.dtype not in (
            np.float64,
            np.float32,
            np.int64,
            np.int32,
            np.int16,
            np.uint8,
            float,
            int,
        ):
            raise ValueError(
                f"The accepted types are f64, f32, i64, i32, u8. Got '{values.dtype}'."
            )
        return cls(values)

    def to_fits(self, path):
        """Write a Skymap in a fits file.

        Parameters
        ----------
        path : `str`, `pathlib.Path`
            The file's path.
        """
        cdshealpix.write_skymap(self.values, str(path))

    def quick_plot(self, *, size=256, convert_to_gal=True, path=None):
        """Preview a skymap in the Mollweide projection.

        Parameters
        ----------
        size : `int`, optional
            The size of the plot in the y-axis in pixels.
            It fixes the resolution of the image. By default 256
        convert_to_gal : `bool`, optional
            Should the image be converted into a galactic frame? by default True
        path : `str` or `pathlib.Path`, optional
            If different from none, the image will not only be displayed, but also saved
            at the given location. By default None
        """
        if _matplotlib_missing:
            raise ModuleNotFoundError(
                "matplotlib is mandatory to use 'quick_plot'. "
                "See https://matplotlib.org/ for installation "
                "instructions."
            )
        img = cdshealpix.pixels_skymap(self.values, size, convert_to_gal)
        fig = plt.imshow(img)
        plt.axis("off")
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        if path:
            plt.savefig(path, bbox_inches="tight", pad_inches=0, transparent=True)
        plt.show()
