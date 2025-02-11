import numpy as np
from astropy import units as u

from .. import mocpy
from ..abstract_moc import AbstractMOC

__author__ = "Matthieu Baumann, Thomas Boch, Manon Marchand, François-Xavier Pineau"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "matthieu.baumann@astro.unistra.fr, thomas.boch@astro.unistra.fr, manon.marchand@astro.unistra.fr, francois-xavier.pineau@astro.unistra.fr"


class FrequencyMOC(AbstractMOC):
    """Multi-order frequency coverage class. Experimental."""

    # Maximum order of F-MOCs
    MAX_ORDER = np.uint8(59)
    # Upper limit (exclusive) on a F-MOC index = 2^60
    MAX_INDEX_EXCLUSIVE = 1152921504606846976
    # Smallest value, in Hz, a F-MOC can encode
    FREQ_MIN_HZ = 5.048_709_793_414_476e-29
    # Largest value, in Hz, a F-MOC can encode
    FREQ_MAX_HZ = 5.846_006_549_323_611e48

    def __init__(self, store_index):
        """Is a Frequency Coverage (F-MOC).

        Args:
            create_key : Object ensure __init__ is called by super-class/class-methods only
            store_index : index of the `FrequencyMOC` in the rust-side storage
        """
        self.store_index = store_index

    @property
    def max_order(self):
        """Depth/order of the F-MOC.

        Returns
        -------
        `np.uint8`

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> fmoc = FrequencyMOC.from_json({8: [12, 14, 16], 22: [120, 121, 122]})
        >>> print(fmoc.max_order)
        22
        """
        depth = mocpy.get_fmoc_depth(self.store_index)
        return np.uint8(depth)

    @classmethod
    def n_cells(cls, depth):
        """Get the number of cells for a given depth.

        Parameters
        ----------
        depth : int
            The depth. It is comprised between 0 and `~mocpy.fmoc.FrequencyMOC.MAX_ORDER`

        Returns
        -------
        int
            The number of cells at the given order

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> FrequencyMOC.n_cells(0)
        2
        """
        if depth < 0 or depth > cls.MAX_ORDER:
            raise ValueError(
                f"The depth should be comprised between 0 and {cls.MAX_ORDER}, but {depth}"
                " was provided.",
            )
        return mocpy.n_cells_fmoc(depth)

    def to_hz_ranges(self):
        """Return the Hertz ranges this `FrequencyMOC` contains, in Hertz.

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> import astropy.units as u
        >>> fmoc = FrequencyMOC.from_frequency_ranges(10, [1, 0.1, 0.01] * u.Hz, [1.5, 0.5, 0.05] * u.Hz)
        >>> print(fmoc.to_hz_ranges())
        [[0.00976562 0.05078125]
         [0.09375    0.5       ]
         [1.         1.5       ]]
        """
        return np.asarray(
            mocpy.to_freq_ranges(self.store_index) * u.Hz,
            dtype=np.float64,
        )

    @property
    def to_depth59_ranges(self):
        """Return the list of ranges this `FrequencyMOC` contains, at the maximum precision.

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> import astropy.units as u
        >>> fmoc = FrequencyMOC.from_frequency_ranges(59, 1 * u.Hz, 1.4 * u.Hz)
        >>> print(fmoc.to_depth59_ranges)
        [[423338364972826624 425139804823774822]]
        """
        return mocpy.to_ranges(self.store_index)

    def degrade_to_order(self, new_order):
        """
        Degrade the `FrequencyMOC` instance to a new, less precise, `FrequencyMOC`.

        The maximum depth (i.e. the depth of the smallest Time cells that can be found in the F-MOC) of the
        degraded F-MOC is set to ``new_order``.

        Parameters
        ----------
        new_order : `int`
            Maximum depth of the output degraded F-MOC.

        Returns
        -------
        `FrequencyMOC`
            The degraded F-MOC.

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> import astropy.units as u
        >>> fmoc = FrequencyMOC.from_frequencies(40, 1 * u.Hz)
        >>> fmoc
        40/807453851648
        >>> fmoc.degrade_to_order(10)
        10/752
        """
        index = mocpy.degrade(self.store_index, new_order)
        return FrequencyMOC(index)

    @classmethod
    def new_empty(cls, max_depth):
        """
        Create a new empty `FrequencyMOC` of given depth.

        Parameters
        ----------
        max_depth : `int`
                The resolution of the FrequencyMOC


        Returns
        -------
        `FrequencyMOC`

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> FrequencyMOC.new_empty(5)
        5/
        """
        index = mocpy.new_empty_fmoc(np.uint8(max_depth))
        return cls(index)

    @classmethod
    def from_depth59_ranges(cls, order, ranges):
        """
        Create a `FrequencyMOC` from a set of `FrequencyMOC` ranges at order 59.

        Parameters
        ----------
        order : `int`, The resolution of the `FrequencyMOC`
        ranges : `numpy.ndarray` or `list`
                 a N x 2 numpy array or list representing the set of depth 61 ranges.

        Returns
        -------
        `~mocpy.fmoc.FrequencyMOC`
            The F-MOC

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> FrequencyMOC.from_depth59_ranges(40, [[0, 10000000]])
        36/0
        38/4
        40/
        """
        ranges = np.zeros((0, 2), dtype=np.uint64) if ranges is None else ranges

        ranges = np.array(ranges)

        if ranges.shape[1] != 2:
            raise ValueError(
                f"Expected a N x 2 `list` or `numpy.ndarray` but second dimension is {ranges.shape[1]}",
            )

        if ranges.dtype is not np.uint64:
            ranges = ranges.astype(np.uint64)

        index = mocpy.from_fmoc_ranges_array2(np.uint8(order), ranges)
        return cls(index)

    @classmethod
    def from_frequencies(cls, order, frequencies):
        """
        Create a `FrequencyMOC` from a `astropy.units.Quantity` that are internally converted in Hertz.

        Parameters
        ----------
        order : `int`
            The resolution of the FrequencyMOC: see `relative_precision_to_order`
        frequencies : `astropy.units.Quantity`
            Quantity converted internally in Hertz

        Returns
        -------
        `FrequencyMOC`

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> import astropy.units as u
        >>> FrequencyMOC.from_frequencies(42, [1e-6, 1e-3, 1] * u.Hz)
        42/2544289697882 2887042656632 3229815406592
        """
        frequencies = frequencies.to(u.Hz)
        frequencies = np.atleast_1d(frequencies)

        store_index = mocpy.from_freq_values(order, frequencies)
        return cls(store_index)

    @classmethod
    def from_frequency_ranges(cls, order, min_freq, max_freq):
        """
        Create a `FrequencyMOC` from a range defined by two `astropy.units.Quantity` converted in Hertz.

        Parameters
        ----------
        order : `int`
            The resolution of the `FrequencyMOC`: see `relative_precision_to_order`
        min_freq : `astropy.units.Quantity`
            astropy quantity converted in Hertz and defining the left part of the intervals
        max_freq : `astropy.units.Quantity`
            astropy quantity converted in Hertz and defining the right part of the intervals

        Returns
        -------
        `FrequencyMOC`

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> import astropy.units as u
        >>> FrequencyMOC.from_frequency_ranges(10, [10, 40] * u.Hz, [20, 60] * u.Hz)
        8/195
        9/389 392 397-398
        10/798
        """
        min_freq = min_freq.to(u.Hz)
        min_freq = np.atleast_1d(min_freq)

        max_freq = max_freq.to(u.Hz)
        max_freq = np.atleast_1d(max_freq)

        if min_freq.shape != max_freq.shape:
            raise ValueError(
                f"Mismatch between min_freq and max_freq of shapes {min_freq.shape} and {max_freq.shape}",
            )

        store_index = mocpy.from_freq_ranges(
            order,
            min_freq,
            max_freq,
        )
        return cls(store_index)

    @property
    def min_freq(self):
        """
        Get the `~astropy.units.Quantity` frequency of the F-MOC smallest frequency.

        Returns
        -------
        min_freq : `astropy.units.Quantity`
            frequency of the first observation

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> import astropy.units as u
        >>> fmoc = FrequencyMOC.from_frequencies(10, [1, 10] * u.Hz)
        >>> print(fmoc.min_freq)
        1.0 Hz
        """
        return mocpy.first_fmoc_hz(self.store_index) * u.Hz

    @property
    def max_freq(self):
        """
        Get the `astropy.units.Quantity` largest frequency of the F-MOC.

        Returns
        -------
        max_freq : `astropy.units.Quantity`
            frequency of the last observation

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> import astropy.units as u
        >>> fmoc = FrequencyMOC.from_frequencies(10, [1, 10] * u.Hz)
        >>> # this order is pretty low, thus the returned max frequency
        >>> # corresponds to the high limit of the cell containing 10 Hz
        >>> # at order 10
        >>> print(fmoc.max_freq)
        11.0 Hz

        """
        return mocpy.last_fmoc_hz(self.store_index) * u.Hz

    def contains(self, frequencies, keep_inside=True):
        """
        Test is a frequency -- or list of frequencies -- is comprised in this `FrequencyMOC`.

        Parameters
        ----------
        frequencies : `astropy.units.Quantity`
            astropy quantity (converted into Hz) to check whether they are contained in the F-MOC or not.
        keep_inside : `bool`, optional
            True by default. If so the filtered table contains only observations that are located the MOC.
            If ``keep_inside`` is False, the filtered table contains all observations lying outside the MOC.

        Returns
        -------
        array : `numpy.ndarray`
            A mask boolean array

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> import astropy.units as u
        >>> # Let's create a FMOC at order 10 for frequencies comprised between
        >>> # 1Hz and 10Hz.
        >>> fmoc = FrequencyMOC.from_frequency_ranges(10, 1 * u.Hz, 10 * u.Hz)
        >>> # We can now test wether fmoc contains a list of frequencies between
        >>> # 1Hz and 15Hz.
        >>> fmoc.contains(range(1, 15, 1) * u.Hz)
        array([ True,  True,  True,  True,  True,  True,  True,  True,  True,
               False, False, False, False, False])
        """
        freq = frequencies.to(u.Hz)
        freq = np.atleast_1d(freq)

        mask = mocpy.filter_freq(self.store_index, freq)

        if keep_inside:
            return mask
        return ~mask

    @staticmethod
    def order_to_relative_precision(order):
        r"""
        Convert a `FrequencyMOC` order to a **relative** precision range.

        Parameters
        ----------
        order : `int`
            order to convert

        Returns
        -------
        `numpy.ndarray`
            lower and upper relative precisions

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> FrequencyMOC.order_to_relative_precision(10)
        array([0.0625, 0.125 ])
        >>> FrequencyMOC.order_to_relative_precision(20)
        array([6.10351562e-05, 1.22070312e-04])

        Notes
        -----
        In FMOCs, the precision of a cell depends on its position
        along the electromagnetic axis. The array returned by
        `order_to_relative_precision` corresponds to the lower and upper
        **relative** precisions.
        These must be multiplied by the value of the observed frequency
        to obtain the **absolute** upper and lower precisions.

        In the code example, we see that at order 10:

        .. math:: F = 10_{- 0.6}^{+ 1}~\mathrm{Hz}

        .. math:: F = 1e3_{- 6e1}^{+ 1e2}~\mathrm{Hz}

        At order 20 these precisions become:

        .. math:: F = 10_{- 6e-4}^{+ 1e-3}~\mathrm{Hz}

        .. math:: F = 1e3_{- 6e-2}^{+ 1e-1}~\mathrm{Hz}

        """
        if order > FrequencyMOC.MAX_ORDER:
            raise ValueError(f"FMOCs have a maximum order of {FrequencyMOC.MAX_ORDER}")

        return np.power([2.0, 2.0], [6 - order, 7 - order])

    @staticmethod
    def relative_precision_to_order(frequency_precision):
        """
        Convert a relative frequency precision into a `FrequencyMOC` order.

        Parameters
        ----------
        frequency_precision : `float`
            precision to be converted in FrequencyMOC order

        Returns
        -------
        order : `int`
            The less precise order fulfilling the given ``frequency_precision``.

        """
        order = int(7 - np.log(frequency_precision) / np.log(2))
        if order > FrequencyMOC.MAX_ORDER:
            order = FrequencyMOC.MAX_ORDER
        elif order < 0:
            order = 0
        return np.uint8(order)

    @classmethod
    def load(cls, path, format="fits"):  # noqa: A002
        """
        Load the `FrequencyMOC` from a file.

        Format can be 'fits', 'ascii', or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        path : `str` or `pathlib.Path`
            The path to the file to load the F-MOC from.
        format : `str`, optional {'ascii', 'fits', 'json'}
            The format from which the F-MOC is loaded.
            Possible formats are "fits", "ascii" or "json".
            By default, ``format`` is set to "fits".
        """
        path = str(path)
        if format == "fits":
            index = mocpy.frequency_moc_from_fits_file(path)
            return cls(index)
        if format == "ascii":
            index = mocpy.frequency_moc_from_ascii_file(path)
            return cls(index)
        if format == "json":
            index = mocpy.frequency_moc_from_json_file(path)
            return cls(index)
        formats = ("fits", "ascii", "json")
        raise ValueError(f"format should be one of {formats}")

    @classmethod
    def _from_fits_raw_bytes(cls, raw_bytes):
        """Load FrequencyMOC from raw bytes of a FITS file."""
        index = mocpy.frequency_moc_from_fits_raw_bytes(raw_bytes)
        return cls(index)

    @classmethod
    def from_string(cls, value, format="ascii"):  # noqa: A002
        """
        Deserialize the `FrequencyMOC` from the given string.

        Format can be 'ascii' or 'json', though the json format is not officially supported by the IVOA.

        WARNING: the serialization must be strict, i.e. **must not** contain overlapping elements

        Parameters
        ----------
        format : `str`, optional
            The format in which the F-MOC will be serialized before being saved.
            Possible formats are "ascii" or "json".
            By default, ``format`` is set to "ascii".

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> FrequencyMOC.from_string("4/4")
        4/4
        """
        if format == "ascii":
            index = mocpy.frequency_moc_from_ascii_str(value)
            return cls(index)
        if format == "json":
            index = mocpy.frequency_moc_from_json_str(value)
            return cls(index)
        formats = ("ascii", "json")
        raise ValueError(f"format should be one of {formats}")

    def plot_frequencies(self, ax, color="blue", frequency_unit="Hz"):
        """Plot a frequency moc.

        This method applies a `matplotlib.collections.PatchCollection`
        to an existing `matplotlib.axes._axes.Axes` object.

        Parameters
        ----------
        ax : `matplotlib.axes._axes.Axes`
        color : `str`, default 'blue'
               any format supported by matplotlib for colors, see
               `matplotlib.colors`
        length_unit : `str` or `astropy.units.core.Unit`, optional
                     any string or astropy.unit of physical type 'frequency', see
                     `astropy.units.physical.get_physical_type`
                     Defaults to Hertz 'Hz'

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> import matplotlib.pyplot as plt
        >>> import astropy.units as u
        >>> fmoc = FrequencyMOC.from_frequencies(10, [1, 0.1, 0.01, 0.001] * u.Hz)
        >>> fig, ax = plt.subplots(figsize=(15, 1))
        >>> fmoc.plot_frequencies(ax, color="pink", frequency_unit="1 / ks")
        """
        if u.get_physical_type(u.Unit(frequency_unit)) != "frequency":
            raise TypeError(
                f"frequency_unit is of type '{u.get_physical_type(u.Unit(frequency_unit))}'"
                " instead of 'frequency', see astropy.units for more information",
            )

        from matplotlib.collections import PatchCollection
        from matplotlib.patches import Rectangle

        min_freq = self.min_freq.to(frequency_unit).value
        max_freq = self.max_freq.to(frequency_unit).value

        patches = []
        for freq_range in self.to_hz_ranges():
            freq0 = (freq_range[0] * u.Hz).to(frequency_unit).value
            freq1 = (freq_range[1] * u.Hz).to(frequency_unit).value
            patches += [Rectangle((freq0, 0), freq1 - freq0, 1, color=color)]
        patches_collection = PatchCollection(patches, match_original=True)
        ax.add_collection(patches_collection)
        ax.tick_params(left=False)
        ax.set(
            yticklabels=[],
            xlim=(min_freq, max_freq),
            xlabel=f"Frequency ({frequency_unit})",
            xscale="log",
        )

    def plot_wavelengths(self, ax, color="blue", length_unit="m"):
        """Plot a `FrequencyMOC` with a conversion to wavelengths.

        This method applies a `matplotlib.collections.PatchCollection`
        to an existing `matplotlib.axes._axes.Axes` object.

        Parameters
        ----------
        ax : `matplotlib.axes._axes.Axes`
        color : `str`, default 'blue'
               any format supported by matplotlib for colors, see
               `matplotlib.colors`.
        length_unit : `str` or `astropy.units.core.Unit`, default 'm'
                     any string or astropy.unit of physical type 'length', see
                     `astropy.units.get_physical_type`
                     Defaults to meters 'm'

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> import matplotlib.pyplot as plt
        >>> import astropy.units as u
        >>> fmoc = FrequencyMOC.from_frequencies(10, [1, 0.1, 0.01, 0.001] * u.Hz)
        >>> fig, ax = plt.subplots(figsize=(15, 1))
        >>> fmoc.plot_wavelengths(ax, color="lightblue", length_unit=u.nm)
        """
        # Tests the physical type of `length_unit`
        if u.get_physical_type(u.Unit(length_unit)) != "length":
            raise TypeError(
                f"length_unit is of type '{u.get_physical_type(u.Unit(length_unit))}'"
                " instead of 'length', see astropy.units for more information",
            )

        from matplotlib.collections import PatchCollection
        from matplotlib.patches import Rectangle

        # get default bonds
        min_lambda = self.max_freq.to(length_unit, equivalencies=u.spectral()).value
        max_lambda = self.min_freq.to(length_unit, equivalencies=u.spectral()).value

        # fetches the patches corresponding to each frequency range
        patches = []
        for freq_range in self.to_hz_ranges():
            # converts to wavelengths
            lambda0 = (
                (freq_range[1] * u.Hz).to(length_unit, equivalencies=u.spectral()).value
            )
            lambda1 = (
                (freq_range[0] * u.Hz).to(length_unit, equivalencies=u.spectral()).value
            )
            patches += [Rectangle((lambda0, 0), lambda1 - lambda0, 1, color=color)]
        patches_collection = PatchCollection(patches, match_original=True)

        # default behaviour, can be overwritten by the user with a call on the `ax` object
        ax.add_collection(patches_collection)
        ax.tick_params(left=False)
        ax.set(
            yticklabels=[],
            xlim=(min_lambda, max_lambda),
            xlabel=f"Wavelength ({length_unit})",
            xscale="log",
        )
