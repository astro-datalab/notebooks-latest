import warnings

import numpy as np
from astropy.time import Time, TimeDelta

from .. import mocpy
from ..abstract_moc import AbstractMOC

__author__ = "Matthieu Baumann, Thomas Boch, Manon Marchand, François-Xavier Pineau"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "matthieu.baumann@astro.unistra.fr, thomas.boch@astro.unistra.fr, manon.marchand@astro.unistra.fr, francois-xavier.pineau@astro.unistra.fr"


DAY_MICRO_SEC = 86400000000.0


def times_to_microseconds(times):
    """
    Convert a `astropy.time.Time` into an array of integer microseconds since JD=0, keeping the microsecond resolution required for `~mocpy.tmoc.TimeMOC`.

    Parameters
    ----------
    times : `astropy.time.Time`
    Astropy observation times

    Returns
    -------
    times_microseconds : `np.array`
    """
    times_jd = np.asarray(times.jd, dtype=np.uint64)
    times_us = np.asarray(
        (times - Time(times_jd, format="jd", scale="tdb")).jd * DAY_MICRO_SEC,
        dtype=np.uint64,
    )

    return times_jd * np.uint64(DAY_MICRO_SEC) + times_us


def microseconds_to_times(times_microseconds):
    """
    Convert an array of integer microseconds since JD=0, to an array of `astropy.time.Time`.

    Parameters
    ----------
    times_microseconds : `np.array`

    Returns
    -------
    times : `astropy.time.Time`
    """
    jd1 = np.asarray(times_microseconds // DAY_MICRO_SEC, dtype=np.float64)
    jd2 = np.asarray(
        (times_microseconds - jd1 * DAY_MICRO_SEC) / DAY_MICRO_SEC,
        dtype=np.float64,
    )

    return Time(val=jd1, val2=jd2, format="jd", scale="tdb")


class TimeMOC(AbstractMOC):
    """Multi-order time coverage class. Experimental."""

    # Maximum order of TimeMOCs
    # (do not remove since it may be used externally).
    MAX_ORDER = np.uint8(61)
    # Number of microseconds in a day
    DAY_MICRO_SEC = 86400000000.0
    # Default observation time : 30 min
    DEFAULT_OBSERVATION_TIME = TimeDelta(30 * 60, format="sec", scale="tdb")

    def __init__(self, store_index):
        """Is a Time Coverage (T-MOC).

        Args:
            store_index: index of the S-MOC in the rust-side storage
        """
        self.store_index = store_index

    @property
    def max_order(self):
        """Depth/order of the T-MOC."""
        depth = mocpy.get_tmoc_depth(self.store_index)
        return np.uint8(depth)

    @classmethod
    def n_cells(cls, depth):
        """Get the number of cells for a given depth.

        Parameters
        ----------
        depth : int
            The depth. It is comprised between 0 and `~mocpy.tmoc.TimeMOC.MAX_ORDER`

        Returns
        -------
        int
            The number of cells at the given order

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> TimeMOC.n_cells(0)
        2
        """
        if depth < 0 or depth > cls.MAX_ORDER:
            raise ValueError(
                f"The depth should be comprised between 0 and {cls.MAX_ORDER}, but {depth}"
                " was provided.",
            )
        return mocpy.n_cells_tmoc(depth)

    def to_time_ranges(self):
        """Return the time ranges this TimeMOC contains."""
        return microseconds_to_times(mocpy.to_ranges(self.store_index))

    @property
    def to_depth61_ranges(self):
        """Return the list of ranges this TimeMOC contains, in microsec since JD=0."""
        return mocpy.to_ranges(self.store_index)

    def degrade_to_order(self, new_order):
        """
        Degrade the MOC instance to a new, less precise, MOC.

        The maximum depth (i.e. the depth of the smallest Time cells that can be found in the MOC) of the
        degraded MOC is set to ``new_order``.

        Parameters
        ----------
        new_order : int
            Maximum depth of the output degraded MOC.

        Returns
        -------
        moc : `~mocpy.tmoc.TimeMOC`
            The degraded MOC.
        """
        index = mocpy.degrade(self.store_index, new_order)
        return TimeMOC(index)

    @classmethod
    def new_empty(cls, max_depth):
        """
        Create a new empty TimeMOC of given depth.

        Parameters
        ----------
        max_depth : int, The resolution of the TimeMOC


        Returns
        -------
        moc : `~mocpy.tmoc.TimeMOC`
            The MOC
        """
        index = mocpy.new_empty_tmoc(np.uint8(max_depth))
        return cls(index)

    @classmethod
    def from_depth61_ranges(cls, max_depth, ranges):
        """
        Create a TimeMOC from a set of Time ranges at order 61 (i.e. ranges of microseconds since JD=0).

        Parameters
        ----------
        max_depth : int, The resolution of the TimeMOC
        ranges: `~numpy.ndarray`
                 a N x 2 numpy array representing the set of depth 61 ranges.

        Returns
        -------
        moc : `~mocpy.tmoc.TimeMOC`
            The MOC
        """
        ranges = np.zeros((0, 2), dtype=np.uint64) if ranges is None else ranges

        if ranges.shape[1] != 2:
            raise ValueError(
                f"Expected a N x 2 numpy ndarray but second dimension is {ranges.shape[1]}",
            )

        if ranges.dtype is not np.uint64:
            ranges = ranges.astype(np.uint64)

        index = mocpy.from_time_ranges_array2(np.uint8(max_depth), ranges)
        return cls(index)

    @classmethod
    def from_times(cls, times, delta_t=DEFAULT_OBSERVATION_TIME):
        """
        Create a TimeMOC from a `astropy.time.Time`.

        Parameters
        ----------
        times : `astropy.time.Time`
            Astropy observation times
        delta_t : `astropy.time.TimeDelta`, optional
            The duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMOC order to represent the observations (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``).

        Returns
        -------
        time_moc : `~mocpy.tmoc.TimeMOC`
        """
        times = times_to_microseconds(times)
        times = np.atleast_1d(times)

        depth = TimeMOC.time_resolution_to_order(delta_t)
        store_index = mocpy.from_time_in_microsec_since_jd_origin(depth, times)
        return cls(store_index)

    @classmethod
    def from_time_ranges(cls, min_times, max_times, delta_t=DEFAULT_OBSERVATION_TIME):
        """
        Create a TimeMOC from a range defined by two `astropy.time.Time`.

        Parameters
        ----------
        min_times : `astropy.time.Time`
            astropy times defining the left part of the intervals
        max_times : `astropy.time.Time`
            astropy times defining the right part of the intervals
        delta_t : `astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMOC order to represent the observations (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``).

        Returns
        -------
        time_moc : `~mocpy.tmoc.TimeMOC`
        """
        # degrade the TimeMoc to the order computed from ``delta_t``
        depth = TimeMOC.time_resolution_to_order(delta_t)

        min_times = times_to_microseconds(min_times)
        min_times = np.atleast_1d(min_times)

        max_times = times_to_microseconds(max_times)
        max_times = np.atleast_1d(max_times)

        if min_times.shape != max_times.shape:
            raise ValueError(
                f"Mismatch between min_times and max_times of shapes {min_times.shape} and {max_times.shape}",
            )

        store_index = mocpy.from_time_ranges_in_microsec_since_jd_origin(
            depth,
            min_times,
            max_times,
        )
        return cls(store_index)

    @classmethod
    def from_time_ranges_approx(
        cls,
        min_times,
        max_times,
        delta_t=DEFAULT_OBSERVATION_TIME,
    ):
        """
        Create a TimeMOC from a range defined by two `astropy.time.Time`.

        Uses the following approximation: simple take the JD time and multiply by the number of microseconds in a day.

        Parameters
        ----------
        min_times : `astropy.time.Time`
            astropy times defining the left part of the intervals
        max_times : `astropy.time.Time`
            astropy times defining the right part of the intervals
        delta_t : `astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMOC order to represent the observations (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``).

        Returns
        -------
        time_moc : `~mocpy.tmoc.TimeMOC`
        """
        # degrade the TimeMoc to the order computed from ``delta_t``
        depth = TimeMOC.time_resolution_to_order(delta_t)

        min_times = np.asarray(min_times.jd)
        min_times = np.atleast_1d(min_times)

        max_times = np.asarray(max_times.jd)
        max_times = np.atleast_1d(max_times)
        if min_times.shape != max_times.shape:
            raise ValueError(
                f"Mismatch between min_times and max_times of shapes {min_times.shape} and {max_times.shape}",
            )

        store_index = mocpy.from_time_ranges(depth, min_times, max_times)
        return cls(store_index)

    @classmethod
    def from_stmoc_space_fold(cls, smoc, stmoc):
        """
        Build a new T-MOC from the fold operation of the given ST-MOC by the given S-MOC.

        Parameters
        ----------
        smoc : `~mocpy.moc.Moc`
        stmoc : `~mocpy.stmoc.STMoc`
        """
        store_index = mocpy.project_on_first_dim(smoc.store_index, stmoc.store_index)
        return cls(store_index)

    def _process_degradation(self, another_moc, order_op):
        """
        Degrade (down-sampling) self and ``another_moc`` to ``order_op`` order.

        Parameters
        ----------
        another_moc : `~mocpy.tmoc.TimeMoc`
        order_op : int
            the order in which self and ``another_moc`` will be down-sampled to.

        Returns
        -------
        result : (`~mocpy.tmoc.TimeMoc`, `~mocpy.tmoc.TimeMoc`)
            self and ``another_moc`` degraded TimeMocs

        """
        max_order = max(self.max_order, another_moc.max_order)
        if order_op > max_order:
            message = (
                "Requested time resolution for the operation cannot be applied.\n"
                f"The TimeMoc object resulting from the operation is of time resolution {TimeMOC.order_to_time_resolution(max_order).sec} sec."
            )
            warnings.warn(message, UserWarning, stacklevel=2)

        self_degradation = self.degrade_to_order(order_op)
        another_moc_degradation = another_moc.degrade_to_order(order_op)
        return self_degradation, another_moc_degradation

    def intersection_with_timeresolution(
        self,
        another_moc,
        delta_t=DEFAULT_OBSERVATION_TIME,
    ):
        """
        Intersection between self and moc.

        ``delta_t`` gives the possibility to the user
        to set a time resolution for performing the tmoc intersection

        Parameters
        ----------
        another_moc : `~mocpy.abstract_moc.AbstractMOC`
            the MOC/TimeMOC used for performing the intersection with self
        delta_t : `~astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMoc order to represent the observations. (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``)

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMOC`
            MOC object whose interval set corresponds to : self & ``moc``

        """
        order_op = TimeMOC.time_resolution_to_order(delta_t)

        self_degraded, moc_degraded = self._process_degradation(another_moc, order_op)
        return super(TimeMOC, self_degraded).intersection(moc_degraded)

    def union_with_timeresolution(self, another_moc, delta_t=DEFAULT_OBSERVATION_TIME):
        """
        Union between self and moc.

        ``delta_t`` gives the possibility to the user
        to set a time resolution for performing the tmoc union

        Parameters
        ----------
        another_moc : `~mocpy.abstract_moc.AbstractMOC`
            the MOC/TimeMoc to bind to self
        delta_t : `~astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMoc order to represent the observations. (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``)

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            MOC object whose interval set corresponds to : self | ``moc``

        """
        order_op = TimeMOC.time_resolution_to_order(delta_t)

        self_degraded, moc_degraded = self._process_degradation(another_moc, order_op)
        return super(TimeMOC, self_degraded).union(moc_degraded)

    def difference_with_timeresolution(
        self,
        another_moc,
        delta_t=DEFAULT_OBSERVATION_TIME,
    ):
        """
        Difference between self and moc.

        ``delta_t`` gives the possibility to the user to set a time resolution for
        performing the tmoc diff.

        Parameters
        ----------
        another_moc : `~mocpy.abstract_moc.AbstractMOC`
            the MOC/TimeMoc to substract from self
        delta_t : `~astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMoc order to represent the observations. (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``)

        Returns
        -------
        result : `~mocpy.moc.MOC` or `~mocpy.tmoc.TimeMoc`
            MOC object whose interval set corresponds to : self - ``moc``

        """
        order_op = TimeMOC.time_resolution_to_order(delta_t)

        self_degraded, moc_degraded = self._process_degradation(another_moc, order_op)
        return super(TimeMOC, self_degraded).difference(moc_degraded)

    @property
    def total_duration(self):
        """
        Get the total duration covered by the temporal moc.

        Returns
        -------
        duration : `~astropy.time.TimeDelta`
            total duration of all the observation times of the tmoc
            total duration of all the observation times of the tmoc

        """
        return TimeDelta(
            mocpy.ranges_sum(self.store_index) / 1e6,
            format="sec",
            scale="tdb",
        )

    @property
    def consistency(self):
        """
        Get a percentage of fill between the min and max time the moc is defined.

        A value near 0 shows a sparse temporal moc (i.e. the moc does not cover a lot
        of time and covers very distant times. A value near 1 means that the moc covers
        a lot of time without big pauses.

        Returns
        -------
        result : float
            fill percentage (between 0 and 1.)

        """
        return self.total_duration.jd / (self.max_time - self.min_time).jd

    @property
    def min_time(self):
        """
        Get the `~astropy.time.Time` time of the tmoc first observation.

        Returns
        -------
        min_time : `astropy.time.Time`
            time of the first observation

        """
        return microseconds_to_times(np.atleast_1d(self.min_index))

    @property
    def max_time(self):
        """
        Get the `~astropy.time.Time` time of the tmoc last observation.

        Returns
        -------
        max_time : `~astropy.time.Time`
            time of the last observation

        """
        return microseconds_to_times(np.atleast_1d(self.max_index))

    def contains(self, times, keep_inside=True):
        """
        Get a mask array (e.g. a numpy boolean array) of times being inside (or outside) the TMOC instance.

        Parameters
        ----------
        times : `astropy.time.Time`
            astropy times to check whether they are contained in the TMOC or not.
        keep_inside : bool, optional
            True by default. If so the filtered table contains only observations that are located the MOC.
            If ``keep_inside`` is False, the filtered table contains all observations lying outside the MOC.

        Returns
        -------
        array : `~numpy.darray`
            A mask boolean array
        """
        # the requested order for filtering the astropy observations table is more precise than the order
        # of the TimeMoc object
        pix_arr = times_to_microseconds(times)

        mask = mocpy.filter_time(self.store_index, pix_arr)

        if keep_inside:
            return mask
        return ~mask

    def contains_with_timeresolution(
        self,
        times,
        keep_inside=True,
        delta_t=DEFAULT_OBSERVATION_TIME,
    ):
        """
        Get a mask array (e.g. a numpy boolean array) of times being inside (or outside) the TMOC instance.

        Parameters
        ----------
        times : `astropy.time.Time`
            astropy times to check whether they are contained in the TMOC or not.
        keep_inside : bool, optional
            True by default. If so the filtered table contains only observations that are located the MOC.
            If ``keep_inside`` is False, the filtered table contains all observations lying outside the MOC.
        delta_t : `astropy.time.TimeDelta`, optional
            the duration of one observation. It is set to 30 min by default. This data is used to compute the
            more efficient TimeMOC order to represent the observations (Best order = the less precise order which
            is able to discriminate two observations separated by ``delta_t``).

        Returns
        -------
        array : `~numpy.darray`
            A mask boolean array
        """
        # the requested order for filtering the astropy observations table is more precise than the order
        # of the TimeMoc object
        current_max_order = self.max_order
        new_max_order = TimeMOC.time_resolution_to_order(delta_t)
        if new_max_order > current_max_order:
            message = (
                "Requested time resolution filtering cannot be applied.\n"
                f"Filtering is applied with a time resolution of {TimeMOC.order_to_time_resolution(current_max_order).sec} sec."
            )
            warnings.warn(message, UserWarning, stacklevel=2)

        rough_tmoc = self.degrade_to_order(new_max_order)
        return rough_tmoc.contains(times, keep_inside)

    @staticmethod
    def order_to_time_resolution(order):
        """
        Convert an TimeMoc order to its equivalent time.

        Parameters
        ----------
        order : int
            order to convert

        Returns
        -------
        delta_t : `~astropy.time.TimeDelta`
            time equivalent to ``order``

        """
        return TimeDelta(2 ** (61 - order) / 1e6, format="sec", scale="tdb")

    @staticmethod
    def time_resolution_to_order(delta_time):
        """
        Convert a time resolution to a TimeMoc order.

        Parameters
        ----------
        delta_time : `~astropy.time.TimeDelta`
            time to convert

        Returns
        -------
        order : int
            The less precise order which is able to discriminate two observations separated by ``delta_time``.

        """
        order = 61 - int(np.log2(delta_time.sec * 1e6))
        return np.uint8(order)

    def plot(self, title="TimeMoc", view=(None, None), figsize=(9.5, 5), **kwargs):
        """
        Plot the TimeMoc in a time window.

        This method uses interactive matplotlib. The user can move its mouse through the plot to see the
        time (at the mouse position).

        Parameters
        ----------
        title : str, optional
            The title of the plot. Set to 'TimeMoc' by default.
        view : (`~astropy.time.Time`, `~astropy.time.Time`), optional
            Define the view window in which the observations are plotted. Set to (None, None) by default (i.e.
            all the observation time window is rendered).

        """
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap

        if self.empty():
            import warnings

            warnings.warn("This time moc is empty", UserWarning, stacklevel=2)
            return

        plot_order = 30
        plotted_moc = (
            self.degrade_to_order(plot_order) if self.max_order > plot_order else self
        )

        min_jd = plotted_moc.min_time.jd if not view[0] else view[0].jd
        max_jd = plotted_moc.max_time.jd if not view[1] else view[1].jd

        if max_jd < min_jd:
            raise ValueError(
                f"Invalid selection: max_jd = {max_jd} must be > to min_jd = {min_jd}",
            )

        fig1 = plt.figure(figsize=figsize)
        ax = fig1.add_subplot(111)

        ax.set_xlabel("iso")
        ax.get_yaxis().set_visible(b=False)

        size = 2000
        delta = (max_jd - min_jd) / size
        min_jd_time = min_jd

        ax.set_xticks([0, size])
        ax.set_xticklabels(
            Time([min_jd_time, max_jd], format="jd", scale="tdb").iso,
            rotation=70,
        )

        y = np.zeros(size)
        for s_time_us, e_time_us in plotted_moc.to_time_ranges():
            s_index = int((s_time_us.jd - min_jd_time) / delta)
            e_index = int((e_time_us.jd - min_jd_time) / delta)
            y[s_index : (e_index + 1)] = 1.0

        # hack in case of full time mocs.
        if np.all(y):
            y[0] = 0

        z = np.tile(y, (int(size // 10), 1))

        plt.title(title)

        color_map = LinearSegmentedColormap.from_list("w2r", ["#fffff0", "#aa0000"])
        color_map.set_under("w")
        color_map.set_bad("gray")

        plt.imshow(z, interpolation="bilinear", **kwargs)

        def on_mouse_motion(event):
            for txt in ax.texts:
                txt.set_visible(b=False)

            text = ax.text(0, 0, "", va="bottom", ha="left")

            time = Time(event.xdata * delta + min_jd_time, format="jd", scale="tdb")

            tx = f"{time.iso}"
            text.set_position((event.xdata - 50, 700))
            text.set_rotation(70)
            text.set_text(tx)

        fig1.canvas.mpl_connect("motion_notify_event", on_mouse_motion)

        plt.show()

    @classmethod
    def load(cls, path, format="fits"):  # noqa: A002
        """
        Load the Time MOC from a file.

        Format can be 'fits', 'ascii', or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        path : str or pathlib.Path
            The path to the file to load the MOC from.
        format : str, optional
            The format from which the MOC is loaded.
            Possible formats are "fits", "ascii" or "json".
            By default, ``format`` is set to "fits".
        """
        path = str(path)
        if format == "fits":
            index = mocpy.time_moc_from_fits_file(path)
            return cls(index)
        if format == "ascii":
            index = mocpy.time_moc_from_ascii_file(path)
            return cls(index)
        if format == "json":
            index = mocpy.time_moc_from_json_file(path)
            return cls(index)
        formats = ("fits", "ascii", "json")
        raise ValueError(f"format should be one of {formats}")

    @classmethod
    def _from_fits_raw_bytes(cls, raw_bytes):
        """Load MOC from raw bytes of a FITS file."""
        index = mocpy.time_moc_from_fits_raw_bytes(raw_bytes)
        return cls(index)

    @classmethod
    def from_string(cls, value, format="ascii"):  # noqa: A002
        """
        Deserialize the Time MOC from the given string.

        Format can be 'ascii' or 'json', though the json format is not officially supported by the IVOA.

        WARNING: the serialization must be strict, i.e. **must not** contain overlapping elements

        Parameters
        ----------
        format : str, optional
            The format in which the MOC will be serialized before being saved.
            Possible formats are "ascii" or "json".
            By default, ``format`` is set to "ascii".
        """
        if format == "ascii":
            index = mocpy.time_moc_from_ascii_str(value)
            return cls(index)
        if format == "json":
            index = mocpy.time_moc_from_json_str(value)
            return cls(index)
        formats = ("ascii", "json")
        raise ValueError(f"format should be one of {formats}")
