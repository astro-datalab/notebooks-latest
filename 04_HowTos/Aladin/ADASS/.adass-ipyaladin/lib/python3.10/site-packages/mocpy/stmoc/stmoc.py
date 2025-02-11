import numpy as np

from .. import MOC, mocpy
from ..abstract_moc import AbstractMOC
from ..tmoc import TimeMOC, microseconds_to_times, times_to_microseconds

__author__ = "Matthieu Baumann, Thomas Boch, Manon Marchand, François-Xavier Pineau"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "matthieu.baumann@astro.unistra.fr, thomas.boch@astro.unistra.fr, manon.marchand@astro.unistra.fr, francois-xavier.pineau@astro.unistra.fr"


class STMOC(AbstractMOC):
    """Time-Spatial Coverage class."""

    def __init__(self, store_index):
        """Is a Time-Spatial Coverage (ST-MOC).

        Args:
            create_key: Object ensure __init__ is called by super-class/class-methods only
            store_index: index of the ST-MOC in the rust-side storage
        """
        self.store_index = store_index

    def __del__(self):
        """Erase STMOC."""
        super().__del__()

    def __eq__(self, other):
        """Assert equality between MOCs."""
        return super().__eq__(other)

    @property
    def max_depth(self):
        """Return max depth of MOC."""
        return mocpy.coverage_2d_depth(self.store_index)

    @property
    def max_order(self):
        """Is a clone of max_depth, to preserve the api between moc types."""
        return self.max_depth

    @property
    def max_time(self):
        """Return STMOC max time."""
        return microseconds_to_times(mocpy.coverage_2d_max_time(self.store_index))

    @property
    def min_time(self):
        """Return STMOC min time."""
        return microseconds_to_times(mocpy.coverage_2d_min_time(self.store_index))

    @classmethod
    def n_cells(cls, depth, dimension):
        """Get the number of cells for a given depth.

        Parameters
        ----------
        depth : int
            The depth. It is comprised between 0 and `~mocpy.moc.MOC.MAX_ORDER` if
            dimension='space' and between 0 and `~mocpy.tmoc.TimeMOC.MAX_ORDER` if
            dimension='time'.

        dimension : str
            Can be either 'time' or 'space'.

        Returns
        -------
        int
            The number of cells at the given order

        Examples
        --------
        >>> from mocpy import STMOC
        >>> STMOC.n_cells(0, dimension='space')
        12
        """
        if dimension == "space":
            return MOC.n_cells(depth)
        if dimension == "time":
            return TimeMOC.n_cells(depth)
        raise ValueError(
            f"Dimension should be either 'time' of 'space' but '{dimension}' was provided.",
        )

    @classmethod
    def new_empty(cls, max_depth_time, max_depth_space):
        """Create a new empty STMOC.

        Parameters
        ----------
        max_depth_time : int
            The time resolution of the STMOC. Should be comprised between 0 and 61.
        max_depth_space : int
            The space resolution of the STMOC. Should be comprised between 0 and 29.

        Returns
        -------
        `~mocpy.moc.STMOC`

        Examples
        --------
        >>> from mocpy import STMOC
        >>> STMOC.new_empty(42, 12)
        t42/ s12/

        """
        index = mocpy.new_empty_stmoc(
            np.uint8(max_depth_time),
            np.uint8(max_depth_space),
        )
        return cls(index)

    def is_empty(self):
        """Check whether the Space-Time coverage is empty."""
        return mocpy.is_empty(self.store_index)

    @classmethod
    def from_times_positions(cls, times, time_depth, lon, lat, spatial_depth):
        """
        Create a 2D Coverage from a set of times and positions associated to each time.

        - Its first dimension refers to `astropy.time.Time` times.
        - Its second dimension refers to lon, lat `astropy.units.Quantity` positions.

        Parameters
        ----------
        time : `astropy.time.Time`
            The times of each sky coordinates.
        time_depth : int
            Time depth.
        lon : `astropy.units.Quantity`
            The longitudes of the sky coordinates observed at a specific time.
        lat : `astropy.units.Quantity`
            The latitudes of the sky coordinates observed at a specific time.
        spatial_depth : int
            Spatial depth.

        Returns
        -------
        result : `~mocpy.stmoc.STMOC`
            The resulting Spatial-Time Coverage map.
        """
        times = times_to_microseconds(times)
        lon = lon.to_value("rad").astype(np.float64)
        lat = lat.to_value("rad").astype(np.float64)

        if times.shape != lon.shape or lon.shape != lat.shape:
            raise ValueError("Times and positions must have the same length.")

        if times.ndim != 1:
            raise ValueError("Times and positions must be 1D arrays.")

        index = mocpy.from_time_lonlat(times, time_depth, lon, lat, spatial_depth)

        return cls(index)

    @classmethod
    def from_time_ranges_positions(
        cls,
        times_start,
        times_end,
        lon,
        lat,
        time_depth=61,
        spatial_depth=29,
    ):
        """
        Create a 2D Coverage from a set of times and positions associated to each time.

        - Its first dimension refers to `astropy.time.Time` times.
        - Its second dimension refers to lon, lat `astropy.units.Quantity` positions.

        Parameters
        ----------
        times_start : `astropy.time.Time`
            The starting times of each observations.
        times_end : `astropy.time.Time`
            The ending times of each observations.
        lon : `astropy.units.Quantity`
            The longitudes of the sky coordinates observed at a specific time.
        lat : `astropy.units.Quantity`
            The latitudes of the sky coordinates observed at a specific time.
        time_depth : int, optional
            Time depth. By default, the time resolution chosen is 1µs.
        spatial_depth : int, optional
            Spatial depth. By default, the space resolution chosen is 393.2μas.

        Returns
        -------
        result : `~mocpy.stmoc.STMOC`
            The resulting Spatial-Time Coverage map.
        """
        # times_start = times_start.jd.astype(np.float64)
        # times_end = times_end.jd.astype(np.float64)

        times_start = times_to_microseconds(times_start)
        times_end = times_to_microseconds(times_end)

        lon = lon.to_value("rad").astype(np.float64)
        lat = lat.to_value("rad").astype(np.float64)

        if (
            times_start.shape != lon.shape
            or lon.shape != lat.shape
            or times_start.shape != times_end.shape
        ):
            raise ValueError("Times and positions must have the same length.")

        if times_start.ndim != 1:
            raise ValueError("Times and positions must be 1D arrays.")

        index = mocpy.from_time_ranges_lonlat(
            times_start,
            times_end,
            time_depth,
            lon,
            lat,
            spatial_depth,
        )

        return cls(index)

    @classmethod
    def from_spatial_coverages(
        cls,
        times_start,
        times_end,
        spatial_coverages,
        time_depth=61,
    ):
        """
        Create a 2D Coverage from a set of time ranges and spatial coverages associated to them.

        Parameters
        ----------
        times_start : `astropy.time.Time`
            The starting times of each observations.
        times_end : `astropy.time.Time`
            The ending times of each observations.
        spatial_coverages : list
            List of `mocpy.MOC` spatial coverage objects.
        time_depth : int, optional
            Time depth. By default, the time resolution chosen is 1µs.

        Returns
        -------
        result : `~mocpy.stmoc.STMOC`
            The resulting Spatial-Time Coverage map.
        """
        # accept also when there is a single spatial moc
        times_start = np.atleast_1d(times_start)
        times_end = np.atleast_1d(times_end)
        spatial_coverages = np.atleast_1d(spatial_coverages)

        times_start = times_to_microseconds(times_start)
        times_end = times_to_microseconds(times_end)

        if times_start.shape != times_end.shape or times_start.shape[0] != len(
            spatial_coverages,
        ):
            raise ValueError(
                "Time ranges and spatial coverages must have the same length",
            )

        if times_start.ndim != 1:
            raise ValueError("Times and spatial coverages must be 1D arrays")

        spatial_coverages_indices = np.fromiter(
            (arg.store_index for arg in spatial_coverages),
            dtype=AbstractMOC._store_index_dtype(),
        )
        index = mocpy.from_time_ranges_spatial_coverages(
            times_start,
            times_end,
            time_depth,
            spatial_coverages_indices,
        )

        return cls(index)

    def query_by_time(self, tmoc):
        """
        Query the ST-MOC by time T-MOC.

        This will perform the union of all the spatial coverages lying in a set of time ranges.

        Parameters
        ----------
        tmoc : ~mocpy.tmoc.TimeMOC``
            Time ranges. Must be a Nx2 shaped astropy time array.

        Returns
        -------
        result : `~mocpy.moc.MOC`
            The spatial coverage being observed within the input time ranges
        """
        # if times.ndim != 2 or times.shape[1] != 2:
        #    raise ValueError(
        #        "Times ranges must be provided. The shape of times must be (_, 2)"
        #    )

        # times = np.asarray(times.jd * 86400000000, dtype=np.uint64)
        # times = utils.times_to_microseconds(times)
        # ranges = mocpy.project_on_second_dim(times, self.store_index)
        # return MOC(IntervalSet(ranges, make_consistent=False))
        # store_index = from_stmoc_time_fold(cls, tmoc, stmoc):
        return MOC.from_stmoc_time_fold(tmoc, self)

    def query_by_space(self, smoc):
        """
        Query the ST-MOC by space coverage.

        This will perform the union of all the time ranges whose associated
        spatial coverages lie in ``moc``.

        Parameters
        ----------
        smoc : `~mocpy.moc.MOC`
            The spatial coverage.

        Returns
        -------
        result : `~mocpy.tmoc.TimeMOC`
            The time ranges observing in the ``spatial_coverage``
        """
        return TimeMOC.from_stmoc_space_fold(smoc, self)

    def contains(self, times, lon, lat, inside=True):
        """
        Return a boolean mask array of the (times, positions) lying inside (or outside) the Space-Time coverage.

        Parameters
        ----------
        times : `astropy.time.Time`
            The times of each sky coordinates.
        lon : `astropy.units.Quantity`
            The longitudes of the sky coordinates observed at a specific time.
        lat : `astropy.units.Quantity`
            The latitudes of the sky coordinates observed at a specific time.
        inside : bool, optional
            True by default. The returned mask array has true values for (time, position)
            lying inside the Space-Time coverage.

        Raises
        ------
        ValueError : If `times`, `lon` and `lat` do not have the same length.

        Returns
        -------
        array : `~np.ndarray`
            A mask boolean array
        """
        # times = times.jd.astype(np.float64)
        times = times_to_microseconds(times)
        lon = lon.to_value("rad").astype(np.float64)
        lat = lat.to_value("rad").astype(np.float64)

        if times.shape != lon.shape or lon.shape != lat.shape:
            raise ValueError("Times and positions must have the same length.")

        if times.ndim != 1:
            raise ValueError("Times and positions must be 1D arrays.")

        result = mocpy.coverage_2d_contains(self.store_index, times, lon, lat)

        if not inside:
            result = ~result

        return result

    @classmethod
    # A002: Argument `format` is shadowing a python function
    def load(cls, path, format="fits"):  # noqa: A002
        """
        Load the Spatial MOC from a file.

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
            index = mocpy.coverage_2d_from_fits_file(path)
            return cls(index)
        if format == "ascii":
            index = mocpy.coverage_2d_from_ascii_file(path)
            return cls(index)
        if format == "json":
            index = mocpy.coverage_2d_from_json_file(path)
            return cls(index)
        formats = ("fits", "ascii", "json")
        raise ValueError(f"format should be one of {formats}")

    @classmethod
    def _from_fits_raw_bytes(cls, raw_bytes):
        """Load MOC from raw bytes of a FITS file."""
        index = mocpy.stmoc_from_fits_raw_bytes(raw_bytes)
        return cls(index)

    @classmethod
    # A002: Argument `format` is shadowing a python function
    def from_string(cls, value, format="ascii"):  # noqa: A002
        """
        Deserialize the Spatial MOC from the given string.

        Format can be 'ascii' or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        format : str, optional
            The format in which the MOC will be serialized before being saved.
            Possible formats are "ascii" or "json".
            By default, ``format`` is set to "ascii".
        """
        if format == "ascii":
            index = mocpy.coverage_2d_from_ascii_str(value)
            return cls(index)
        if format == "json":
            index = mocpy.coverage_2d_from_json_str(value)
            return cls(index)
        formats = ("ascii", "json")
        raise ValueError(f"format should be one of {formats}")
