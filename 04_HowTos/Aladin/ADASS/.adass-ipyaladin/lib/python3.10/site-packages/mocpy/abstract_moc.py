import abc
from collections.abc import Iterable
from functools import reduce
from io import BytesIO
from pathlib import Path

import numpy as np

from . import mocpy, serializer

__author__ = "Matthieu Baumann, Thomas Boch, Manon Marchand, François-Xavier Pineau"
__copyright__ = "CDS, Centre de Données astronomiques de Strasbourg"

__license__ = "BSD 3-Clause License"
__email__ = "matthieu.baumann@astro.unistra.fr, thomas.boch@astro.unistra.fr, manon.marchand@astro.unistra.fr, francois-xavier.pineau@astro.unistra.fr"


class AbstractMOC(serializer.IO, metaclass=abc.ABCMeta):
    """Is an abstract coverage MOC."""

    def __repr__(self):
        return self.to_string(format="ascii", fold=80)

    def __del__(self):
        """Erase MOC."""
        mocpy.drop(self.store_index)

    def __eq__(self, other):
        """
        Test equality between this a MOC instance and ``another_moc`` of same type.

        Parameters
        ----------
        other : `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The moc object to test the equality with

        Returns
        -------
        bool
            True if the self and ``other`` are equal.
        """
        if not isinstance(other, AbstractMOC):
            raise TypeError(
                f"Cannot compare a {type(self)} with a {type(other)}",
            )

        return mocpy.check_eq(self.store_index, other.store_index)

    def __copy__(self):
        mocpy.copy(self.store_index)
        return self.__class__(self.store_index)

    def __deepcopy__(self, memo):
        return self.__copy__()

    def __getstate__(self):
        return self.serialize(format="json")

    def __setstate__(self, state):
        # this is called when a MOC is unpickled
        # we create a new ref count with copy
        moc = self.from_json(state)
        mocpy.copy(moc.store_index)
        self.__dict__ = moc.__dict__

    @staticmethod
    def _store_index_dtype():
        """Store the datatype of the index.

        This function is used for portability between 32 and 64 bits systems.

        Returns
        -------
        _type_
            Type of the index, can only take the values np.uint64 or np.uint32
        """
        usize_n_bits = mocpy.usize_n_bits()
        if usize_n_bits == 64:
            return np.uint64
        if usize_n_bits == 32:
            return np.uint32
        raise ValueError("Unsupported store index usize type!")

    def __add__(self, moc):
        """Compute the union of self with another MOC.

        Operator +

        Parameters
        ----------
        moc :  `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            Another MOC to compute the union with.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The union of self and moc.
        """
        return self.union(moc)

    def __radd__(self, moc):
        """Compute the union of self with another MOC.

        Operator + (right side, here commutative operation)

        Parameters
        ----------
        moc :  `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            Another MOC to compute the union with.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The union of self and moc.
        """
        if moc == 0:
            # allows to use the sum() method on a list of MOCs
            return self
        return self.union(moc)

    def __or__(self, moc):
        """Compute the union of self with another MOC.

        Operator | (bitwise or)

        Parameters
        ----------
        moc :  `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            Another MOC to compute the union with.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The union of self and moc.
        """
        return self.union(moc)

    def __sub__(self, moc):
        """
        Operator - definition.

        Computes the difference of self with another MOC

        Parameters
        ----------
        moc :  `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            Another MOC to compute the difference with.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The difference of self with moc.
        """
        return self.difference(moc)

    def __and__(self, moc):
        """
        Operator & definition.

        Computes the intersection of self with another MOC

        Parameters
        ----------
        moc :  `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            Another MOC to compute the intersection with.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The intersection of self and moc.
        """
        return self.intersection(moc)

    def __invert__(self):
        """
        Unary operator ~ definition.

        Computes the complement of self

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The complement MOC of self.
        """
        return self.complement()

    def empty(self):
        """(e.g. a numpy boolean array).

        Returns
        -------
        bool
            True if the MOC instance is empty.
        """
        return mocpy.is_empty(self.store_index)

    @abc.abstractproperty
    def max_order(self):
        """Depth of the MOC instance."""

    @property
    def min_index(self):
        """Return the smallest index (at the deepest absolute resolution) the MOC contains."""
        return mocpy.first_index(self.store_index)

    @property
    def max_index(self):
        """Return the largest index (at the deepest absolute resolution) the MOC contains."""
        return mocpy.last_index(self.store_index)

    @property
    def uniq_gen(self):
        """
        Return a `np.array` of the generic uniq indices of the cell in the MOC.

        Warnings
        --------
        This is not defined in the MOC standard and is not HEALPix scpecific.

        Notes
        -----
        * It consists on the regular index with a sentinel bit placed at the immediate left
          of the index's MSB. At a given depth, the sentinel bit is always put o the same bit.
        * Because the uniq HEALPix encoding is not adapted for non-HEALPIx indices.
        """
        return mocpy.to_uniq_gen(self.store_index)

    @abc.abstractclassmethod
    def n_cells(cls, depth, dimension=None):
        """Give the number of cells for the given depth. This is defined in children classes."""

    @property
    def uniq_zorder(self):
        """
        Return a `np.array` of the zorder uniq indices of the cell in the MOC.

        Warnings
        --------
        This is not defined in the MOC standard and is not HEALPix specific.

        Notes
        -----
        * It consists on a regular index shifted on the left so that indices at all level have the same MSB.
          Plus a sentinel bit placed at the immediate right of the LSB.
        * Because the uniq HEALPix encoding is not adapted for non-HEALPIx indices
          AND because the natural ordering of such indices follow the same order as the zorder indices
          (which is very useful for streaming processing, e.g. when dealing with multi-order maps)
        """
        return mocpy.to_uniq_zorder(self.store_index)

    def flatten(self):
        """Return the list of indices of all cells in the MOC at the MOC depth."""
        return mocpy.flatten_to_moc_depth(self.store_index)

    def complement(self):
        """Return the complement of the MOC instance.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The resulting MOC.
        """
        index = mocpy.complement(self.store_index)

        return self.__class__(index)

    def intersection(self, another_moc, *mocs):
        """Intersection between the MOC instance and other MOCs.

        Parameters
        ----------
        another_moc : `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The MOC to do the intersection with.
        mocs : `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            Other additional MOCs to perform the intersection with.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The resulting MOC.

        Examples
        --------
        >>> from mocpy import FrequencyMOC
        >>> import astropy.units as u
        >>> fmoc_large_band = FrequencyMOC.from_frequency_ranges(order=42,
        ...                                                      min_freq=0.1*u.Hz,
        ...                                                      max_freq=100*u.Hz)
        >>> fmoc_sharp_band = FrequencyMOC.from_frequency_ranges(order=42,
        ...                                                      min_freq=10*u.Hz,
        ...                                                      max_freq=20*u.Hz)
        >>> fmoc_sharp_band.intersection(fmoc_large_band) == fmoc_sharp_band
        True
        """
        if mocs:
            store_indices = np.array(
                [self.store_index, another_moc.store_index]
                + [moc.store_index for moc in mocs],
                dtype=AbstractMOC._store_index_dtype(),
            )
            if isinstance(self.max_order, Iterable):  # is a composite MOC, ex: STMOCs
                # https://github.com/cds-astro/cds-moc-rust/blob/main/src/storage/u64idx/opn.rs#L96
                index = reduce(mocpy.intersection, store_indices)
            else:
                index = mocpy.multi_intersection(store_indices)
        else:  # case with only two mocs
            index = mocpy.intersection(self.store_index, another_moc.store_index)

        return self.__class__(index)

    def union(self, another_moc, *mocs):
        """Union between the MOC instance and other MOCs.

        Parameters
        ----------
        another_moc : `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The MOC to do the union with.
        mocs : `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            Other additional MOCs to perform the union with.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The resulting MOC.

        Examples
        --------
        >>> from mocpy import TimeMOC
        >>> from astropy.time import Time, TimeDelta
        >>> older = TimeMOC.from_time_ranges(min_times=Time('1999-01-01T00:00:00.123456789'),
        ...                                  max_times=Time('2005-01-01T00:00:00'),
        ...                                  delta_t = TimeDelta(1, format='jd')
        ...                                 )
        >>> newer = TimeMOC.from_time_ranges(min_times=Time('2000-01-01T00:00:00'),
        ...                                  max_times=Time('2010-01-01T00:00:00'),
        ...                                  delta_t = TimeDelta(1, format='jd')
        ...                                 )
        >>> union = older.union(newer) # == older + newer
        >>> print(union.min_time.jyear, union.max_time.jyear)
        [1998.99847987] [2010.00183614]
        """
        if mocs:
            store_indices = np.array(
                [self.store_index, another_moc.store_index]
                + [moc.store_index for moc in mocs],
                dtype=AbstractMOC._store_index_dtype(),
            )
            if isinstance(self.max_order, Iterable):  # is a composite MOC, ex: STMOCs
                index = reduce(mocpy.union, store_indices)
            else:
                index = mocpy.multi_union(store_indices)
        else:  # case with only two mocs
            index = mocpy.union(self.store_index, another_moc.store_index)

        return self.__class__(index)

    def symmetric_difference(self, another_moc, *mocs):
        """Symmetric difference (XOR) between the MOC instance and other MOCs.

        a XOR b == (a and not b) or (not a and b)
        It is not implemented yet for STMOCs

        Parameters
        ----------
        another_moc : `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The MOC used that will be subtracted to self.
        mocs : `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            Other additional MOCs to perform the difference with.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The resulting MOC.

        Examples
        --------
        >>> from mocpy import MOC
        >>> moc1 = MOC.from_string("3/0-1 362-363")
        >>> moc2 = MOC.from_string("3/0 2 277 279")
        >>> moc1.symmetric_difference(moc2)
        3/1-2 277 279 362-363
        """
        if isinstance(self.max_order, Iterable):  # is a composite MOC, ex: STMOCs
            raise NotImplementedError(
                "Symmetric difference is not implemented yet for Space-Time MOCs",
            )
        if mocs:
            store_indices = np.array(
                [self.store_index, another_moc.store_index]
                + [moc.store_index for moc in mocs],
                dtype=AbstractMOC._store_index_dtype(),
            )
            index = mocpy.multi_symmetric_difference(store_indices)
        else:  # case with only two mocs
            index = mocpy.symmetric_difference(
                self.store_index,
                another_moc.store_index,
            )
        return self.__class__(index)

    def difference(self, another_moc, *mocs):
        """Difference between the MOC instance and other MOCs.

        Parameters
        ----------
        another_moc : `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The MOC used that will be subtracted to self.
        mocs : `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            Other additional MOCs to perform the difference with.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The resulting MOC.

        Examples
        --------
        >>> from mocpy import MOC
        >>> moc1 = MOC.from_string("3/0-7")
        >>> moc2 = MOC.from_string("3/0-3")
        >>> moc3 = MOC.from_string("3/4-7")
        >>> moc1.difference(moc2, moc3) # should the empty MOC of order 3 (3/)
        3/
        """
        if mocs:
            store_indices = np.array(
                [self.store_index, another_moc.store_index]
                + [moc.store_index for moc in mocs],
                dtype=AbstractMOC._store_index_dtype(),
            )
            if isinstance(self.max_order, Iterable):  # is a composite MOC, ex: STMOCs
                index_union = reduce(mocpy.union, store_indices)
            else:  # case with only two mocs
                index_union = mocpy.multi_union(store_indices)
            index = mocpy.difference(self.store_index, index_union)
            mocpy.drop(index_union)
        else:
            index = mocpy.difference(self.store_index, another_moc.store_index)

        return self.__class__(index)

    def extended(self):
        """Return the MOC extended by the external border made of cells at the MOC maximum depth.

        The only difference with respect to `add_neighbours` is that `extended` returns a new MOC
        instead of modifying the existing one.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The extended MOC
        """
        index = mocpy.extend(self.store_index)
        return self.__class__(index)

    def contracted(self):
        """Return the MOC contracted by removing the internal border made of cells at the MOC maximum depth.

        The only difference with respect to `remove_neighbours` is that `contracted` returns a new MOC
        instead of modifying the existing one.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The extended MOC
        """
        index = mocpy.contract(self.store_index)
        return self.__class__(index)

    def add_neighbours(self):
        """Extend the MOC instance so that it includes the HEALPix cells touching its border.

        The depth of the HEALPix cells added at the border is equal to the maximum depth of the MOC instance.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            self extended by one degree of neighbours.
        """
        prevstore_index = self.store_index
        self.store_index = mocpy.extend(prevstore_index)
        return self

    def remove_neighbours(self):
        """Remove from the MOC instance the HEALPix cells located at its border.

        The depth of the HEALPix cells removed is equal to the maximum depth of the MOC instance.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            self minus its HEALPix cells located at its border.
        """
        prevstore_index = self.store_index
        self.store_index = mocpy.contract(prevstore_index)
        return self

    @abc.abstractclassmethod
    def load(cls):
        """Load a MOC, has to be defined in children classes."""

    @classmethod
    def from_json(cls, json_moc):
        """Create a MOC from a dictionary of HEALPix cell arrays indexed by their depth.

        Parameters
        ----------
        json_moc : dict(str : [int]
            A dictionary of HEALPix cell arrays indexed by their depth.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            the MOC.
        """
        import json

        return cls.from_string(
            json.dumps(json_moc, sort_keys=True, indent=2),
            format="json",
        )

    @classmethod
    def from_fits(cls, path_or_url, timeout=1000):
        """Load a MOC from a FITS file.

        The specified FITS file must store the MOC
        (i.e. the list of HEALPix cells it contains)
        in a binary HDU table.

        Parameters
        ----------
        path : str
            The path to the FITS file.
        timeout : float
            Timeout for the query, defaults to 1000s

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The resulting MOC.

        """
        if isinstance(path_or_url, BytesIO):
            return cls._from_fits_raw_bytes(path_or_url.read())

        # this try except clause is there to support
        # Windows users with python 3.7 and should be dropped
        # when we remove support of python 3.7
        try:
            if Path(path_or_url).is_file():
                return cls.load(path_or_url, format="fits")
        except OSError:
            pass

        try:
            import requests

            response = requests.get(
                path_or_url,
                headers={"User-Agent": "MOCPy"},
                timeout=timeout,
            )
            if response:
                raw_bytes = BytesIO()
                raw_bytes.write(response.content)
                raw_bytes.seek(0)
                return cls.from_fits(raw_bytes)
            response.raise_for_status()
        except ModuleNotFoundError:
            raise UserWarning(
                "The `requests` module is required to fetch FITS files from an url",
            ) from ModuleNotFoundError

    @classmethod
    def from_str(cls, value):
        """Create a MOC from a string.

        This grammar is expressed is the `MOC IVOA <http://ivoa.net/documents/MOC/20190215/WD-MOC-1.1-20190215.pdf>`__
        specification at section 2.3.2.

        Parameters
        ----------
        value : str
            The MOC as a string following the grammar rules.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The resulting MOC

        Examples
        --------
        >>> from mocpy import MOC
        >>> moc = MOC.from_str("2/2-25 28 29 4/0 6/")

        See Also
        --------
        from_string: a duplicate of this method, with added ``fold`` option
        """
        # TODO : decide if we want to remove this duplicated method
        return cls.from_string(value, format="ascii")

    def degrade_to_order(self, new_order):
        """Degrade the MOC instance to a new, less precise, MOC.

        The maximum depth (i.e. the depth of the smallest cells that can be found in the MOC) of the
        degraded MOC is set to ``new_order``.

        Parameters
        ----------
        new_order : int
            Maximum depth of the output degraded MOC.

        Returns
        -------
        `mocpy.MOC`, `mocpy.TimeMOC`, `mocpy.FrequencyMOC`, `mocpy.STMOC`
            The degraded MOC.
        """
        raise NotImplementedError("Method degrade_to_order not implemented")

    def to_string(self, format="ascii", fold=0):  # noqa: A002
        """Write the MOC into a string.

        Format can be 'ascii' or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        format : str, optional
            The format in which the MOC will be serialized before being saved.
            Possible formats are "ascii" or "json".
            By default, ``format`` is set to "ascii".
        fold: int
            if >0, print ascii or json output with a maximum line width
        """
        if format == "ascii":
            if fold > 0:
                return mocpy.to_ascii_str_with_fold(self.store_index, fold)[:-1]
            return mocpy.to_ascii_str(self.store_index)[:-1]
        if format == "json":
            if fold > 0:
                return mocpy.to_json_str_with_fold(self.store_index, fold)
            return mocpy.to_json_str(self.store_index)
        formats = ("ascii", "json")
        raise ValueError(f"format should be one of {formats}")

    def save(
        self,
        path,
        format="fits",  # noqa: A002
        overwrite=False,
        pre_v2=False,
        fold=0,
        fits_keywords=None,
    ):
        """Write the MOC to a file.

        Format can be 'fits', 'ascii', or 'json', though the json format is not officially supported by the IVOA.

        Parameters
        ----------
        path : str or pathlib.Path
            The path to the file to save the MOC in.
        format : str, optional
            The format in which the MOC is saved.
            Possible formats are "fits", "ascii" or "json".
            By default, ``format`` is set to "fits".
        overwrite : bool, optional
            If the file already exists and you want to overwrite it, then set the  ``overwrite`` keyword.
            Default to False.
        fold: int
            if >0, print ascii or json output with a maximum line width
        fits_keywords: dict, optional
            Additional keywords to add to the FITS header.
        """
        path = Path(path)
        if path.is_file() and not overwrite:
            raise OSError(
                f"File '{path}' already exists! Set ``overwrite`` to "
                "True if you want to replace it.",
            )

        if format == "fits":
            if fold != 0:
                raise ValueError(
                    "The ``fold`` argument is only applied for the ``json`` and ``ascii`` output formats.",
                )
            from astropy.io import fits

            hdulist = fits.HDUList.fromstring(
                mocpy.to_fits_raw(self.store_index, pre_v2),
            )
            hdu = hdulist[1]
            if fits_keywords:
                for key in fits_keywords:
                    hdu.header[key] = fits_keywords[key]
            hdulist.writeto(path, overwrite=overwrite)
            return

        if fits_keywords:
            raise ValueError(
                "The ``fits_keyword`` argument is only valid when ``fits`` is the output format.",
            )

        if format == "ascii":
            if fold <= 0:
                mocpy.to_ascii_file(self.store_index, str(path))
            else:
                mocpy.to_ascii_file_with_fold(self.store_index, str(path), fold)
            return

        if format == "json":
            if fold <= 0:
                mocpy.to_json_file(self.store_index, str(path))
            else:
                mocpy.to_json_file_with_fold(self.store_index, str(path), fold)
            return

        formats = ("fits", "ascii", "json")
        raise ValueError(
            f"'format' should be one of {formats} but '{format}' was encountered.",
        )
