from astropy.io import fits

from . import mocpy


class IO:
    """Input and outputs for MOCs."""

    def serialize(self, format="fits", optional_kw_dict=None, pre_v2=False):
        """
        Serialize the MOC into a specific format.

        Possible formats are FITS, JSON and STRING

        Parameters
        ----------
        format : str
            'fits' by default. The other possible choice is 'json' or 'str'.
        optional_kw_dict : dict
            Optional keywords arguments added to the FITS header. Only used if ``format`` equals to 'fits'.

        Returns
        -------
        result : `astropy.io.fits.HDUList` or JSON dictionary
            The result of the serialization.
        """
        formats = ("fits", "json", "str")
        if format not in formats:
            raise ValueError(f"format should be one of {formats}")

        if format == "fits":
            hdulist = fits.HDUList.fromstring(
                mocpy.to_fits_raw(self.store_index, pre_v2),
            )
            hdu = hdulist[1]
            if optional_kw_dict:
                for key in optional_kw_dict:
                    hdu.header[key] = optional_kw_dict[key]
            return hdulist

        if format == "str":
            return self.to_string(format="ascii", fold=0)

        import json

        json_str = self.to_string(format="json")
        return json.loads(json_str)

    def write(
        self,
        path,
        format="fits",
        overwrite=False,
        optional_kw_dict=None,
        pre_v2=False,
    ):
        """
        Write the MOC to a file.

        Format can be 'fits' or 'json', though only the fits format is officially supported by the IVOA.

        Parameters
        ----------
        path : str
            The path to the file to save the MOC in.
        format : str, optional
            The format in which the MOC will be serialized before being saved. Possible formats are "fits" or "json".
            By default, ``format`` is set to "fits".
        overwrite : bool, optional
            If the file already exists and you want to overwrite it, then set the  ``overwrite`` keyword. Default to False.
        optional_kw_dict : optional
            Optional keywords arguments added to the FITS header. Only used if ``format`` equals to 'fits'.
        """
        import warnings

        warnings.warn(
            'This method is deprecated. Use MOC.save(path, "fits") instead!',
            DeprecationWarning,
            stacklevel=2,
        )
        self.save(path, format, overwrite, pre_v2, fits_keywords=optional_kw_dict)
