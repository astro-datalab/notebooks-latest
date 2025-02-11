from astropy import version as astropy_version


def check_astropy_version():
    """
    Check if the installed version of astropy is compatible with MIVOT.
    """
    if not astropy_version.version:
        return False
    if astropy_version.version < "6.0":
        print(f"Astropy version {astropy_version.version} is below "
              f"the required version 6.0 for the use of MIVOT.")
        return False
    return True
