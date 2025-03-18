# DESI Example Notebooks

The Astro Data Lab includes a collection of notebooks showcasing the DESI data. They are distributed in a few locations as follows.

## This folder: `03_ScienceExamples/DESI/`

- `01_Intro_to_DESI_DR1.ipynb` shows how to access the redshift catalog from the Astro Data Lab database, how to separate objects based on the DESI targeting information, how to access all the available spectra for a given object using [SPARCL (SPectra Analysis and Retrievable Catalog Lab)](https://astrosparcl.datalab.noirlab.edu), and finally how to plot the "best" spectrum.

- `01_Intro_to_DESI_EDR.ipynb` is the previous Early Data Release (EDR) version with similar functionality: how to access the redshift catalog from the Astro Data Lab database, how to separate objects based on the DESI targeting information, how to access all the available spectra for a given object using SPARCL, and finally how to plot the "best" spectrum. We recommend using the latest (DR1) version instead as it supersedes the EDR data.

- `01a_Intro_to_DESI_DR1-Py3.ipynb` is adapted from the `01_Intro_to_DESI_DR1.ipynb` to work without any DESI software. It requires a Python-3 environment with the datalab and sparclclient installed (both can be pip-installed locally by users if not working in the Jupyter server).

- `02_DESI_SDSS_Comparison.ipynb` shows how to use SPARCL data discovery to find available SDSS DR16 and DESI DR1 spectra of sources in a specific region of the sky with redshift and spectype constraints, how to retrieve and compare spectra for the same galaxy observed with both SDSS and DESI.

## How-to folder: `04_HowTos/`

- `QueryClient/How_to_query_DESI_DR1_Data.ipynb` demonstrates a variety of queries to the Astro Data Lab `desi_dr1` database.

- `QueryClient/How_to_query_DESI_EDR_Data.ipynb` demonstrates a variety of queries to the Astro Data Lab `desi_edr` database (note: superseded by the DR1 version).

- `SPARCL/How_to_use_SPARCL.ipynb` provides a basic introduction to using the SPARCL client (or sparclclient) to find and retrieve spectroscopic data within a Python notebook context. [(SPARCL = SPectra Analysis and Retrievable Catalog Lab)](https://astrosparcl.datalab.noirlab.edu)

- `SPARCL/Plot_SPARCL_Spectra_with_Jdaviz.ipynb` shows how to retrieve spectra from SPARCL and display them using the [Jdaviz](https://jdaviz.readthedocs.io/en/latest/index.html) data analysis visualization tool.

- `SPARCL/Plot_SPARCL_Spectra_with_Prospect.ipynb` shows how to retrieve spectra from SPARCL and display them using the [prospect](https://desi-prospect.readthedocs.io/en/latest/) interactive spectral visualization tool.

## Other useful references

### DESI Data Documentation Websites

The [DESI Landing Page](https://datalab.noirlab.edu/desi/index.php) at Astro Data Lab includes a brief overview of DESI and the early data release. The [Data access page](https://datalab.noirlab.edu/desi/access.php) describes multiple ways of accessing data through the Astro Data Lab or SPARCL.

The [DESI Data Documentation website](https://data.desi.lbl.gov/doc/) describes data access, data format, data releases, and includes links to technical papers as well as information on data license and acknowledgments. It is the overriding reference for official DESI information.

### Astro Data Lab & DESI Help

If you have an Astro Data Lab question, please visite the [Data Lab User Forum](https://datalab.noirlab.edu/help/).

If you have a DESI question, please visit the [DESI User Forum](https://help.desi.lbl.gov). 