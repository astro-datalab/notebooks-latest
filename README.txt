::::::::::::::::'###:::::'######::'########:'########:::'#######::::::::::::::::
:::::::::::::::'## ##:::'##... ##:... ##..:: ##.... ##:'##.... ##:::::::::::::::
::::::::::::::'##:. ##:: ##:::..::::: ##:::: ##:::: ##: ##:::: ##:::::::::::::::
:::::::::::::'##:::. ##:. ######::::: ##:::: ########:: ##:::: ##:::::::::::::::
::::::::::::: #########::..... ##:::: ##:::: ##.. ##::: ##:::: ##:::::::::::::::
::::::::::::: ##.... ##:'##::: ##:::: ##:::: ##::. ##:: ##:::: ##:::::::::::::::
::::::::::::: ##:::: ##:. ######::::: ##:::: ##:::. ##:. #######::::::::::::::::
:::::::::::::..:::::..:::......::::::..:::::..:::::..:::.......:::::::::::::::::
:'########:::::'###::::'########::::'###:::::::'##::::::::::'###::::'########:::
: ##.... ##:::'## ##:::... ##..::::'## ##:::::: ##:::::::::'## ##::: ##.... ##::
: ##:::: ##::'##:. ##::::: ##:::::'##:. ##::::: ##::::::::'##:. ##:: ##:::: ##::
: ##:::: ##:'##:::. ##:::: ##::::'##:::. ##:::: ##:::::::'##:::. ##: ########:::
: ##:::: ##: #########:::: ##:::: #########:::: ##::::::: #########: ##.... ##::
: ##:::: ##: ##.... ##:::: ##:::: ##.... ##:::: ##::::::: ##.... ##: ##:::: ##::
: ########:: ##:::: ##:::: ##:::: ##:::: ##:::: ########: ##:::: ##: ########:::
:........:::..:::::..:::::..:::::..:::::..:::::........::..:::::..::........::::


           Welcome to the Astro Data Lab Jupyter Notebook repository


                       web: https://datalab.noao.edu
                    github: https://github.com/noaodatalab


                        Version of this file: 20191204
                        

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

This file contains information on:

- The set of default Jupyter notebooks contained in user accounts
- How to check for updated versions of the notebooks or for new notebooks
- How to contribute interesting notebooks to Data Lab for public use 

You can follow the order below if you are just getting started.

Generally, all notebooks should work with Python 3 and some may still 
work with Python 2 (we no longer develop for Python 2). Furthermore, an 
html version of the notebooks is included in order to show them fully
rendered.

DEFAULT DATALAB NOTEBOOKS
=========================

01- GETTING STARTED

The notebooks in "01_GettingStartedWithDatalab/" provide a 101 intro to
Python, Jupyter and SQL, and show for Data Lab basic steps such as
loading modules, authenticating, making a list of available datasets,
an example query, and an example image cutout. It is meant as a
reference for these general basic steps. It shows how to obtain the
statistics of catalog tables in order to determine the row and column
counts while avoiding to count all rows on large tables, which can be
slow.

02- DATA ACCESS OVERVIEW

The notebook in "02_DataAccessOverview/" provides users with examples
of typical functions and commands to explore and use some of the main
datasets hosted by the Data Lab. It is a reference for scientific
applications, though not as detailed as the specific science examples
given below (item 4).

03- SCIENCE EXAMPLES

The "03_ScienceExamples" folder contains notebooks that showcase
scientific applications using the datasets hosted at Data Lab. Each
science application contains at least one notebook, and each
survey/dataset is featured in at least one notebook. In some
instances, the same science case is featured with two or more surveys.

- DwarfGalaxies: discover dwarf galaxies as stellar overdensities in
  the DES DR1, NSC DR1 and SMASH datasets

- ExploringM31: explore the M31 galaxy with the PHAT dataset

- GalacticStructure: probe stellar populations in different parts of
  the Galactic Plane using the DECaPS dataset

- LargeScaleStructure: inspect large-scale structures using
  spectroscopic information from SDSS combined with photometric
  information from the DESI pre-imaging Legacy Survey (LS)
		     
- Pal5TidalTails: identify tidal tails of the globular cluster Palomar
  5 in the NSC DR1 catalog, as well as a jointly with Gaia DR2 to
  explore the proper motion of the cluster and its tails

- SpectralEnergyDistributions: (1) use narrow-band filters to construct
  SEDs of objects from the S-PLUS DR1 dataset, and (2) compare the mid-
  infrared photometry from unWISE DR1 and AllWISE (3.4 & 4.6 micron)

- StarGalQSOSeparation: use photometric properties (colors,
  morphology/shape parameters, etc.) to distinguish between stars,
  galaxies, and QSOs in the DES DR1 and LS DR7 datasets

- TimeSeriesAnalysisRrLyraeStar: analyze time-series to measure the
  period of RR Lyrae stars using photometry from SMASH

The ScienceExamples notebooks are located here:

   https://datalab.noao.edu/notebooks/web/ScienceExamples/

04- HOW-TOS

The "04_HowTos/" folder contains sub-folders with notebooks that show
how to use Data Lab services with more detail than the brief examples
included in the GettingStarted and DataAccessOverview notebooks. The
functionality is shown for the full set of keywords and options for
the following:

- AuthClient: authenticating with the Data Lab
- CrossmatchTables: crossmatching a user-provided table against a table hosted by Data Lab
- FileService: using files rather than database tables, including spectra from SDSS/BOSS
- QueryClient: sending queries to the databases and retrieving results
- RowVsCstore: using row-stored versus column-stored database tables
- SIA service: obtaining cutouts using a Simple Image Access service [future]
- StoreClient: storing data in virtual storage (vospace or mydb)

The How-To notebooks are located here:

   https://datalab.noao.edu/notebooks/web/HowTos/

05- CONTRIB

The "05_Contrib/" directory holds user-contributed notebooks to Data
Lab. Please see ./CONTRIBUTING file for detailed instructions.

06 - EPO

The "06_EPO/" directory provides "Education and Public Outreach"
notebooks, aimed at school students and teachers interested in
astronomical research and in teaching astronomy. The notebooks were
originally developed for the "Teen Astronomy Café" activities
organized jointly by the LSST and NOAO EPO departments.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

HOW TO UPDATE NOTEBOOKS
=======================

Data Lab never modifies the notebooks that were placed in your
notebooks/ directory during account creation. Over time, as the
default notebooks evolve, they will diverge from those in notebooks/.

To obtain a full copy of the newest default notebooks, click in the
top-right corner of the Jupyter dashboard on "New", then on "Terminal",
and use the `getlatest` function:

# without argument: copies to a directory named with current date and time
username@datalab>getlatest
Copied /dlusers/username/notebooks-latest/ to notebooks_20180709_212650/

# with target directory as argument
username@datalab>getlatest mydir
Copied /dlusers/username/notebooks-latest/ to mydir/

All notebooks have a __version__ variable defined in the very first
cell. Simply running 'grep version foofile.ipynb' will show the
version of the given file.

Finally, copies of this README.txt file as well as the latest
notebooks are kept on:

- the Data Lab website: https://datalab.noao.edu/notebooks/web/
- the Data Lab Github account: https://github.com/noaodatalab/notebooks-latest

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

DOCUMENTATION & RESOURCES
=========================

The User Manual includes a tutorial on using Jupyter Notebooks with the Data Lab:
https://datalab.noao.edu/docs/manual/UsingTheNOAODataLab/JupyterNotebooks/JupyterNotebooks.html

The User Manual also includes additional information on the Science Examples 
featured in the notebooks:
https://datalab.noao.edu/docs/manual/UsingTheNOAODataLab/ScienceExamples/

Helpful advice on using SQL and writing queries can be found here: 
https://datalab.noao.edu/docs/manual/UsingTheNOAODataLab/SQLGotchas/SQLGotchas.html

Lastly, please visit the Helpdesk to see the FAQs or ask your questions: 
https://datalab.noao.edu/help/
