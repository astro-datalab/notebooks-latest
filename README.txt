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


                       web: https://datalab.noirlab.edu
                    github: https://github.com/astro-datalab


                        Version of this file: 20230728
                        

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

This file contains information on:

- The set of default Jupyter notebooks contained in user accounts
- How to check for updated versions of the notebooks or for new notebooks
- How to contribute interesting notebooks to Data Lab for public use 

You can follow the order below if you are just getting started.

All notebooks are developed for Python 3. Furthermore, an HTML version
of the notebooks is included in order to show them fully rendered.

DEFAULT DATALAB NOTEBOOKS
=========================

01- GETTING STARTED

The notebooks in "01_GettingStartedWithDatalab/" provide a 101 intro
to Python, Jupyter and SQL, and show, for Data Lab, some basic steps
such as loading modules, authenticating, making a list of available
datasets, an example query, and an example image cutout. It also shows
how to obtain the statistics of catalog tables in order to determine
approximate row and column counts.

02- DATA ACCESS OVERVIEW

The notebook in "02_DataAccessOverview/" provides users with examples
of typical functions and commands to explore and use some of the main
datasets hosted by Astro Data Lab. It is a reference for scientific
applications, though not as detailed as the specific science examples
given below (item 03).

03- SCIENCE EXAMPLES

The "03_ScienceExamples/" folder contains notebooks that showcase
scientific applications using the datasets hosted at Data Lab. Each
science application contains at least one notebook, and each
survey/dataset is featured in at least one notebook. In some
instances, the same science case is featured with two or more surveys.

- DESI: introduction to the DESI EDR dataset at Data Lab and a comparison
  between SDSS and DESI spectra.

- DwarfGalaxies: discover dwarf galaxies as stellar overdensities in
  the DELVE DR1 and DR2, DES DR1, NSC DR1 and DR2, and SMASH datasets.

- EmLineGalaxies: two notebooks highlight how to obtain and stack
  spectra using the Data Lab spectro service, and how to detect
  outliers in the BPT diagnostic diagram.

- ExploringM31: explore the M31 galaxy with the PHAT dataset.

- GNIRS_DQS_SpectralInventory: show how to access the Gemini Near Infrared
  Spectrograph - Distant Quasar Survey (GNIRS-DQS) at Data Lab and
  example spectra plots.

- GOGREEN_GalaxiesInRichEnvironments: two notebooks showcase data
  access and image cutout services with the GOGREEN and GCLASS first
  data release, the first Gemini Large and Long program whose
  high-level science products are hosted at Data Lab.

- GalacticStructure: probe stellar populations in different parts of
  the Galactic Plane using the DECaPS dataset, and in the SMASH
  fields. Another notebook explores star clusters in Gaia,
  including animated visualizations.

- LargeScaleStructure: inspect large-scale structures using
  spectroscopic information from SDSS combined with photometric
  information from the DESI pre-imaging Legacy Surveys (LS).

- MagellanicClouds: examine the stellar substructures that surround
  the Magellanic Clouds using the VHS and Gaia datasets.
		     
- Pal5TidalTails: identify tidal tails of the globular cluster Palomar
  5 in the NSC catalog, as well as a jointly with Gaia to
  explore the proper motion of the cluster and its tails.

- SpectralEnergyDistributions: (1) use narrow-band filters to construct
  SEDs of objects from the S-PLUS dataset, and (2) compare the mid-
  infrared photometry from unWISE and AllWISE (3.4 & 4.6 micron).

- StarGalQSOSeparation: use photometric properties (colors,
  morphology/shape parameters, etc.) to distinguish between stars,
  galaxies, and QSOs in the DES and LS datasets.

- TimeSeriesAnalysisRrLyraeStar: analyze time-series to measure the
  period of RR Lyrae stars using photometry from SMASH.

- WhiteDwarfs: search for and analyze white dwarfs and other peculiar
  objects possibly ejected from the Galactic disk at very high velocities
  (> 400 km/s).

The ScienceExamples notebooks are located here:

   https://github.com/astro-datalab/notebooks-latest/tree/master/03_ScienceExamples/

04- HOW-TOS

The "04_HowTos/" folder contains sub-folders with notebooks that show
how to use Data Lab services with more detail than the brief examples
included in the GettingStarted and DataAccessOverview notebooks. The
functionality is shown for the full set of keywords and options for
the following:

- AuthClient: authenticating with the Data Lab
- CrossmatchTables: crossmatching a user-provided table against a table hosted by Data Lab, and using pre-crossmatched tables
- DataReduction: showing how to perform GMOS imaging data reduction and GMOS longslit spectroscopy data reduction using the Gemini DRAGONS and Gemini Pyraf packages
- FileService: using files rather than database tables, including spectra from SDSS/BOSS
- QueryClient: sending queries to the databases and retrieving results
- SPARCL: discovering, retrieving, analyzing, and plotting spectra (SDSS, BOSS, DESI) using the SPARCL service
- SiaService: obtaining cutouts using a Simple Image Access service
- StoreClient: storing data in virtual storage (VOSpace or MyDB)

The How-To notebooks are located here:

   https://github.com/astro-datalab/notebooks-latest/tree/master/04_HowTos/

05- CONTRIB

The "05_Contrib/" directory holds community-contributed notebooks to
Data Lab, including ANTARES example notebooks. Please see
./CONTRIBUTING file for detailed instructions.

06 - EPO

The "06_EPO/" directory provides "Education and Public Outreach"
notebooks, aimed at school students and teachers interested in
astronomical research and in teaching astronomy. The currently three sub-directories contain:

- TeenAstronomyCafe: notebooks originally developed for the
  "TeenAstronomy Café" activities organized jointly by the LSST and
  NOIRLab's outreach and engagement departments (middle/high school or
  undergraduate astronomy)

- e-TeenAstronomyCafe: same, but these notebooks are executable at
  Colab

- LaSerenaSchoolForDataScience: notebooks developed for the La Serena
  School for Data Science (undergraduate & graduate level), including
  machine learning, classification problems, and more.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

HOW TO UPDATE NOTEBOOKS
=======================

Data Lab never modifies the notebooks that were placed in your
notebooks/ directory during account creation. Over time, as the
default notebooks evolve, they will diverge from those in notebooks/.

To obtain a full copy of the newest default notebooks, click in the
top-left corner of the JupyterLab interface dashboard on the "+" icon,
which opens a new launcher page. Then, in the "Other" section, click
on the "Terminal" option, where you can use the `getlatest` function:

# without argument: copies the latest notebooks to a directory named with current date and time
username@datalab>getlatest
Copied /dlusers/username/notebooks-latest/ to notebooks_20211118_212650/

# with target directory as argument
username@datalab>getlatest mydir
Copied /dlusers/username/notebooks-latest/ to mydir/

All notebooks have a __version__ variable defined in the very first
cell. Simply running 'grep version foofile.ipynb' will show the
version of the given file.

Finally, copies of this README.txt file as well as the latest
notebooks are kept at the Data Lab GitHub account:
https://github.com/astro-datalab/notebooks-latest/ from where you can
freely clone them.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

DOCUMENTATION & RESOURCES
=========================

The User Manual includes a tutorial on using Jupyter Notebooks with the Data Lab:
https://datalab.noirlab.edu/docs/manual/UsingAstroDataLab/JupyterNotebooks/JupyterNotebooks.html

The User Manual also includes additional information on the Science Examples 
featured in the notebooks:
https://datalab.noirlab.edu/docs/manual/UsingAstroDataLab/ScienceExamples/index.html

Helpful advice on using SQL and writing queries can be found here: 
https://datalab.noirlab.edu/docs/manual/UsingAstroDataLab/SQLGotchas/SQLGotchas/SQLGotchas.html

Lastly, please visit the Helpdesk to see the FAQs or ask your questions: 
https://datalab.noirlab.edu/help/
