::::::::::::::::::'##::: ##::'#######:::::'###:::::'#######:::::::::::::::::::
:::::::::::::::::: ###:: ##:'##.... ##:::'## ##:::'##.... ##::::::::::::::::::
:::::::::::::::::: ####: ##: ##:::: ##::'##:. ##:: ##:::: ##::::::::::::::::::
:::::::::::::::::: ## ## ##: ##:::: ##:'##:::. ##: ##:::: ##::::::::::::::::::
:::::::::::::::::: ##. ####: ##:::: ##: #########: ##:::: ##::::::::::::::::::
:::::::::::::::::: ##:. ###: ##:::: ##: ##.... ##: ##:::: ##::::::::::::::::::
:::::::::::::::::: ##::. ##:. #######:: ##:::: ##:. #######:: ::::::::::::::::
..::::..:::.......:::..:::::..:::.......::::::::::::::::::::::::::::::::::::::
'########:::::'###::::'########::::'###:::::::'##::::::::::'###::::'########::
 ##.... ##:::'## ##:::... ##..::::'## ##:::::: ##:::::::::'## ##::: ##.... ##:
 ##:::: ##::'##:. ##::::: ##:::::'##:. ##::::: ##::::::::'##:. ##:: ##:::: ##:
 ##:::: ##:'##:::. ##:::: ##::::'##:::. ##:::: ##:::::::'##:::. ##: ########::
 ##:::: ##: #########:::: ##:::: #########:::: ##::::::: #########: ##.... ##:
 ##:::: ##: ##.... ##:::: ##:::: ##.... ##:::: ##::::::: ##.... ##: ##:::: ##:
 ########:: ##:::: ##:::: ##:::: ##:::: ##:::: ########: ##:::: ##: ########::
........:::..:::::..:::::..:::::..:::::..:::::........::..:::::..::........:::


           Welcome to the NOAO Data Lab Jupyter Notebook server!


                       web: https://datalab.noao.edu
                  github: https://github.com/noao-datalab


::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

This file contains information about the set of default Jupyter notebooks 
contained in user accounts, how to check for updated versions of the notebooks 
or for new notebooks. You can follow the order below if you are just getting 
started!

Generally, all notebooks should work with Python 3 (preferred) and Python 2. 
Furthermore, an html version of the notebooks is included in order to show 
them fully rendered.


1- GETTING STARTED

The notebook "GettingStartedWithDatalab" shows basic steps such as loading 
modules, authenticating, making a list of available datasets, an example query, 
and an example image cutout. It is meant as a reference for these general basic 
steps. It shows how to obtain the statistics of catalog tables in order to 
determine the row and column counts while avoiding to count all rows on large 
tables, which can be slow.


2- DATA ACCESS OVERVIEW

The notebook "DataAccessOverview" provides users with examples of typical 
functions and commands to explore and use some of the main datasets hosted by 
the Data Lab. It is a reference for scientific applications, though not as 
detailed as the specific science examples given below (item 4).


3- HOW-TOS

The "HowTos" folder contains sub-folders with notebooks that show how to use 
Data Lab services with more detail than the brief examples included in the 
Getting Started and Data Access Overview notebooks. The functionality is 
shown for the full set of keywords and options for the following:

- AuthClient: authenticating with the Data Lab
- QueryClient: sending queries to the databases and retrieving results
- RowVsCstore: using row-stored versus column-stored database tables
- SIA service: obtaining cutouts using a Simple Image Access service
- StoreClient: storing data in virtual storage (vospace or mydb)

The How-To notebooks are located here: 
   https://datalab.noao.edu/notebooks/web/HowTos/


4- SCIENCE EXAMPLES

The "ScienceExamples" folder contains notebooks that showcase scientific 
applications using the datasets hosted at Data Lab. Each science application 
contains at least one notebook, and each survey/dataset is featured in 
at least one notebook. In some instances, the same science case is featured 
with two or more surveys.

- DwarfGalaxies: discover dwarf galaxies as stellar overdensities in the 
                 DES DR1, NSC DR1 and SMASH datasets

- ExploringM31: explore the M31 galaxy with the PHAT dataset

- GalacticStructure: probe stellar populations in different parts of the 
                     Galactic Plane using the DECaPS dataset

- LargeScaleStructure: inspect large-scale structures using spectroscopic 
                       information from SDSS combined with photometric 
                       information from the DESI pre-imaging Legacy Survey (LS)

- StarGalQSOSeparation: use photometric properties (colors, morphology/shape 
                        parameters, etc.) to distinguish between stars, galaxies, 
                        and QSOs in the DES DR1 and LS DR3 datasets

- TimeSeriesAnalysisRrLyraeStar: analyze time-series to measure the period of 
                                 RR Lyrae stars using photometry from SMASH

- DECaPSBasicAccess: investigate Galactic structure by creating CMDs for 
		     stellar populations within the Galactic bulge, star 
		     clusters, and field.

The ScienceExamples notebooks are located here: 
   https://datalab.noao.edu/notebooks/web/ScienceExamples/


::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

HOW TO UPDATE NOTEBOOKS?

This README file as well as the notebooks listed herein are kept up-to-date 
on this Data Lab webpage: https://datalab.noao.edu/notebooks/web/

There will also be a copy on the Data Lab Github account (in progress): 
https://github.com/noao-datalab

The notebook naming convention typically includes the date in the filename, such 
as "ScienceCaseSurvey_YYYYMMDD.ipynb". Users can compare the content of their 
folders with the list of notebooks at either location above to determine whether 
there are new notebooks. If there are updated or new notebooks of interest, the 
current set-up requires that the users manually download them locally, and then 
upload them in their user accounts. We are developing an automated service that 
will allow users to check for updates and select notebooks, which will then be 
automatically transferred over. This service will be released at a later time.


::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

DOCUMENTATION & RESOURCES

The User Manual includes a tutorial on using Jupyter Notebooks with the Data Lab:
http://datalab.noao.edu/docs/manual/UsingTheNOAODataLab/JupyterNotebooks/JupyterNotebooks.html

The User Manual also includes additional information on the Science Examples 
featured in the notebooks:
http://datalab.noao.edu/docs/manual/UsingTheNOAODataLab/ScienceExamples/index.html

Helpful advice on using SQL and writing queries can be found here: 
http://datalab.noao.edu/docs/manual/UsingTheNOAODataLab/SQLGotchas/SQLGotchas.html

Lastly, please visit the Helpdesk to see the FAQs or ask your questions: 
http://datalab.noao.edu/help/

