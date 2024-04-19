
## DRAGONS data reduction (Jupyter Notebooks):

This repository has Jupyter Notebook examples of data reduction for the Gemini Observatory instruments. Usually, you need DRAGONS installed on your computer to run these notebooks, but the Astro Data Lab has a custom kernel called DRAGONS (Py3.7) that will allow you to run them. The notebooks are set to open the DRAGONS kernel by default, but if this doesn’t happen, you can click ‘kernel’ at the top of the notebook, then ‘change notebook,’ and select DRAGONS (Py3.7). These notebooks were written using the [DRAGONS' Application Program Interface (API)](https://dragons-recipe-system-users-manual.readthedocs.io/en/v2.1.1/appendices/full_api_example.html) for Python, based on the examples provided in the [DRAGONS Documentation](https://dragons.readthedocs.io/).

---
## Before you try the notebooks:

- [x] Install the Python-based DRAGONS platform by creating an anaconda environment. Detailed instructions can be found on the [Gemini Data Processing Software](https://www.gemini.edu//observing/phase-iii/understanding-and-processing-data/Data-Processing-Software) webpage.

- [x] At this point, you will have DRAGONS installed. To run the notebooks in the US NGO repository, you will have to allow the Jupyter notebooks to run on your Conda DRAGONS environment. You can do that by installing the packages listed below. From your bash terminal:

   ```
   conda install -n dragons nb_conda_kernels ipykernel
   ```
   
- [x] The Jupyter Notebooks are configured to download the raw (public) data directly from the [Gemini Observatory Archive](https://archive.gemini.edu/searchform). For that, you'll need the `wget` python package installed. From your bash terminal:

   ```
   conda install -n dragons wget
   ```

- [x] Download one of the files listed on the [Current notebooks available](#current-notebooks-available) section. Once the download is completed, start a bash terminal, go to the folder where you have downloaded the file, and type:

   ```
   conda activate dragons
   jupyter notebook [name].ipynb
   ```
   
*Note: at this stage, a browser window will open and the Jupyter notebook will be loaded. On the top right part of the window, you should see* `Python [conda env:dragons]` *or something similar, ensuring that the correct DRAGONS environment is loaded.*

---
## Current notebooks available:

### Flamingos2_Imaging_BrownDwarf.ipynb

Flamingos-2 imaging (Y-band) of the brown dwarf WISE J041358.14-475039.3. This is extracted from the Gemini/DRAGONS F2 tutorial, Section 3. Dataset includes Flats, Darks, and Science frames. Link to the [Jupyter notebook](notebooks-latest/04_HowTos/DataReduction/DRAGONS_reduction_examples
/Flamingos2_Imaging_BrownDwarf.ipynb) or [direct download](https://raw.githubusercontent.com/usngo/DRAGONS/main/Flamingos2_Imaging_BrownDwarf.ipynb) (right-click &#8594; Save Link As...).

### GMOS_Imaging_StarryField.ipynb

GMOS imaging (i-band) of a stellar field. This is extracted from the Gemini/DRAGONS GMOS tutorial, Section 3. Dataset includes Biases, Twilight Flats, and Science frames. Link to the [Jupyter notebook](GMOS_Imaging_StarryField.ipynb) or [direct download](https://raw.githubusercontent.com/usngo/DRAGONS/main/GMOS_Imaging_StarryField.ipynb) (right-click &#8594; Save Link As...).

### GMOS_Imaging_Galaxy.ipynb

GMOS imaging (g-band) of the elliptical galaxy NGC5018.  Dataset includes Biases, Twilight Flats, and Science frames. Link to the [Jupyter notebook](GMOS_Imaging_Galaxy.ipynb) or [direct download](https://raw.githubusercontent.com/usngo/DRAGONS/main/GMOS_Imaging_Galaxy.ipynb) (right-click &#8594; Save Link As...).

### GNIRS_Imaging_GammaRayBurst.ipynb

GNIRS imaging (J-band - point source through keyhole) of GRB120116A. This is extracted from the Gemini/DRAGONS GNIRS tutorial, Example 1-B. Dataset includes Flats, Darks, and Science frames. Link to the [Jupyter notebook](GNIRS_Imaging_GammaRayBurst.ipynb) or [direct download](https://raw.githubusercontent.com/usngo/DRAGONS/main/GNIRS_Imaging_GammaRayBurst.ipynb) (right-click &#8594; Save Link As...).

### GSAOI_Imaging_EllipticalGalaxy.ipynb

GSAOI imaging (K-short) of a field around NGC5128. This is extracted from the Gemini/DRAGONS GSAOI tutorial, Section 3. Dataset includes Flats, Standard Star, and Science frames. Link to the [Jupyter notebook](GSAOI_Imaging_EllipticalGalaxy.ipynb) or [direct download](https://raw.githubusercontent.com/usngo/DRAGONS/main/GSAOI_Imaging_EllipticalGalaxy.ipynb) (right-click &#8594; Save Link As...).

### NIRI_Imaging_Supernova.ipynb

NIRI imaging (K-prime) of SN2014J. This is extracted from the Gemini/DRAGONS NIRI tutorial, Section 4. Dataset includes Flats, Standard Star, Darks, and Science frames. Link to the [Jupyter notebook](NIRI_Imaging_Supernova.ipynb) or [direct download](https://raw.githubusercontent.com/usngo/DRAGONS/main/NIRI_Imaging_Supernova.ipynb) (right-click &#8594; Save Link As...).

### GMOS_longslit_WhiteDwarf.ipynb
GMOS longslit data of a candidate DB white dwarf J2145+0031. This example is based on the Gemini DRAGONS GMOS longslit tutorial, Section 3. This dataset includes arcs, biases, and flats for the standard star and science target. This tutorial contains commented-out code that enables interactive mode. These cells are not required for this example, but feel free to uncomment them to test the interactive features. Link to the [Jupyter notebook](GMOS_longslit_WhiteDwarf.ipynb) or [direct download] (https://raw.githubusercontent.com/usngo/DRAGONS/main/GMOS_longslit_WhiteDwarf.ipynb) (right-click &#8594; Save Link As...).

---
## Additional resources

- **DRAGONS Documentation:** Data reduction with DRAGONS is currently available for all imaging data from current Gemini instruments. New modes and instruments are being added. Detailed instructions on the current stable version can be found on the [documentation page](https://dragons.readthedocs.io/en/stable/).

- **Gemini DR Forum:** The Gemini Observatory maintains a [Data Reduction User Forum](http://drforum.gemini.edu/) with materials on data reduction, software installation, and much more!

---
## Need help?

Problems, comments, suggestions, and/or need help with setting up and running the Jupyter notebooks? You can contact the US NGO members via our [Portal](http://ast.noao.edu/csdc/usngo), [Twitter](https://twitter.com/usngo), or submit a *New issue* through github.

For assistance with DRAGONS installation and procedures, please submit a ticket to the [Gemini Helpdesk](https://www.gemini.edu/observing/helpdesk/submit-general-helpdesk-request) (Partner Country: US; Topic: DRAGONS).
---
