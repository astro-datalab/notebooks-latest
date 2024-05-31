
## DRAGONS data reduction (Jupyter Notebooks):

This repository has Jupyter Notebook examples of data reduction for the Gemini Observatory instruments. Usually, you need DRAGONS installed on your computer to run these notebooks, but the Astro Data Lab has a custom kernel called DRAGONS (Py3.7) that will allow you to run them. The notebooks are set to open the DRAGONS kernel by default, but if this doesn’t happen, click on the name of the current kernel in the top-right corner and select DRAGONS (Py3.7). These notebooks were written using the [DRAGONS' Application Program Interface (API)](https://dragons.readthedocs.io/projects/recipe-system-users-manual/en/release-3.2.x/appendices/full_api_example.html) for Python, based on the examples provided in the [DRAGONS Documentation](https://dragons.readthedocs.io/).

---
## Current notebooks available:

### Flamingos2_Imaging_BrownDwarf.ipynb

Flamingos-2 imaging (Y-band) of the brown dwarf WISE J041358.14-475039.3. This is extracted from the Gemini/DRAGONS F2 tutorial, Section 3. Dataset includes flats, darks, and science frames. Link to the [Jupyter notebook](https://github.com/astro-datalab/notebooks-latest/tree/master/04_HowTos/DataReduction/DRAGONS_reduction_examples/Flamingos2_Imaging_BrownDwarf/Flamingos2_Imaging_BrownDwarf.ipynb).

### GMOS_Imaging_StarryField.ipynb

GMOS imaging (i-band) of a stellar field. This is extracted from the Gemini/DRAGONS GMOS tutorial, Section 3. Dataset includes biases, twilight flats, and science frames. Link to the [Jupyter notebook](https://github.com/astro-datalab/notebooks-latest/tree/master/04_HowTos/DataReduction/DRAGONS_reduction_examples/GMOS_Imaging_StarryField/GMOS_Imaging_StarryField.ipynb).

### GMOS_Imaging_Galaxy.ipynb

GMOS imaging (g-band) of the elliptical galaxy NGC5018.  Dataset includes biases, twilight flats, and science frames. Link to the [Jupyter notebook](https://github.com/astro-datalab/notebooks-latest/tree/master/04_HowTos/DataReduction/DataReduction/DRAGONS_reduction_examples/GMOS_Imaging_Galaxy/GMOS_Imaging_Galaxy.ipynb).

### GNIRS_Imaging_GammaRayBurst.ipynb

GNIRS imaging (J-band - point source through keyhole) of GRB120116A. This is extracted from the Gemini/DRAGONS GNIRS tutorial, Example 1-B. Dataset includes flats, darks, and science frames. Link to the [Jupyter notebook](https://github.com/astro-datalab/notebooks-latest/tree/master/04_HowTos/DataReduction/DRAGONS_reduction_examples/GNIRS_Imaging_GammaRayBurst/GNIRS_Imaging_GammaRayBurst.ipynb).

### GSAOI_Imaging_EllipticalGalaxy.ipynb

GSAOI imaging (K-short) of a field around NGC5128. This is extracted from the Gemini/DRAGONS GSAOI tutorial, Section 3. Dataset includes flats, standard star, and science frames. Link to the [Jupyter notebook](https://github.com/astro-datalab/notebooks-latest/tree/master/04_HowTos/DataReduction/DRAGONS_reduction_examples/GSAOI_Imaging_EllipticalGalaxy/GSAOI_Imaging_EllipticalGalaxy.ipynb).

### NIRI_Imaging_Supernova.ipynb

NIRI imaging (K-prime) of SN2014J. This is extracted from the Gemini/DRAGONS NIRI tutorial, Section 4. Dataset includes flats, standard star, darks, and science frames. Link to the [Jupyter notebook](https://github.com/astro-datalab/notebooks-latest/tree/master/04_HowTos/DataReduction/DRAGONS_reduction_examples/NIRI_Imaging_Supernova/NIRI_Imaging_Supernova.ipynb).

### GMOS_Longslit_WhiteDwarf.ipynb
GMOS longslit data of a candidate DB white dwarf J2145+0031. This example is based on the Gemini DRAGONS GMOS longslit tutorial, Section 3. This dataset includes arcs, biases, and flats for the standard star and science target. This tutorial contains commented-out code that enables interactive mode. These cells are not required for this example, but feel free to uncomment them to test the interactive features. Link to the [Jupyter notebook](https://github.com/astro-datalab/notebooks-latest/tree/master/04_HowTos/DataReduction/DRAGONS_reduction_examples/GMOS_Longslit_WhiteDwarf/GMOS_Longslit_WhiteDwarf.ipynb).

---
## Additional resources

- **DRAGONS Documentation:** Data reduction with DRAGONS is currently available for all imaging data from current Gemini instruments. New modes and instruments are being added. Detailed instructions on the current stable version can be found on the [documentation page](https://dragons.readthedocs.io/en/stable/).

---
## Need help?

Do you have problems, comments, suggestions, or need help running DRAGONS? You can contact the US NGO members via our [Portal](http://ast.noao.edu/csdc/usngo).

If you are experiencing trouble with the DRAGONS (Py3.7) kernel or running the Jupyter Notebooks, please reach out to us via the Astro Data Lab Helpdesk: 
https://datalab.noirlab.edu/help/

---
