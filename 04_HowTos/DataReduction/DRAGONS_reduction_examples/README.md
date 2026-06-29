
## DRAGONS data reduction (Jupyter Notebooks):

This repository has Jupyter Notebook examples of data reduction for the Gemini Observatory instruments. Usually, you need DRAGONS installed on your computer to run these notebooks, but the Astro Data Lab has a custom kernel called DRAGONS-4.2.1 (DL,Py3.12) that will allow you to run them. The notebooks are set to open the DRAGONS kernel by default, but if this doesn’t happen, click on the name of the current kernel in the top-right corner and select DRAGONS-4.2.1 (DL,Py3.12). These notebooks were written using the [DRAGONS' Application Program Interface (API)](https://dragons.readthedocs.io/projects/recipe-system-users-manual/en/stable/appendices/full_api_example.html) for Python, based on the examples provided in the [DRAGONS Documentation](https://dragons.readthedocs.io/).

---
## Available Notebooks

### Imaging

| Notebook | Instrument | Target | Band | Description |
|---|---|---|---|---|
| [Flamingos2_Imaging_BrownDwarf.ipynb](Flamingos2_Imaging_BrownDwarf/Flamingos2_Imaging_BrownDwarf.ipynb) | Flamingos-2 | WISE J041358.14-475039.3 | Y | Brown dwarf imaging. Includes Flats, Darks, and Science frames. |
| [GMOS_Imaging_StarryField.ipynb](GMOS_Imaging_StarryField/GMOS_Imaging_StarryField.ipynb) | GMOS | Stellar field | i | Imaging of a stellar field. Includes Biases, Twilight Flats, and Science frames. |
| [GMOS_Imaging_Galaxy.ipynb](GMOS_Imaging_Galaxy/GMOS_Imaging_Galaxy.ipynb) | GMOS | NGC 5018 | g | Imaging of an elliptical galaxy. Includes Biases, Twilight Flats, and Science frames. |
| [GNIRS_Imaging_GammaRayBurst.ipynb](GNIRS_Imaging_GammaRayBurst/GNIRS_Imaging_GammaRayBurst.ipynb) | GNIRS | GRB120116A | J | Point-source imaging through the keyhole. Includes Flats, Darks, and Science frames. |
| [GSAOI_Imaging_EllipticalGalaxy.ipynb](GSAOI_Imaging_EllipticalGalaxy/GSAOI_Imaging_EllipticalGalaxy.ipynb) | GSAOI | NGC 5128 field | K-short | AO-assisted imaging of an elliptical galaxy field. Includes Flats, Standard Star, and Science frames. |
| [NIRI_Imaging_Supernova.ipynb](NIRI_Imaging_Supernova/NIRI_Imaging_Supernova.ipynb) | NIRI | SN 2014J | K-prime | Supernova imaging. Includes Flats, Standard Star, Darks, and Science frames. |

### Longslit Spectroscopy

| Notebook | Instrument | Target | Band | Description |
|---|---|---|---|---|
| [GMOS_Longslit_WhiteDwarf.ipynb](GMOS_longslit_WhiteDwarf/GMOS_Longslit_WhiteDwarf.ipynb) | GMOS | J2145+0031 | Optical | Candidate DB white dwarf longslit reduction. Includes arcs, biases, and flats for both the standard star and science target. Contains optional commented-out cells for interactive mode. |
| [GNIRS_longslit_Bestar.ipynb](GNIRS_longslit_Bestar/GNIRS_longslit_Bestar.ipynb) | GNIRS | HD 41335 | L | Longslit spectroscopy of a Be-star. Includes science, flats, tellurics, and a bad pixel mask. |
| [Flamingos2_longslit_JH_HK_point_source.ipynb](Flamingos2_longslit_JH_HK_point_source/Flamingos2_longslit_JH_HK_point_source.ipynb) | Flamingos-2 | J1344+0005 | JH HK | Longslit point-source spectroscopy. Includes science, flats, arcs, tellurics, and a bad pixel mask. |

### Cross-dispersed Spectroscopy

| Notebook | Instrument | Target | Band | Description |
|---|---|---|---|---|
| [GNIRS_XD_Short_Blue.ipynb](GNIRS_XD_Short_Blue/GNIRS_XD_Short_Blue.ipynb) | GNIRS | SN 2016ija / DLT16am | XD Short Blue | Cross-dispersed spectroscopy of a Type II supernova. Includes science, flats, arcs, pinholes, tellurics, and a bad pixel mask. |

### IFU Spectroscopy

| Notebook | Instrument | Target | Band | Description |
|---|---|---|---|---|
| [GHOST_IFU_Star.ipynb](GHOST_IFU_Star/GHOST_IFU_Star.ipynb) | GHOST | XX Oph | Optical | Standard-resolution, single-target IFU spectroscopy.Covers full reduction of a stellar target including processing for biases, flat fields, and arc calibration frames.  |


---
## Additional resources

- **DRAGONS Documentation:** Data reduction with DRAGONS is currently available for all imaging data from current Gemini instruments. New modes and instruments are being added. Detailed instructions on the current stable version can be found on the [documentation page](https://dragons.readthedocs.io/en/stable/).

---
## Need help?

Do you have problems, comments, suggestions, or need help running DRAGONS? You can contact the US NGO members via our [Portal](https://noirlab.edu/science/programs/csdc/usngo).

If you are experiencing trouble with the DRAGONS-4.2.1 (DL,Py3.12) kernel or running the Jupyter Notebooks, please reach out to us via the Astro Data Lab Helpdesk: 
https://datalab.noirlab.edu/help/

---
