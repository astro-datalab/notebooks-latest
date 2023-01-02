*Version:* 20221121 | *Author:* Giada Pastorelli <gpastorelli.astro@gmail.com>

# LSST simulations with TRILEGAL: Jupyter notebooks


## Description

Series of tutorial contributed notebooks to explore the LSST simulations performed with the TRILEGAL code (Dal Tio et al., 2021).
The notebooks are included in the directory `05_Contrib/LSST_Sim`.

For a detailed description of the simulations, please refer to [Dal Tio et al., (2021)](https://doi.org/10.3847/1538-4365/ac7be6).

Below is a summary of the main feautures of the available synthetic catalogs. 

The synthetic catalogs are produced with the TRILEGAL code that generates mock stellar samples through a population synthesis approach. The catalogs include all Milky Way stars that are expected to be visible on the ﬁnal stacked LSST images, and are complete down to the r = 27.5 mag depth of the co-added Wide–Fast–Deep survey images.
The Magellanic Clouds are also included. 

We provide a complete simulation for single stars (`lsst_sim.simdr2` dataset), and a second simulation for a fraction of the binaries, including the interacting ones, as derived with the BinaPSE module of TRILEGAL (`lsst_sim.simdr2_binary` dataset).

This project is part of the in-kind LSST contribution for the SMWLV Science Collaboration *Population models of the LSST stellar content*, PI Dr. L. Girardi, INAF-Padova (Italy).

## Content of the simulations

### Stellar and population properties

1. Photometry in LSST u, g, r, i, z, y AB magnitudes and Gaia G, GBP, GRP Vega magnitudes;
2. Stellar properties, such as luminosity, temperature, surface composition, gravity, pulsation periods;
3. Positional and kinematic properties, such as distances, proper motions, space velocities;
4. Additional quantities for binary systems, i.e stellar types, orbital parameters, radial velocity amplitudes, maximum depth of primary and secondary eclipses

### Evolutionary phases and specific phenomena

TRILEGAL simulations include the following evolutionary phases for single stars: Pre-main sequence (PMS), Main sequence (MS), Hertzsprung gap, Red giant branch (RGB), Core helium burning (CHeB, Early asymptotic giant branch (EAGB), Thermally pulsing AGB (TPAGB), Post-AGB, CO-WD. Products of binary evolution are also included in the binary catalogs: Helium main sequence, Helium Hertzsprung gap, Helium giant branch, He-WD, ONe-WD, Neutron star, Black hole.

The current simulated datasets include detached and interacting binaries, Cepheids, Long Period Variables, and eclipsing binaries. 

## Access

The synthetic catalogs are accessible through the NOIRLab Astro Data Lab:

1. Single star catalog: <https://datalab.noirlab.edu/query.php?name=lsst_sim.simdr2>
2. Binary star catalog: <https://datalab.noirlab.edu/query.php?name=lsst_sim.simdr2_binary>

Stellar density maps are incorporated in the LSST metrics analysis framework (MAF).

Please refer to the series of notebooks included here for detailed instructions and recommendations about the use of these datasets. 

## Feedback and contacts 

We welcome your feedback about the available notebooks, as well as requests for specific tutorials and examples. 

Main contact:
Giada Pastorelli <gpastorelli.astro@gmail.com>

PI of the project:
Leo Girardi <leo.girardi@inaf.it>
