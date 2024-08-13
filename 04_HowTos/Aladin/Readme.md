## Aladin Lite test
This repository contains my notes and test notebooks for Aladin Lite v3. 

----

##Aladin Lite v3
- Download latest Aladin Lite v3 package
  - https://aladin.cds.unistra.fr/AladinLite/doc/release/
- I found a tutorial on Aladin’s website on ‘Building an interactive sky map with Aladin Lite’
  - https://aladin.cds.unistra.fr/AladinLite/doc/tutorials/interactive-finding-chart/

In order to run the Aladin Lite viewer on your machine, you will need to create a blank file called index.html. To initialize the viewer, add the following lines to index.html.

<!DOCTYPE>
<html>
  <body>
    <h1>Trifid interactive map</h1>
	<!-- our code needs jQuery library -->
	<script type="text/javascript" src="http://code.jquery.com/jquery-1.9.1.min.js" charset="utf-8"></script>

	<!-- Aladin Lite container at requested dimensions -->
	<div id="aladin-lite-div" style="width:700px;height:400px;"></div>

	<!-- Aladin Lite JS code -->
	<script type="text/javascript" src="https://aladin.cds.unistra.fr/AladinLite/api/v3/latest/aladin.js" charset="utf-8"></script>

	<!-- Creation of Aladin Lite instance with initial parameters -->
	<script type="text/javascript">
	    let aladin;
            A.init.then(() => {
	        aladin = A.aladin('#aladin-lite-div', {survey: "P/DSS2/color", fov:1.5, target: "trifid nebula"});
            });
	</script>
  </body>
</html>

Then run 'python -m SimpleHTTPServer' in the same directory as index.html. Finally, enter this url http://0.0.0.0:8000/index.html into your browser to access the interactive Aladin viewer.

This example is very basic. The version of index.html that is stored in this repo is the test file I was using that adds HiPS and MOCs to the viewer. 

----

## Aladin Multi-Order Coverage map (MOC) test
- MOCPy v0.15.0
  - https://cds-astro.github.io/mocpy/_collections/notebooks/filtering_astropy_table.html

Create_MOC_test.ipynb:
    - This notebook will take in a catalog of sources (must contain ra and dec) and create a MOC.
    - For this example, I am using a small sampling of the des_dr2 catalog which can be found in this repo as des_dr2.csv.
    - Towards the end of the notebook, ipyaladin will be used to bring up the Aladin viewer within the notebook and display the newly created MOC. 

NOTE: Before using MOCPy on a DataLab notebook, you must add the following line to the top:
    !pip install ipyaladin mocpy

This notebook still needs documentation, but is fully working and is essentially complete. 

----

##Aladin Hierarchical Progressive Survey Catalogue (HiPS) test
- Use this link to download HiPSgen-cat
  - https://aladin.cds.unistra.fr/hips/Hipsgen-cat.jar
- The HiPS catalog tool has a small tutorial page here
  - https://aladin.cds.unistra.fr/hips/HipsCat.gml

To generate a HiPS cat from a catalog that contains ra and dec (all values must be in decimal degrees and cannot contain empty values), use the following command in the same directory as Hipsgen-cat.jar and the catalog:

java -jar Hipsgen-cat.jar -cat <catalog name, no extension> -in <catalog name, with extension> -f ASCII -af CSV -out <desired output directory name> -ra ra -dec dec -score <scored value>

    - <catalog name, no extension> means if the files is called catalog_test.csv, you would enter catalog_test
    - <catalog name, with extension> means you should enter the file name with its extension, so catalog_test.csv
    - <desired output directory name>, this command will output a directory containing many files. What you enter here will eb the name of the parent directory that contains all of the files.
    - <scored value> can be a different column in the input catalog that can be used to rank the sources. For example, if I were to use 'mag_auto_g', the produced HiPS cat will prioritize showing the sources with bright g-band magnitudes first. Then as you zoom in, the next brightest sources will be displayed.

As of 08/05/24, I am able to sucessfully create HiPS catalogs, and display them on my Aladin Lite viewer and the desktop version of Aladin, but I am unable to display it in a Jupyter Notebook becasue ipyaladin is not compatable with HiPS cats yet. According to the Issues page on the ipyaladin GitHub, they are working on including this functionality, but it is not available yet.
