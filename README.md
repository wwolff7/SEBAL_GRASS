# SEBAL_GRASS

## Requirements

* GRASS GIS 7.X [https://grasswiki.osgeo.org/wiki/GRASS-Wiki] (https://grasswiki.osgeo.org/wiki/GRASS-Wiki)
* Python 2.7.X [https://www.python.org/] (https://www.python.org/)

## Usage

### Organization of datasets 

1. Make the download of Landsat 8 scene (LS8 OLI/TIRS) at [Earth Explorer] (http://earthexplorer.usgs.gov/)
    - For a scene without clouds, select the option **Level 1 GeoTIFF Data
    Product**

3. Reproject LS8 images for the coordinate system of interest
   - e.g, WGS 84 24N to SIRGAS 2000 24S, using gdal recursively in command line:
   - `mkdir rep && for i in *.TIF ; do gdalwarp -s_srs EPSG:32624 -t_srs EPSG:31984 -of GTiff $i rep\$i; done` 

4. Remove the null values (black borders) of LS8 images 
   - e.g, using gdal recursively in command line:
   - `mkdir nodata && for i in *.TIF; do gdal_translate -a_nodata 0 $i nodata/$i; done` 

5. For the same region of LS8 scene make the download of the ASTER GLOBAL DEM at [Earth Explorer] (http://earthexplorer.usgs.gov/)
   
6. Reproject the DEM for coordinate system of interest and rename to **MDT_Sebal.TIF**
   - e.g, WGS 84 24N to SIRGAS 2000 24S, using gdal in command line:
   - `gdalwarp -s_srs EPSG:32624 -t_srs EPSG:31984 -of GTiff ASTGTM2_S23W048.tif MDT_Sebal.TIF` 

7. Launch a GRASS-GIS 7.X session
    - Select GRASS GIS database directory
    - Define a new GRASS location 
      - Read a projection and datum terms from a georeferenced data file
      - Select the raster **MDT_Sebal.TIF** 
    - Define a new GRASS mapset  

8. Place **Sebal70.py** script in directory where is located LS8 images

9. In Terminal Linux navigate into the directory where is located the **Sebal70.py** and the images LS8
   - Run python in command line:
   - `python Sebal70.py`
   - Follow the instrustructions indicated in Terminal
   - use query tool to visualise cold and hot pixels in GRASS GIS display

# Remarks

* Tested on Ubuntu 14.04 and 16.04 LTS 

# References

ALLEN, R. G.; TASUMI, M.; TREZZA, R. Satellite-Based Energy Balance for
Mapping Evapotranspiration with Internalized Calibration (METRIC)-Model. **Journal of
Irrigation and Drainage Engineering**, v. 133, n. 4, p. 380–394, 2007. Available at:
<http://ascelibrary.org/doi/abs/10.1061/(ASCE)0733-9437(2007)133%3A4(380)>.

BASTIAANSSEN, W. G. M.; MENENTI, M.; FEDDES, R. A.; HOLTSLAG, A.
A. M. A remote sensing surface energy balance algorithm for land (SEBAL). 1.
Formulation.**Journal of Hydrology**, v. 212–213, n. 1–4, p. 198–212, 1998. Available at:
<http://linkinghub.elsevier.com/retrieve/pii/S0022169498002534>.

SILVA, B. B.; BRAGA, A. C.; BRAGA, C. C.; OLIVEIRA, L. M. M. D.; MONTENEGRO,
S. M. G. L.; BARBOSA JR., B. Procedures for calculation of the albedo with
OLI-Landsat 8 images: Application to the Brazilian semi-arid. **Revista Brasileira
de Engenharia Agrícola e Ambiental**, v. 20, n. 1, p. 3–8, 2016. Available at:
<http://dx.doi.org/10.1590/1807-1929/agriambi.v20n1p3-8>.

TASUMI, M.; ALLEN, R. G.; TREZZA, R.; WRIGHT, J. L. Satellite-Based Energy
Balance to Assess Within-Population Variance of Crop Coefficient Curves. **Journal of
Irrigation and Drainage Engineering**, v. 131, n. 1, p. 94–109, 2005. Available at:
<http://ascelibrary.org/doi/10.1061/(ASCE)0733-9437(2005)131:1(94)>.





