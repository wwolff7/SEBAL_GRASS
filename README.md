# SEBAL_GRASS [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.167350.svg)](https://doi.org/10.5281/zenodo.167350)

<div data-badge-popover="right" data-badge-type="medium-donut" data-doi="10.5281/zenodo.167350" data-hide-no-mentions="true" class="altmetric-embed"></div>

## Requirements

* GRASS GIS 7.X [https://grasswiki.osgeo.org/wiki/GRASS-Wiki] (https://grasswiki.osgeo.org/wiki/GRASS-Wiki)
* Python 2.7.X [https://www.python.org/] (https://www.python.org/)

## Usage

### Organization of datasets

1. Download at [Earth Explorer] (http://earthexplorer.usgs.gov/) Landsat 8 scene (LS8 - OLI/TIRS)
    - For a scene without clouds, select the option **Level 1 GeoTIFF Data
    Product**

3. Reproject LS8 images for the coordinate system of interest
   - e.g, WGS 84 24N to SIRGAS 2000 24S, using gdal recursively in command line:
   - `mkdir rep && for i in *.TIF ; do gdalwarp -s_srs EPSG:32624 -t_srs EPSG:31984 -of GTiff $i rep\$i; done`

4. Remove null values (black borders) of LS8 images
   - e.g, using gdal recursively in command line:
   - `mkdir nodata && for i in *.TIF; do gdal_translate -a_nodata 0 $i nodata/$i; done`

5. Download at [Earth Explorer] (http://earthexplorer.usgs.gov/) the Digital Elevation Model (DEM) from ASTER (remember to choose the same region of LS8 scene)

6. It is necessary to reproject the DEM for the coordinate system of interest and rename to **MDT_Sebal.TIF**
   - e.g, WGS 84 24N to SIRGAS 2000 24S, using gdal in command line:
   - `gdalwarp -s_srs EPSG:32624 -t_srs EPSG:31984 -of GTiff ASTGTM2_S23W048.tif MDT_Sebal.TIF`

7. Launch a GRASS-GIS 7.X session
    - Select GRASS GIS database directory
    - Define a new GRASS location
      - Read the projection and datum terms from a georeferenced data file
      - Select the raster **MDT_Sebal.TIF**
    - Define a new GRASS mapset

8. Place **Sebal70.py** script in the directory where LS8 images are located

9. In Terminal Linux or Command Prompt Windows opened by GRASS, navigate to the directory where the **Sebal70.py** and LS8 images are located
   - Run python in command line:
     - `python Sebal70.py`
   - Follow the instrustructions indicated in Terminal
   - Use query tool to visualize cold and hot pixels in GRASS GIS display

## Remarks

* Tested on Ubuntu 14.04, 16.04 LTS, and Windows 8.

## References

ALLEN, R. G.; TASUMI, M.; TREZZA, R. Satellite-Based Energy Balance for
Mapping Evapotranspiration with Internalized Calibration (METRIC)-Model. **Journal of
Irrigation and Drainage Engineering**, v. 133, n. 4, p. 380–394, 2007. Available at:
<http://ascelibrary.org/doi/abs/10.1061/(ASCE)0733-9437(2007)133%3A4(380)>.

BASTIAANSSEN, W. G. M.; MENENTI, M.; FEDDES, R. A.; HOLTSLAG, A.
A. M. A remote sensing surface energy balance algorithm for land (SEBAL). 1.
Formulation. **Journal of Hydrology**, v. 212–213, n. 1–4, p. 198–212, 1998. Available at:
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

## Acknowledgements

To the grant 2016/15342-2, São Paulo Research Foundation (FAPESP) by the financial support.
