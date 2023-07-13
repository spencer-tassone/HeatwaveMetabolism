Dataset title: Metabolism estimates for 356 U.S. rivers (2007-2017)
Source: USGS Science Data Catalog
Publication date: 2018-10-15
URL: https://www.sciencebase.gov/catalog/item/59bff507e4b091459a5e0982
DOI: 10.5066/F70864KX

Citation:
Appling, A.P., Read, J.S., Winslow, L.A., Arroita, M., Bernhardt, E.S., Griffiths, N.A., Hall, R.O., Jr., Harvey, J.W., Heffernan, J.B., Stanley, E.H., Stets, E.G., and Yackulic, C.B., 2018, Metabolism estimates for 356 U.S. rivers (2007-2017): U.S. Geological Survey data release, https://doi.org/10.5066/F70864KX. 

Summary:
This data release provides modeled estimates of gross primary productivity, ecosystem respiration, and gas exchange coefficients for 356 streams and rivers across the United States. The release also includes the model input data and alternative input data, model fit and diagnostic information, spatial data for the modeled sites (catchment boundaries and site point locations), and potential predictors of metabolism such as discharge and light availability.

The data are organized into these items:

    Site data - Site identifiers, details, and quality indicators - Table with 1 row per site (tab-delimited file)
    Spatial data
        a. Site coordinates - One shapefile of points for all sites combined (.shp, .shx, .dbf, and .prj files)
        b. Site catchment boundaries - One shapefile of polygons for all sites combined (.shp, .shx, .dbf, and .prj files)
    Timeseries data - Data on water quality and quantity, collected or computed from outside sources - Tables with one row per time series observation (1 tab-delimited file per site-variable combination, 1 zip file per site)
    Model inputs - Data formatted for use in estimating metabolism - Tables of prepared time series inputs (1 tab-delimited file per site, in 1 zip file per site)
    Model configurations - Model specifications used to estimate metabolism - Table with 1 row per model (1 tab-delimited file, compressed into zip file)
    Model outputs - Complete fits from metabolism estimation models - Text and 4 tables for each model (tab-delimited files, 1 zip file per model)
    Model diagnostics - Key diagnostics and overall assessments of model performance - Table with 1 row per model (1 tab-delimited file, compressed into zip file)
    Metabolism estimates and predictors - Daily metabolism estimates and potential predictor variables to support further exploration - Table with 1 row per site-date combination (1 tab-delimited file, compressed into zip file)

This work was funded by the USGS Powell Center (Working Group title: "Continental-scale overview of stream primary productivity, its links to water quality, and consequences for aquatic carbon biogeochemistry"). Additional financial support came from the USGS NAWQA program and Office of Water Information. NSF grants DEB-1146283 and EF1442501 partially supported ROH. The USGS Wisconsin Modeling Center provided access to the HTCondor cluster for generating metabolism estimates.

Associated manuscript:
Appling, A.P., Read, J.S., Winslow, L.A., Arroita, M., Bernhardt, E.S., Griffiths, N.A., Hall, R.O., Jr., Harvey, J.W., Heffernan, J.B., Stanley, E.H., Stets, E.G., and Yackulic, C.B., 2018, The metabolic regimes of 356 rivers in the United States. Scientific Data, Volume 5, Article Number 180292. https://doi.org/10.1038/sdata.2018.292. 

------------------------------------
site_data.tsv

Column		Description

site_name	USGS site number as "nwis_<site number>"
nwis_id		USGS site number, leading 0 may be removed
long_name	USGS descriptive site name
lat		Decimal latitude
lon		Decimal longitude
coord_datum	GPS datum
alt		Site altitude (feet above sea level?)
alt_datum	Site altitude datum
site_type	USGS site type
dvqcoefs.c	hydraulic geometry coefficients describing the at-a-station relationships among river discharge,
dvqcoefs.f		velocity, width, and depth. These were obtained from an analysis by Gomez-Velez et al. 
dvqcoefs.a		2015 (see Source Citation) and are regionalized at the HUC2 (USGS Hydrologic Unit Code 2) 
dvqcoefs.b		level based on several thousand measurements of instantaneous low-flow and bankful depths 
dvqcoefs.k		and widths in the conterminous United States.
dvqcoefs.m
struct.canal_flag	Probability flag for whether site is downstream of a canal, where values are either 95, 80, 50, or 0, where 95 indicates 
			the distance to structure is within the 95th percentile of O2 turnover distance, and has the least probable interference
			on metabolism estimates. 0 indicates the most probable interference on metabolism estimates.
struct.dam_flag		"	"	"	""	"	"	a dam
struct.npdes_flag   	"	"	"	""	"	"	a point-source discharge point



