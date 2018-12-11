README file for the DMS_P_SAT database and related data analysis

Martí Galí Tàpias, 2018-12-11

Questions and requests can be addressed to:
marti.gali.tapias@gmail.com

Content and data sources
========================
This curated dataset is based on the public global sea-surface dimethylsulfide (DMS) database. The original data was downloaded from http://saga.pmel.noaa.gov/dms/ on 14 Nov 2013. PMEL data contributions are described in “ContribNum_DMSdb_all.docx”. The PMEL database was then enlarged with more recent data as described in detail in the supporting online material of Galí et al. 2015. The in situ dataset was further enhanced by addition of climatological data (MLD, nutrients…), phytoplankton size classes, biogeochemical provinces, and satellite matchup data (SeaWiFS, MODIS-Aqua and AVHRR). Finally, a number of quality control flags were produced for the in situ data. Detailed information on all these procedures can be found in Galí et al. 2015 (especially the supporting online materials).


References
==========
Galí, M., Devred, E., Levasseur, M., Royer, S. J., & Babin, M. (2015). A remote sensing algorithm for planktonic dimethylsulfoniopropionate (DMSP) and an analysis of global patterns. Remote Sensing of Environment, 171, 171-184. https://doi.org/10.1016/j.rse.2015.10.012

Galí, M., Levasseur, M., Devred, E., Simó, R., & Babin, M. (2018). Sea-surface dimethylsulfide (DMS) concentration from satellite data at global and regional scales.  Biogeosciences, 15(11), 3497-3519. https://doi.org/10.5194/bg-15-3497-2018


Data files and code
===================
The curated dataset (43080 rows, 103 columns) is stored in matlab *.mat format and in *.csv format:

standard_QCed_DMS_P_db_11-Dec-2018.mat		data and headers
standard_QCed_DMS_P_db_11-Dec-2018.csv		data
headers_simple.txt				headers (variable names)
headers_explained.txt				detailed description of variables

Please check carefully headers_explained.txt for a more detailed description of variables and their units.

I produced a matlab *.m file with lots of tips for data analysis:
data_analysis_tips.m

Following these database tweaks you should be able to reproduce the results found in the two papers mentioned above. Some of them took some (or a long) time to figure out! SO this should be a very useful resource.


How to cite
===========
This dataset can be freely distributed, but please cite it using its DOI on Zenodo:
https://doi.org/10.5281/zenodo.2205131

When relevant, please cite also the journal articles mentioned above.


Acknowledgments
================
I acknowledge NASA’s Ocean Biology Processing Group (OBPG) and PMEL for maintaining the respective databases and making data freely available.


Note
====
By publishing this curated dataset, I want by no means to appropriate myself of the datasets obtained by other researchers and research programs. Distributing this database with its own DOI is a contribution to to the traceability of my past research, and to future research efforts.


Enjoy!
======
Martí Galí Tàpias
