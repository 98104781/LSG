# LSG
Work in progress

This program is intended to enable the user to generate on-demand (mass)spectral libraries for the identification of lipid species.
Lipids are generated with respect to class and fatty acid composition, spectra are then generated with respect to their adducts with the help of a simple script.

Though there are some repositories available online which enable analysis of some common lipid species, such repositories are limited in various aspects:
* They offer a limited number of spectra/species
* Files are often large, including many species which are not of interest
* Spectra can be incomplete, lowering the confidence of identifications
* Others (eg, lipidmaps) may provide only structure files and requiring the user to generate, limiting their applicability


Fragmentation patterns used to generate the spectra are based on peer-reviewed studies.
Currently only supports Lyso PI, PI, PIP, PI2P under negative ESI mode.

User-interface is incomplete.
To specify the range of fatty-acids used to generate lipids, modification of the 'Tails_to_Generate' variable in 'GenerateLipids.py' is required.
The script is run via executing 'Main', selecting 'File' in the toolbar, followed by 'Export'.
A save location and name will be required.

The spectral library can be exported with the file extension '.MSP' selected.
Otherwise, an Excalibur compatible DDA precursor list (for DDA analysis via orbitrap) may be exported by selecting '.CSV'.
