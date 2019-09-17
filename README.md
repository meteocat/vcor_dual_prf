[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/meteocat/vcor_dual_prf/master?filepath=.%2Fexamples%2F01_correct_and_display.ipynb)

Correction of dual-PRF velocity dealiasing errors
=================================================

This repository includes a function for correcting the dealiasing errors that arise in weather radar Doppler velocity data that have been collected using the dual-PRF technique. 

The function integrates four image processing methods developed in the literature, which apply field continuity for identification and correction of the errors. 

The function allows the user to tailor the dual-PRF error correction by specifying the neighbour kernel and selecting the statistic used for the estimation of the reference velocity. The correction procedures proposed in the literature may be reproduced through the particular selection of these statistics: the mean velocity (Joe and May, 2003), the median velocity (Holleman and Beekhuis, 2003), the circular mean of the PRF-scaled phase (Altube et al., 2017) and the circular mean of the phase (Hengstebeck et al., 2018).

[Try usage example](https://mybinder.org/v2/gh/meteocat/vcor_dual_prf/master?filepath=.%2Fexamples%2F01_correct_and_display.ipynb)
