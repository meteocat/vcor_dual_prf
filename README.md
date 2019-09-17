[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/meteocat/vcor_dual_prf/master?filepath=.%2Fexamples%2F01_correct_and_display.ipynb)

Correction of dual-PRF velocity dealiasing errors
=================================================

This repository includes a function for correcting the dealiasing errors that arise in weather radar Doppler velocity data that have been collected using the dual-PRF technique. 

The function integrates four image processing methods developed in the literature, which apply field continuity for identification and correction of the errors. 

The function allows the user to tailor the dual-PRF error correction by specifying the neighbour kernel and selecting the statistic used for the estimation of the reference velocity. The correction procedures proposed in the literature may be reproduced through the particular selection of these statistics: the mean velocity (Joe and May, 2003), the median velocity (Holleman and Beekhuis, 2003), the circular mean of the PRF-scaled phase (Altube et al., 2017) and the circular mean of the phase (Hengstebeck et al., 2018).

[Usage example](https://mybinder.org/v2/gh/meteocat/vcor_dual_prf/master?filepath=.%2Fexamples%2F01_correct_and_display.ipynb)

## Detection/correction methods

- 'mean' : local mean velocity, [Joe and May, 2003](https://journals.ametsoc.org/doi/full/10.1175/1520-0426%282003%2920%3C429%3ACODPVE%3E2.0.CO%3B2)
- 'median' : local median velocity, [Holleman and Beekhuis, 2003](https://journals.ametsoc.org/doi/full/10.1175/1520-0426%282003%2920%3C443%3AAACODP%3E2.0.CO%3B2)
- 'cmean_sc' : local circular mean velocity, (with PRF-based scaling) [Altube et al., 2017](https://journals.ametsoc.org/doi/10.1175/JTECH-D-16-0065.1)
- 'cmean' : local circular mean velocity, [Hegstebeck et al., 2018](https://journals.ametsoc.org/doi/full/10.1175/JTECH-D-16-0230.1)


## Function parameters
```
correct_dualprf (radar, method_det, [...])
```
- **radar**: Py-ART radar instance

    Contains RAW file data and metadata. Must contain a dual-PRF radial velocity field.

- **method_det**: str

    Method used for detection of outliers.

- **[...]**:   
    
    
 Parameter | Class | Default value | Description
| ----------- | ----------- | ----------- | ----------- |
| vel_field | str | 'velocity' | Input velocity field name (dual-PRF) |
| kernel_det | array | None (7x7 array) | Neighbour kernel for local statistics (detection), 1/0 values |
| min_valid_det | int | 1 | Minimum number of valid neighbours (detection) |
