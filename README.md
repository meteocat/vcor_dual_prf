[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/meteocat/vcor_dual_prf/master?filepath=.%2Fexamples%2F01_correct_and_display.ipynb)

Correction of dual-PRF velocity dealiasing errors
=================================================

This repository includes a function for correcting the dealiasing errors that arise in weather radar Doppler velocity data that have been collected using the dual-PRF technique. 

The function integrates four image processing methods developed in the literature, all of which require spatial continuity by comparison of the gate velocity with a reference estimated from the velocities in the surrounding gates.

The function allows the user to tailor the dual-PRF error correction by specifying the neighbour kernel and selecting the statistic used for the estimation of the reference velocity. The correction procedures proposed in the literature may be reproduced through the particular selection of these statistics.

[TRY! Usage example](https://mybinder.org/v2/gh/meteocat/vcor_dual_prf/master?filepath=.%2Fexamples%2F01_correct_and_display.ipynb) (Takes a good while to load the environment, excuses for the inconvenience)

## Detection/correction methods

- 'mean' : local mean velocity, [Joe and May, 2003](https://journals.ametsoc.org/doi/full/10.1175/1520-0426%282003%2920%3C429%3ACODPVE%3E2.0.CO%3B2)
- 'median' : local median velocity, [Holleman and Beekhuis, 2003](https://journals.ametsoc.org/doi/full/10.1175/1520-0426%282003%2920%3C443%3AAACODP%3E2.0.CO%3B2)
- 'cmean_sc' : local circular mean velocity, (with PRF-based scaling) [Altube et al., 2017](https://journals.ametsoc.org/doi/10.1175/JTECH-D-16-0065.1)
- 'cmean' : local circular mean velocity, [Hegstebeck et al., 2018](https://journals.ametsoc.org/doi/full/10.1175/JTECH-D-16-0230.1)

## Common framework

The function is developed in Python language and builds upon the [Py-ART module](https://arm-doe.github.io/pyart/)<sup>[1]</sup>, as it takes a Py-ART radar object as input. The radar object provides a common framework for the application of the function to the wide range of data formats used by the growing community of Py-ART users.

<sup>[1]</sup><sub>Helmus, J.J. & Collis, S.M., (2016). The Python ARM Radar Toolkit (Py-ART), a Library for Working with Weather Radar Data in the Python Programming Language. Journal of Open Research Software. 4(1), p.e25. DOI: http://doi.org/10.5334/jors.119</sub>

## Positional function parameters
```
correct_dualprf (radar, method_det, ...)
```
- **radar**: Py-ART radar instance

    Contains RAW file data and metadata. Must contain a dual-PRF radial velocity field.

- **method_det**: str

    Method used for detection of outliers.

## Optional parameters/settings
    
 Parameter | Class | Default value | Description | Stage |
| ----------- | ----------- | ----------- | ----------- | ----------- |
| vel_field | str | 'velocity' | Input velocity field name (dual-PRF) | NA |
| kernel_det | array | None (7x7 array) | Neighbour kernel for local statistics, 1/0 values | detection |
| min_valid_det | int | 1 | Minimum number of valid neighbours | detection |
| max_dev | float | 1 | Maximum deviation fraction from 2xNyquist velocity | detection |
| two_step | bool | 1.0 | Separate detection and correction stages? (detected outliers are removed before correction) | NA |
| method_cor | str | same as *method_det* | Method used for correction of outliers. Only applied if *two_step*=True | correction |
| kernel_cor | array | same as *kernel_det* | Neighbour kernel for local statistics, 1/0 values | correction |
| min_valid_cor | int | 1 | Minimum number of valid neighbours |  correction |
| new_field | str | 'velocity_cor' | Corrected velocity field name in the output radar object | NA |
| new_field_lname | str | 'Outlier corrected dual-PRF velocity'' | Corrected velocity field LONG name in the output radar object | NA |

## Future/ work

- Implement a subroutine that identifies the PRF at which each radial has been scanned (see e.g. [Holleman and Beekhuis, 2003](https://journals.ametsoc.org/doi/full/10.1175/1520-0426%282003%2920%3C443%3AAACODP%3E2.0.CO%3B2))

- Test on different input data formats

- Speed median filtering with masked arrays?

