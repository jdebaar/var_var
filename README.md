+++++++++++++++++++++++++++++++++++++++++++++

Manuscript: The elephant in the room: Improved kriging posterior variance, applied to
spatial climate data

Code: Method illustration for Kriging with non-Euclidian distance and multi-objective hyperparameter estimation

Source: https://github.com/jdebaar/var_var

Date: 18 March 2022

Contact: jouke.de.baar@knmi.nl

+++++++++++++++++++++++++++++++++++++++++++++

1 Purpose of the code

In our manuscript, we present methods and results for kriging with non-Euclidian distance and multi-objective hyperparameter estimation. The purpose of this code is to illustrate the methods proposed.

2 Running the code

The code can simply be executed by running 'main.R' in R

3 Directories

The following directories are available:
- input_data: The data required to run the illustration
- output_figures: The resulting figures
- licences: The licens files 

4 Contents of the code

4.1 Contents of 'main.R'

This script is fully annotated, including references to the pseudo-algorithm in the manuscript. The script contains the following steps:
- Settings, with a description.
- Loading the data
- Sampling the synthetic station data
- Pre-computing the lags
- Brute force single-objective optimization (reference method, Algorithm 1)
- Brute force multi-objective optimizatin (proposed methdod, Algorithm 2)

4.2 Contents of 'krigFunctions.R'

This file contains the following functions:
- corrMatrix: Compute correlation matrix
- corrLag: Compute non-Euclidian lags, to be used for correlation matrix
- corrMatrixFast: More efficient implementation of corrMatrix
- krigPost: Compute the kriging posterior
- krigLOOCVfast: Leave-one-out cross validation
