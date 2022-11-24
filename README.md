# Modelling novelty detection in the thalamocortical loop
Here we provide the Matlab code to reproduce the results of the paper "Modelling novelty detection in the thalamocortical loop".

1. download this repository into your local directory <local_path> and unzip it.
2. create two folders named  `Simulation Results` and `Figures` in the directory <local_path>/SomatosensorySSA_Model-master to save the generated data and figures.

## Run our model on different protocols
The stimulus-specific adaptation and true deviance detection results can be reproduced by running `main.m` with `Cond_Arr = [1 4]` specified in the script, which represents the oddball low and many-standards protocols, respectively.

Similarly, the novelty-predicting effect results can be reproduced by running `main.m` with `Cond_Arr = [7 8]`, which represents the sequenced and randomized protocols, respectively.

The `main.m` script in turn runs `model_init.m` that is for initializing our model and stimulus parameters and `model_run.m` that is for simulating our model on each protocol. All simulation results are saved in `Simulation Results` folder.

## Visualize simulation results
We provide the `Visualization.m` to reproduce most figures in the papers. The `Visualization.m` script loads the results saved `Simulation Results` folder and generates the figures saved in `Figures` folder. 

## License
The provided Matlab code is licenced under the MIT license, see the LICENSE file.
