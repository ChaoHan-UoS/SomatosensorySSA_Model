# Modelling novelty detection in the thalamocortical loop
Here we provide the Matlab code to reproduce the results of the paper "Modelling novelty detection in the thalamocortical loop".

## Run our model on different protocols
The stimulus-specific adaptation and true deviance detection results can be reproduced by running `main.m` with `Cond_Arr = [1 4]` specified in the script, which represents the oddball low and many-standards protocols.

Similarly, the novelty-predicting effect results can be reproduced by running `main.m` with `Cond_Arr = [7 8]` specified in the script, which represents the sequenced and randomized protocols.

The `main.m` script in turn runs `model_init.m` that is for initializing our model and stimulus parameters and `model_run.m` that is for simulating our model on each protocol. All simulation results are saved in `/Simulation Results` folder that is needed to be created by the user. 

## Visualize simulation results
We provide the `Visualization.m` to reproduce most figures in the papers. The `Visualization.m` script loads the results saved `/Simulation Results` folder and generates the figures saved in `/Figures` folder created by the user. 

## Demo for generating Figure S1 in the paper
1. download this repository into your local directory <local_path>
2. download the simulation results of oddball low and many-standards protocols: 
https://figshare.com/articles/dataset/Simulation_Results/21614145 and unzip it in the directory <local_path>/SomatosensorySSA_Model-master
3. run `Plotting_S1.m` to generate Figure S1 in the paper

## License
The provided Matlab code is licenced under the MIT license, see the LICENSE file.
