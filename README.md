Spike density component analysis decomposes EEG or MEG waveforms into neural components 
modelled with temporal probability density functions.

For further information, please refer to: 

Haumann, N T; Hansen, B; Huotilainen, M; Vuust, P; Brattico, E;
"Applying Stochastic Spike train theory for high-accuracy human MEG/EEG",
Journal of Neuroscience Methods (2020), doi: https://doi.org/10.1016/j.jneumeth.2020.108743 

This version of the SCA code runs in the Matlab environment and requires the FieldTrip toolbox to be installed (see https://github.com/fieldtrip). 
For EEGLab users it is possible to apply the 'eeglab2fieldtrip' function to analyze data processed in EEGLab (https://github.com/fieldtrip/fieldtrip/blob/master/external/eeglab/eeglab2fieldtrip.m)
'Runsca' also requires the Curve Fitting Toolbox to be installed for Gaussian or sine models. 
For experimental comparison with gamma models, the Statistics and Machine Learning Toolbox needs to be installed. 
