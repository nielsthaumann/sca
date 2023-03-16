--- High-accuracy MEG/EEG diagnostics with SCA ---

The accuracy of electroencephalography (EEG) and magnetoencephalography (MEG) is challenged by overlapping neural sources. 
This lack of accuracy is a severe limitation to the application of MEG/EEG to clinical diagnostics. 
As a solution, we here introduce a spike density component analysis (SCA) method for isolating specific neural sources.

Spike density component analysis (SCA) decomposes EEG or MEG waveforms into neural components modeled with temporal probability density functions (use runsca.m).
Statistical thresholding of individual-level evoked responses can be applied to identify reliable evoked response SCA components (use runsca_stats.m).

For further information, please refer to: 

Haumann, N T; Hansen, B; Huotilainen, M; Vuust, P; Brattico, E;
"Applying Stochastic Spike train theory for high-accuracy human MEG/EEG",
Journal of Neuroscience Methods (2020), doi: https://doi.org/10.1016/j.jneumeth.2020.108743 

Bruzzone, S E P; Haumann, N T; Kliuchko, M; Vuust, P, Brattico, E;
"Applying Spike-density Component Analysis for high-accuracy auditory event-related potentials in children",
Clinical Neurophysiology (2021), https://doi.org/10.1016/j.clinph.2021.05.007

Haumann, N T; Petersen, B; Friis Andersen, A S; Faulkner, K S; Brattico, E; Vuust, P;
"Mismatch negativity as a marker of music perception in individual cochlear implant users: A spike density component analysis study",
Clinical Neurophysiology (2023), https://doi.org/10.1016/j.clinph.2023.01.015

This version of the SCA code runs in the Matlab environment and requires the FieldTrip toolbox to be installed (see https://github.com/fieldtrip). 
For EEGLab users it is possible to apply the 'eeglab2fieldtrip' function to analyze data processed in EEGLab (https://github.com/fieldtrip/fieldtrip/blob/master/external/eeglab/eeglab2fieldtrip.m)
'Runsca' also requires the Curve Fitting Toolbox to be installed for Gaussian or sine models. 
For experimental comparison with gamma models, the Statistics and Machine Learning Toolbox needs to be installed. 
