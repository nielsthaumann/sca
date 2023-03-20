<h1>High-accuracy MEG/EEG diagnostics with SCA </h1>

The accuracy of electroencephalography (EEG) and magnetoencephalography (MEG) is challenged by overlapping neural sources. 
This lack of accuracy is a severe limitation to the application of MEG/EEG to clinical diagnostics. 
As a solution, we here introduce a spike density component analysis (SCA) method for isolating specific neural sources.

Spike density component analysis (SCA) decomposes EEG or MEG waveforms into neural components modeled with temporal probability density functions (with the <a href="https://github.com/nielsthaumann/sca/blob/master/runsca.m">runsca.m</a> function). 
Statistical thresholding of individual-level evoked responses can be applied to identify reliable evoked response SCA components (with the <a href="https://github.com/nielsthaumann/sca/blob/master/runsca_stats.m">runsca_stats.m</a> function).

</br>

<p align="center">
  <image width="510" height="1110" src="https://repository-images.githubusercontent.com/251616840/e1517d29-b7a4-40c2-8019-6f40f0fc85b1">
</p>
<i>Overview of the steps in the spike density component analysis (SCA) statistics procedure. (A) Single-participant EEG difference waveforms decomposed into neural activity with independent scalp distributions (left) and time-courses (right). (B) Constraining the peak electrodes, polarity, and latency reveals candidate MMN components. (C) The spatiotemporal filter (black color) efficiently focuses the across-trial EEG variance on an MMN candidate component to be statistically tested (interquartile range shown in shaded color; comparison to no filtering shown in blue color). (D) Across-trial amplitude histograms. (E) An automatic detection of an SCA component clearly representative of individual-level MMN (shown with shaded 99% confidence intervals uncorrected for multiple testing).</i>
  (https://doi.org/10.1016/j.clinph.2023.01.015)
</br></br>

For further information, please refer to: 

Haumann, N T; Hansen, B; Huotilainen, M; Vuust, P; Brattico, E;
<b>"Applying Stochastic Spike train theory for high-accuracy human MEG/EEG"</b>,
Journal of Neuroscience Methods (2020), https://doi.org/10.1016/j.jneumeth.2020.108743 

Bruzzone, S E P; Haumann, N T; Kliuchko, M; Vuust, P, Brattico, E;
<b>"Applying Spike-density Component Analysis for high-accuracy auditory event-related potentials in children"</b>,
Clinical Neurophysiology (2021), https://doi.org/10.1016/j.clinph.2021.05.007

Haumann, N T; Petersen, B; Friis Andersen, A S; Faulkner, K S; Brattico, E; Vuust, P;
<b>"Mismatch negativity as a marker of music perception in individual cochlear implant users: A spike density component analysis study"</b>,
Clinical Neurophysiology (2023), https://doi.org/10.1016/j.clinph.2023.01.015

</br>

This version of the SCA code runs in the Matlab environment and requires the FieldTrip toolbox to be installed (see https://github.com/fieldtrip). 
For EEGLab users it is possible to apply the 'eeglab2fieldtrip' function to analyze data processed in EEGLab (https://github.com/fieldtrip/fieldtrip/blob/master/external/eeglab/eeglab2fieldtrip.m)
'Runsca' also requires the Curve Fitting Toolbox to be installed for Gaussian or sine models. 
For experimental comparison with gamma models, the Statistics and Machine Learning Toolbox needs to be installed. 
