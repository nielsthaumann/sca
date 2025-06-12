<h1>High-accuracy MEG/EEG diagnostics with SCA </h1>

The accuracy of electroencephalography (EEG) and magnetoencephalography (MEG) is challenged by overlapping neural sources. 
This lack of accuracy is a severe limitation to the application of MEG/EEG to clinical diagnostics. 
As a solution, we here introduce a spike density component analysis (SCA) method for isolating specific neural sources.

Spike density component analysis (SCA) decomposes EEG or MEG waveforms into neural components modeled with temporal probability density functions (with the <a href="https://github.com/nielsthaumann/sca/blob/master/runsca.m">runsca.m</a> function). </br>
Statistical thresholding of individual-level evoked responses can be applied to identify reliable evoked response SCA components (with the <a href="https://github.com/nielsthaumann/sca/blob/master/runsca_stats.m">runsca_stats.m</a> function). </br>
Statistical thresholding of individual-level evoked responses can also be performed with a more liberal half-split average consistency (SCA-HSAC) procedure [runsca_hsac](https://github.com/nielsthaumann/sca/blob/master/runsca_hsac.m) or a simple template match (SCA-TM) approach [runsca_tm](https://github.com/nielsthaumann/sca/blob/master/runsca_tm.m) . 

</br>

<p align="center">
  <image width="510" height="1110" src="https://repository-images.githubusercontent.com/251616840/e1517d29-b7a4-40c2-8019-6f40f0fc85b1">
</p>
<i>Overview of the steps in the spike density component analysis (SCA) statistics procedure. (A) Single-participant EEG difference waveforms decomposed into neural activity with independent scalp distributions (left) and time-courses (right). (B) Constraining the peak electrodes, polarity, and latency reveals candidate MMN components. (C) The spatiotemporal filter (black color) efficiently focuses the across-trial EEG variance on an MMN candidate component to be statistically tested (interquartile range shown in shaded color; comparison to no filtering shown in blue color). (D) Across-trial amplitude histograms. (E) An automatic detection of an SCA component clearly representative of individual-level MMN (shown with shaded 99% confidence intervals uncorrected for multiple testing).</i>
</br></br>
<p align="center">
  <a href=https://doi.org/10.1016/j.clinph.2023.01.015> (Haumann et al., 2023, "Mismatch negativity as a marker of music perception in individual cochlear implant users: <br>A spike density component analysis study",<i> Clinical Neurophysiology</i>)</br>
  </a>
</br></br>


# Explanation of the SCA algorithm

Imagine the human brain as a megacity with ten million people all signaling each other messages at once. EEG electrodes and MEG sensors would be like microphones set up outside the city walls, used to pick up a message from just one of the ten million people:
<p align="center">
  <image width="563" height="353" src="https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs40069-018-0319-7/MediaObjects/40069_2018_319_Fig1_HTML.png">
</p>
<p align="center">
  <a href=https://doi.org/10.1186/s40069-018-0319-7> (Hwang et al., 2019. "Experimental verification of modal identification of a high-rise building using independent component analysis." <i>  International Journal of Concrete Structures and Materials</i>, 13, 1-15)</br>
  </a>
</p>
</br>
Suppose some news breaks in the megacity. The exact time news reaches an individual is uncertain, but the arrival time varies systematically around a predicted time. </br>
    Similarly, the images below show white flashes at variable times of arrival during the influx of Ca2+ ions after 'neuronal spiking' potentials have pushed out Na+ ions of neuron cell bodies as part of the electrochemical signaling cascade in the brain:
<br></br>

![CaImAn Calcium Imaging Analysis 3-4](https://github.com/user-attachments/assets/0851bfa7-91fe-4fef-aa96-61e5eebda7a3)
<p>
  <a href=https://youtu.be/hs0TIO-8NWQ?si=rtYYkDvMJF_Oi1Sw> (Jérémie Kalfon. "A Computational toolbox for large scale Calcium Imaging Analysis* and behavioral analysis" </br>
  </a>
</p>
<br></br>
'Neuronal spiking' with systematic noise and influx of Ca2+ ions are a part of the electrochemical signaling cascade that triggers post-synaptic 'local field potentials'. The post-synaptic 'local field potentials' from a large number (>10,000) of neurons are the main contributors to the brain activity measured with MEG/EEG: 
<br></br>
<p align="left">
  <image width="563" height="482" src="https://upload.wikimedia.org/wikipedia/commons/8/8d/SimulationNeuralOscillations.png">
   </br>
  <a href=https://en.wikipedia.org/wiki/Neural_oscillation#/media/File:SimulationNeuralOscillations.png> (TjeerdB. "Neuronal spiking (firing of neurons) is simulated by a rate-modulated Poisson process (upper panel). Local field potential is simulated by the low-pass filtered sum of a number of these processes, representing the mean activity of a large number of neurons (lower panel)." </br>
</a>
</p>
<br></br>
A popular assumption is that brain activity measured with MEG/EEG is shaped as ideal sine waves. However, in line with the observations of systematic noise in the timing of 'neuronal spiking' potentials, it turns out that the large-scale neuronal electromagnetic activity measured with MEG/EEG is accurately modeled by temporal Gaussian density functions with systematic noise. The SCA modeling with temporal Gaussian density functions enables accurate analysis of the timing of evoked responses in the brain (ER) and suppression of ongoing interfering brain activity: 
<br></br>
<p align="center">
  <image width="563" height="819" src="https://ars.els-cdn.com/content/image/1-s2.0-S0165027020301667-gr1_lrg.jpg">
</p>
<br></br>
<p align="center">
  <image width="750" height="772" src="https://ars.els-cdn.com/content/image/1-s2.0-S0165027020301667-gr3.jpg">
</p>
<br></br>
<p align="center">
  <image width="750" height="772" src="https://ars.els-cdn.com/content/image/1-s2.0-S0165027020301667-gr10_lrg.jpg">
</p>
<p align="center">
  <a href=https://doi.org/10.1016/j.jneumeth.2020.108743> (Haumann et al., 2020, "Applying Stochastic Spike train theory for high-accuracy human MEG/EEG",<i> Journal of Neuroscience Methods</i>)</br>
  </a>
</p>
<br></br>



# For further information, please refer to: 

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
