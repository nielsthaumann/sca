% Spike density component analysis decomposes EEG or MEG waveforms into neural components 
% modelled with temporal probability density functions.
% 
% It is assumed that SNIR > 1 and components differ in latency, width, or topography.
% 
% Use as
% 
%   comp = runsca(cfg, data)
% 
% where cfg is a configuration structure
% and the data is average evoked response MEG or EEG waveforms 
% obtained with FT_TIMELOCKANALYSIS or FT_TIMELOCKGRANDAVERAGE. 
% 
% The configuration or part of the configuration can simply be empty (e.g., cfg = []), 
% in which case default settings are applied (see below). 
% 
% NB: Scaling of the channel data to the correct standard unit of measurement (e.g., micro-Volt, femto-Tesla)
% is crucial for the accuracy of the component fit results! If no configuration is defined 
% a default scaling estimate will be applied, which might be inappropriate. 
% 
% 
% Settings for the input data
% 
% cfg.channel             = cell-array with channel selection (default = 'all'), see FT_CHANNELSELECTION for details
% cfg.scaling_factor      = multiply by defined number to correct measurements to appropriate unit (e.g., 10^15 to correct from T to fT)
% cfg.scaling_unit        = text-string defining the measurement unit (e.g., 'fT') 
% cfg.baseline_correct    = baseline correction based on 'mean' (default), 'median' or 'none'
% cfg.baseline_window     = time range applied for baseline correction [begin end] in seconds (default = [min 0])
% 
% Settings for the SCA analysis
% 
% cfg.search_time         = constrain the analysis to [begin end] time in seconds (default = [min max])
% cfg.model               = decompose the data into Gaussian components (default = 'gauss'), or 'gamma' or 'sine' components for comparison
% cfg.amplitude_threshold = forces SCA analysis to stop when the residual amplitude passes below defined threshold (default = 0)
% cfg.sort_by             = sort SCA components by 'amplitude' (peak amplitude), 'latency' (default) or 'variance' (variance explained)
% cfg.spatial_correlation = correlation threshold for combining spatially similar components (Pearson's r) (Default = 1. R-values <1 can be applied for simulating less accurate methods)
% 
% Visualization settings
% 
% cfg.show_components     = calls FT_DATABROWSER to show the SCA components, 'yes' or 'no' (default)
% cfg.demomode            = detailed visual demonstration of SCA analysis, 'on' or 'off' (default)
% cfg.layout              = is needed if components are shown or demo mode is on. Use the output of FT_PREPARE_LAYOUT as input
% cfg.show_residuals      = also show residuals when components are shown or demo mode is on, 'yes' (default) or 'no'
% 
% 'Runsca' requires the Curve Fitting Toolbox to be installed for Gaussian or sine models. 
% For experimental comparison with gamma models, the Statistics and Machine Learning Toolbox needs to be installed. 
% 
% Beta version 20200430.
% 
% When applying this function the following publication must be cited:
% 
% Haumann, N T; Hansen, Brian; Huotilainen, M; Vuust, P; Brattico, E;
% "Applying Stochastic Spike train theory for high-accuracy human MEG/EEG"
% Journal of Neuroscience Methods (2020), doi: https://doi.org/10.1016/j.jneumeth.2020.108743
% 
