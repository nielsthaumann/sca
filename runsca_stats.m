function data_scastats = runsca_stats(cfg, comps, data, standard)

% Individual-level statistical thresholding of spike density component analysis (SCA) components.  
%
% Performs statistical thresholding of individual-level evoked responses measured with EEG or MEG.
% 
% Applies SCA components as spatiotemporal filters to suppress variance in interfering SCA components. 
% Uses t-test statistics optimal for the typical large samples of trials for evoked responses. 
% Implements Benjamini-Yekutieli false discovery rate for multiple SCA component testing correction. 
%
% The output data can be visualized and further processed with FieldTrip functions.
% 
% Use as
%
%   data_scastats = runsca_stats(cfg, comps, data)
%
%   or
%
%   data_scastats = runsca_stats(cfg, comps, data, standard)
%
% where cfg is a configuration structure,
% comps are SCA components obtained with RUNSCA.
% data are epoched/segmented trials
% obtained with FT_PREPROCESSING or FT_REDEFINETRIAL.
%
% When comp is an output of RUNSCA based on average difference waves,
% e.g., for MMN studies, it is also necessary to define
% trial-level difference waves: deviant trials minus average standard trials.
% In this case, data are epoched/segmented deviant trials obtained with FT_PREPROCESSING or FT_REDEFINETRIAL,
% and standard are epoched/segmented standard trials obtained with FT_PREPROCESSING or FT_REDEFINETRIAL.
%
% 
% The configuration or part of the configuration can simply be empty (e.g., cfg = []),
% in which case default settings are applied (see below).
% 
% SCA statistics setting
% 
% cfg.p_thres        = p-values below the defined threshold indicate statistical significance (default is .01)
% 
% Region of interest settings
%
% SCA components are excluded if they do not...
% cfg.channels       = peak in any of the defined channels (e.g., {'F3','Fz','F4'})
%                          (see FT_CHANNELSELECTION) (default is all channels)
% cfg.polarity       = show a negative peak = -1, positive peak = +1, or peak of any polarity = 0 (default is any polarity)
% cfg.latency        = peak in the latency range [min max] in milliseconds (e.g., [100 200]) (default is total peristimulus time range after 0 ms)
% 
% Visualization settings
%
% cfg.showsteps      = visually inspect the steps in the SCA statistics procedure, true or false (default)
% cfg.ncomps         = maximum number of SCA components to show when visualizing the steps in the SCA statistics procedure (default is 15) (use inf to show all tested SCA components)
% cfg.showresult     = visually inspect the SCA statistics result, true or false (default)
% cfg.layout         = is needed for inspecting the topographies. Use the output of FT_PREPARE_LAYOUT as input
%
% 
% Output data from the SCA statistics:  
% 
% data_scastats.testedcomps            = labels for identification of the tested SCA components (by peak latency, peak channel label, and peak amplitude)
% data_scastats.p_value_uncor          = uncorrected p-values for the tested SCA components
% data_scastats.mean                   = mean amplitude area measures for the tested SCA components
% data_scastats.stdev                  = standard deviations of the amplitude area measures for the tested SCA components
% data_scastats.n_trials               = number of trials (sample size) for the tested SCA components
% data_scastats.t_value                = t-statistics for the tested SCA components
% data_scastats.ci_ratio               = upper and lower confidence interval ratios in relation to the mean area measures for the tested SCA components
% data_scastats.significant_uncor      = uncorrected significance of the tested SCA components
% data_scastats.significant_byfdr      = Benjamini-Yekutieli false discovery rate corrected significance of the tested SCA components
% 
% data_scastats.time                   = time course for the summed statistically thresholded SCA components (enables visualization and further processing with FieldTrip functions)
% data_scastats.avg                    = average waveform (channel x time matrix) for the summed statistically thresholded SCA components (enables visualization and further processing with FieldTrip functions)
% data_scastats.label                  = channel labels for the summed statistically thresholded SCA components (enables visualization and further processing with FieldTrip functions)
% data_scastats.dimord                 = order of dimentions for the summed statistically thresholded SCA components (enables visualization and further processing with FieldTrip functions)
% data_scastats.ci_upper               = upper confidence intervals for the summed statistically thresholded SCA components (enables visualization and further processing with custom code)
% data_scastats.ci_lower               = lower confidence intervals for the summed statistically thresholded SCA components (enables visualization and further processing with custom code)
% 
% 
% This function runs in the Matlab environment and requires the FieldTrip toolbox to be installed (see https://github.com/fieldtrip).
%
%
% Beta version 20250113.
%
% When applying this function the following publication must be cited:
%
% Haumann, N T; Petersen, B; Friis Andersen, A S; Faulkner, K S; Brattico, E; Vuust, P;
% "Mismatch negativity as a marker of music perception in individual cochlear implant users: 
% A spike density component analysis study",
% Clinical Neurophysiology (2023), https://doi.org/10.1016/j.clinph.2023.01.015
%

%% Prepare the input data and settings

fprintf('\nRunning SCA statistics...\n\n')

% Verify required inputs are provided
if nargin<3
    error('Please provide the required inputs (cfg, comps, data, ...). Type help runsca_stats for more information.')
end

% Verify second input is SCA components
if ~isfield(comps,'sigma')
    error('Please provide SCA components as the second input. Type help runsca_stats for more information.')
end

% Verify third input is epoched data
if isfield(data,'trial')
    if ~(length(data.trial)>1)
        error('Please epoched/segmented trials as the third input. Type help runsca_stats for more information.')
    end
elseif ~isfield(data,'trial')
    error('Please epoched/segmented trials as the third input. Type help runsca_stats for more information.')
end

% Verify fourth input is epoched data
if nargin==4
    if isfield(standard,'trial')
        if ~(length(standard.trial)>1)
            error('Please epoched/segmented standard trials as the fourth input. Type help runsca_stats for more information.')
        end
    elseif ~isfield(standard,'trial')
        error('Please epoched/segmented standard trials as the fourth input. Type help runsca_stats for more information.')
    end
end

% Verify FieldTrip is installed
if ~exist('ft_getopt','file')
    error('Could not find the function ''ft_getopt''. Please ensure that the FieldTrip Toolbox is installed, and related functions are added to the paths with ''ft_defaults''.')
end

% Ensure the sampling rate is the same for comps and data (and standard)
if 1/diff(comps.time{1,1}(1:2)) ~= 1/diff(data.time{1,1}(1:2))
    error(['The sampling rate of the comps (',num2str(1/diff(comps.time{1,1}(1:2))),' Hz) and data (',num2str(1/diff(data.time{1,1}(1:2))),' Hz) must be the same.'])
end
if nargin==4
    if 1/diff(comps.time{1,1}(1:2)) ~= 1/diff(data.time{1,1}(1:2)) || 1/diff(comps.time{1,1}(1:2)) ~= 1/diff(standard.time{1,1}(1:2)) || 1/diff(data.time{1,1}(1:2)) ~= 1/diff(standard.time{1,1}(1:2))
        error(['The sampling rate of the comps (',num2str(1/diff(comps.time{1,1}(1:2))),' Hz), data (',num2str(1/diff(data.time{1,1}(1:2))),' Hz), and standard (',num2str(1/diff(standard.time{1,1}(1:2))),' Hz) must be the same.'])
    end
end

% Ensure to process only the shared comps and data (and standard) channels
channel_labels = intersect(comps.topolabel, data.label); % Shared channel labels
if nargin==3
    if isempty(channel_labels)
        error('The comps and data must have shared channels.')
    end
    if length(comps.topolabel) ~= length(data.label)
        warning(['The comps and data have different channels. Using only the ',num2str(length(channel_labels)),' shared channels.'])
    end
    
    % Correct the order of the component channels according to the shared channels list
    comps_channels = [];
    for i=1:length(channel_labels)
        comps_channels(i) = find(ismember(comps.topolabel,channel_labels{i}));
    end
    comps.topolabel = comps.topolabel(comps_channels); 
    comps.topo = comps.topo(comps_channels, :);
    comps.unmixing = comps.unmixing(:, comps_channels);
    clear('comps_channels'); % Cleanup memory
    
    % Correct the order of the data channels according to the shared channels list
    data_channels = [];
    for i=1:length(channel_labels)
        data_channels(i) = find(ismember(data.label,channel_labels{i}));
    end
    data.label = data.label(data_channels);
    for i=1:length(data.trial)
        data.trial{i} = data.trial{i}(data_channels, :);
    end
    clear('data_channels'); % Cleanup memory
    
elseif nargin==4
    
    channel_labels = intersect(channel_labels, standard.label); % Shared channel labels, also constrained by standard channels
    if isempty(channel_labels)
        error('The comps, data, and standard must have shared channels.')
    end
    if ( length(comps.topolabel) ~= length(data.label) ) || ( length(comps.topolabel) ~= length(standard.label) ) || ( length(data.label) ~= length(standard.label)  )
        warning(['The comps, data, and standard have different channels. Using only the ',num2str(length(channel_labels)),' shared channels.'])
    end
    
    % Correct the order of the component channels according to the shared channels list
    comps_channels = [];
    for i=1:length(channel_labels)
        comps_channels(i) = find(ismember(comps.topolabel,channel_labels{i}));
    end
    comps.topolabel = comps.topolabel(comps_channels); 
    comps.topo = comps.topo(comps_channels, :);
    comps.unmixing = comps.unmixing(:, comps_channels);
    clear('comps_channels'); % Cleanup memory
    
    % Correct the order of the data channels according to the shared channels list
    data_channels = [];
    for i=1:length(channel_labels)
        data_channels(i) = find(ismember(data.label,channel_labels{i}));
    end
    data.label = data.label(data_channels);
    for i=1:length(data.trial)
        data.trial{i} = data.trial{i}(data_channels, :);
    end
    clear('data_channels'); % Cleanup memory
    
    % Correct the order of the standard channels according to the shared channels list
    standard_channels = [];
    for i=1:length(channel_labels)
        standard_channels(i) = find(ismember(standard.label,channel_labels{i}));
    end
    standard.label = standard.label(standard_channels);
    for i=1:length(standard.trial)
        standard.trial{i} = standard.trial{i}(standard_channels, :);
    end
    clear('standard_channels'); % Cleanup memory
    
end
clear('channel_labels'); % Cleanup memory

% Ensure to process only the shared comps and data (and standard) time samples
time = intersect(round(comps.time{1,1}*1000), round(data.time{1,1}*1000))/1000; % Shared time points
if nargin==3
    
    if isempty(time)
        error('The comps and data must have shared time points.')
    end
    if length(comps.time{1,1}) ~= length(data.time{1,1})
        warning(['The comps and data have different time points. Using only the ',num2str(length(time)),' shared time points.'])
    end
    
    % Trim the component time samples to the shared time points
    comps_samples = find(ismember(round(comps.time{1,1}*1000), round(time*1000)));
    comps.time{1} = comps.time{1}(comps_samples);
    comps.trial{1} = comps.trial{1}(:, comps_samples);
    clear('comps_samples'); % Cleanup memory
    
    % Trim the data time samples to the shared time points
    data_samples = find(ismember(round(data.time{1,1}*1000), round(time*1000)));
    for i=1:length(data.trial)
        data.trial{i} = data.trial{i}(:, data_samples); 
    end
    for i=1:length(data.time)
        data.time{i} = data.time{i}(:, data_samples); 
    end
    if isfield(data, 'sampleinfo')
        for i=1:size(data.sampleinfo,1)
            data.sampleinfo(i,:) = data.sampleinfo(i,1) + [data_samples(1)-1, data_samples(end)];
        end
    end
    clear('data_samples'); % Cleanup memory
    
elseif nargin==4
    
    time = intersect(round(time*1000), round(standard.time{1,1}*1000))/1000; % Shared time points, also contrained by standard time points
    if isempty(time)
        error('The comps, data, and standard must have shared time points.')
    end
    if ( length(comps.time{1,1}) ~= length(data.time{1,1}) ) || ( length(comps.time{1,1}) ~= length(standard.time{1,1}) ) || ( length(data.time{1,1}) ~= length(standard.time{1,1})  )
        warning(['The template, comps, and standard have different time points. Using only the ',num2str(length(time)),' shared time points.'])
    end
    
    % Trim the component time samples to the shared time points
    comps_samples = find(ismember(round(comps.time{1,1}*1000), round(time*1000)));
    comps.time{1} = comps.time{1}(comps_samples);
    comps.trial{1} = comps.trial{1}(:, comps_samples);
    clear('comps_samples'); % Cleanup memory
    
    % Trim the data time samples to the shared time points
    data_samples = find(ismember(round(data.time{1,1}*1000), round(time*1000)));
    for i=1:length(data.trial)
        data.trial{i} = data.trial{i}(:, data_samples); 
    end
    for i=1:length(data.time)
        data.time{i} = data.time{i}(:, data_samples); 
    end
    if isfield(data, 'sampleinfo')
        for i=1:size(data.sampleinfo,1)
            data.sampleinfo(i,:) = data.sampleinfo(i,1) + [data_samples(1)-1, data_samples(end)];
        end
    end
    clear('data_samples'); % Cleanup memory
    
    % Trim the standard time samples to the shared time points
    standard_samples = find(ismember(round(standard.time{1,1}*1000), round(time*1000)));
    for i=1:length(standard.trial)
        standard.trial{i} = standard.trial{i}(:, standard_samples); 
    end
    for i=1:length(standard.time)
        standard.time{i} = standard.time{i}(:, standard_samples); 
    end
    if isfield(standard, 'sampleinfo')
        for i=1:size(standard.sampleinfo,1)
            standard.sampleinfo(i,:) = standard.sampleinfo(i,1) + [standard_samples(1)-1, standard_samples(end)];
        end
    end
    clear('standard_samples'); % Cleanup memory
    
end

% Check whether any time steps are non-consecutive (>=1.5 sampling interval)
percent_consecutive_time = 100 - ( sum( diff(time) >=  1.5*diff(comps.time{1,1}(1:2)) )/length(time)*100 ); 
if percent_consecutive_time ~= 100
    warning([num2str(100 - percent_consecutive_time),'% of the shared time points are non-consecutive. There is up to ',num2str( max(abs(diff(time) - diff(comps.time{1,1}(1:2)))) * 1000),' ms jitter in the time steps within the applied time range. (This issue might be related to diverging preprossesing procedures for the inputs.)'])
end
clear('time','percent_consecutive_time'); % Cleanup memory

% Calculate trial-level difference waves
if nargin==4
    
    avg_standard = nanmean(cat(3,standard.trial{:}),3);
    disp('Subtracting the average standard from each deviant trial.')
    for i=1:length(data.trial)
        data.trial{i} = data.trial{i} - avg_standard;
    end
    clear('standard', 'avg_standard'); % Cleanup memory
    
end

% Convert epochs into data matrix (channel, time, trial)
data_matrix = cat(3,data.trial{:});
clear('data'); % Cleanup memory


%% Region of interest analysis

% Apply the default settings if no ROI constraints are defined
if ~isfield(cfg,'channels') 
    cfg.channels = comps.topolabel; % Default setting is all channels
end
if ~isfield(cfg,'polarity') 
    cfg.polarity = 0; % default seting is any polarity
end
if ~isfield(cfg,'latency') 
    cfg.latency = [0 max(comps.time{1})]*1000; % default is total peristimulus time range after 0 ms
end

% Find components with matching peak channels and polarity
if cfg.polarity == 0
    peak_channel_id = [];
    [~,peak_channel_id(:,1)] = min(comps.topo,[],1); % Peak amplitude that is negative...
    [~,peak_channel_id(:,2)] = max(comps.topo,[],1); % or positive
    peak_channel_id = unique(peak_channel_id);
else
    [~,peak_channel_id] = max(cfg.polarity*comps.topo,[],1); % Peak amplitude of specified polarity
end
peak_channels = comps.topolabel(peak_channel_id);
roi_match_peak_channel = find(ismember(peak_channels,cfg.channels));

% Find components with matching peak latencies
roi_match_peak_latency = find(comps.latency >= cfg.latency(1) & comps.latency <= cfg.latency(2));

roi_match = intersect(roi_match_peak_channel,roi_match_peak_latency);

clear('peak_channel_id', 'peak_channels', 'roi_match_peak_channel', 'roi_match_peak_latency'); % Cleanup memory

fprintf('\nROI channels: Peak in')
for i=1:length(cfg.channels)
    fprintf([' ', cfg.channels{i}])
end
fprintf('.\n')
if cfg.polarity==-1
    disp('ROI polarity: Negative peak.')
elseif cfg.polarity==+1
    disp('ROI polarity: Positive peak.')
else
    disp('ROI polarity: Any negative or positive peak.')
end
disp(['ROI latency : ',num2str(cfg.latency(1)),' - ',num2str(cfg.latency(2)),' ms peak latency.'])

disp(['Excluded ',num2str(size(comps.trial{1,1},1)-length(roi_match)),' SCA components outside the region of interest.'])
disp(['Keeping ',num2str(length(roi_match)), ' SCA components for statistical analysis.'])
fprintf('\n')


%% Visualize the SCA components to be tested statistically

% Apply the default settings if no visualization settings are defined
if ~isfield(cfg,'showsteps') 
    cfg.showsteps = false; % Default setting is no inspection
end
if ~isfield(cfg,'showresult') 
    cfg.showresult = false; % Default setting is no inspection
end
if ~isfield(cfg,'ncomps')
    cfg.ncomps = 15; % Default setting is 15
end
if cfg.ncomps==inf
    cfg.ncomps = length(comps.label); % Maximum is number of SCA components
end
if ~isfield(cfg,'layout')
    cfg.layout = []; % Default is no layout provided
end

if cfg.showsteps && ~isempty(roi_match)
    
    % Prepare the required resolution of subplots for SCA component visualizations
    inspect_res = ceil(sqrt(cfg.ncomps)); 
    
    % Topographies before ROI constraints
    figure('name','SCA topographies before ROI constraints','color','w');
    if ~isempty(cfg.layout)
        disp('Showing SCA topographies before ROI constraints.')
        cfg_inspect = struct;
        cfg_inspect.channel = comps.topolabel;
        cfg_inspect.colorbar = 'no';
        cfg_inspect.comment = ' ';
        cfg_inspect.baseline = 'no';
        cfg_inspect.layout = cfg.layout;
        cfg_inspect.zlim = 'maxabs';
        for i=1:min(cfg.ncomps, length(comps.label))
            cfg_inspect.component = i; 
            subplot(inspect_res,inspect_res,i)
            evalc('ft_topoplotER(cfg_inspect, comps)');
            title(comps.label(i),'fontsize',13)
        end
    else
        warning('No channel layout was provided! Skipping topography plotting.')
    end
    
    % Waveforms before ROI constraints
    disp('Showing SCA waveforms before ROI constraints.')
    figure('name','SCA waveforms before ROI constraints','color','w');
    plot(comps.time{1}*1000, comps.trial{1}(1:min(cfg.ncomps, length(comps.label)),:)) % Butterfly plot
    hold on
    plot(comps.time{1}*1000, zeros(length(comps.time),1),'black') % Draw baseline
    legend(comps.label(1:min(cfg.ncomps, length(comps.label))),'location','NorthEastOutside','fontsize',13)
    ylabel('Amplitude (at peak channel)','fontsize',13)
    xlabel('Time (ms)','fontsize',13)
    set(gca,'fontsize',13)
    title('SCA waveforms before ROI constraints','fontsize',13)
    
    % Topographies after ROI constraints
    figure('name','SCA topographies after ROI constraints','color','w');
    if ~isempty(cfg.layout)
        disp('Showing SCA topographies after ROI constraints.')
        cfg_inspect = struct;
        cfg_inspect.channel = comps.topolabel;
        cfg_inspect.colorbar = 'no';
        cfg_inspect.comment = ' ';
        cfg_inspect.baseline = 'no';
        cfg_inspect.layout = cfg.layout;
        cfg_inspect.zlim = 'maxabs';
        for i=1:min(length(roi_match), cfg.ncomps)
            cfg_inspect.component = roi_match(i); 
            subplot(inspect_res,inspect_res,i)
            evalc('ft_topoplotER(cfg_inspect, comps)');
            title(comps.label(roi_match(i)),'fontsize',13)
        end
    else
        warning('No channel layout was provided! Skipping topography plotting.')
    end
    
    % Waveforms after ROI
    fprintf('Showing SCA waveforms after ROI constraints.\n\n')
    figure('name','SCA waveforms after ROI constraints','color','w');
    plot(comps.time{1}*1000, comps.trial{1}(roi_match(1:min(cfg.ncomps, length(roi_match))),:)) % Butterfly plot
    hold on
    plot(comps.time{1}*1000, zeros(length(comps.time),1),'black') % Draw baseline
    legend(comps.label(roi_match(1:min(length(roi_match), cfg.ncomps))),'location','NorthEastOutside','fontsize',13)
    ylabel('Amplitude (at peak channel)','fontsize',13)
    xlabel('Time (ms)','fontsize',13)
    set(gca,'fontsize',13)
    title('SCA waveforms after ROI constraints','fontsize',13)    
    
    clear('cfg_inspect', 'inspect_res'); % Cleanup memory
    
end


%% Perform the SCA statistics

% Apply the default setting if none is defined
if ~isfield(cfg,'p_thres') 
    cfg.p_thres = .01; % Default setting is .01
end

% Prepare the output
data_scastats = struct; 
data_scastats.testedcomps        = comps.label(roi_match)'; % Store the labels for the tested SCA components
data_scastats.p_value_uncor      = [];
data_scastats.mean               = [];
data_scastats.stdev              = [];
data_scastats.n_trials           = [];
data_scastats.t_value            = [];
data_scastats.ci_ratio           = [];
data_scastats.significant_uncor  = [];
data_scastats.significant_byfdr  = [];
data_scastats.time               = [];
data_scastats.avg                = [];
data_scastats.label              = [];
data_scastats.dimord             = [];
data_scastats.ci_upper           = [];
data_scastats.ci_lower           = [];

if cfg.showsteps && ~isempty(roi_match)
    
    % Prepare the required resolution of subplots for SCA component visualizations
    inspect_res = ceil(sqrt(min(cfg.ncomps,length(roi_match))+1)); % (+1 for a legend)
    
    % Prepare visualizing the spatio-temporal filtering and statistical testing
    spattempfilt = figure('name','Spatio-temporal filtering','color','white');
    distribution = figure('name','Trial area measure distributions','color','white');
    
end

% Prepare variable for storing the area measure for each tested SCA component and trial (component, trial)
area_measure = zeros( length(roi_match) , size(data_matrix, 3) );

% Main loop across the tested SCA components
for i=1:length(roi_match)
    
    fprintf(['Testing SCA component ',num2str(i),' of ',num2str(length(roi_match)),'.\n'])
    
    disp('Suppressing inter-trial variance originating from interfering SCA components and channel noise.')
    
    % Project the tested SCA component into a channel x time x trial matrix
    component_filter = permute( repmat( comps.trial{1}(roi_match(i),:)' , [ 1 , length(comps.topolabel) , size(data_matrix,3) ] ) , [2, 1, 3] ) .* repmat( comps.topo(:,roi_match(i)) , [ 1 , length(comps.time{1}) , size(data_matrix,3) ]); % (channel, time, trial)
    
    % Apply the spatiotemporal filter on the data matrix
    filtered_comp = sign(data_matrix.*component_filter) .* sqrt(abs(data_matrix.*component_filter)); % (channel, time, trial)
    
    if cfg.showsteps && i<=cfg.ncomps
        
        % Show the effects of suppresing interfering SCA components and their variance
        figure(spattempfilt)
        subplot(inspect_res,inspect_res,i)
        
        % Show the unprocessed trials
        peak_channel_index = ismember(comps.topolabel, comps.channel{roi_match(i)});
        plot_data = squeeze(median(data_matrix(peak_channel_index,:,:),3)); 
        error_data = []; 
        error_data(1,:) = squeeze( prctile(data_matrix(peak_channel_index,:,:),75,3) ) - squeeze( median(data_matrix(peak_channel_index,:,:),3) ); % Upper quartile difference to median
        error_data(2,:) = squeeze( median(data_matrix(peak_channel_index,:,:),3)) - squeeze( prctile(data_matrix(peak_channel_index,:,:),25,3) ); % Lower quartile difference to median
        plot(comps.time{1}*1000, plot_data, 'b')
        patch([comps.time{1}, comps.time{1}(end:-1:1)]*1000, [plot_data + error_data(1,:), plot_data(end:-1:1) - error_data(2,end:-1:1)], 'b', 'FaceAlpha',.2, 'EdgeColor', 'none')
        box off
        
        % Show the trials after the SCA component spatiotemporal filter
        hold on
        plot_data = comps.sign(roi_match(i)) * sum( filtered_comp .* abs( repmat(comps.topo(:,roi_match(i)), [1, length(comps.time{1}), size(data_matrix,3)])), 1) / sum( max(comps.trial{1}(roi_match(i),:)) * comps.topo(:,roi_match(i)).^2 ) * max(comps.trial{1}(roi_match(i),:)); % Convert into component space corresponding to peak channel (amplitude corresponds to peak channel)
        error_data = []; 
        error_data(1,:) = squeeze( prctile(plot_data,75,3)) - squeeze(median(plot_data,3) ); % Upper quartile difference to median
        error_data(2,:) = squeeze( median(plot_data,3)) - squeeze(prctile(plot_data,25,3) ); % Lower quartile difference to median
        plot(comps.time{1}*1000, squeeze(median(plot_data,3)), 'color', 'black', 'LineWidth',2)
        patch([comps.time{1}, comps.time{1}(end:-1:1)]*1000, [squeeze(median(plot_data,3)) + error_data(1,:), squeeze(median(plot_data(1,end:-1:1,:),3)) - error_data(2,end:-1:1)], 'black', 'FaceAlpha',.2, 'EdgeColor', 'none')
        box off
        
        title([comps.label{roi_match(i)},' median (+/-quartile) of trials'],'fontsize',13)
        xlabel('Time (ms)','fontsize',13)
        ylabel('Amplitude','fontsize',13)
        set(gca,'fontsize',13)
        xlim([min(comps.time{1}) max(comps.time{1})]*1000)
        if comps.sign(roi_match(i))==-1
            set(gca,'Ydir','reverse')
        end
        
        clear('peak_channel_index', 'plot_data', 'error_data'); % Cleanup memory
        
    end
    
    % Filtered SCA component area measure over expected filtered SCA component area measure multiplied by the expected amplitude
    disp('Estimating the summed area measure for the SCA component for each trial.')
    area_measure(i,:) = max(comps.trial{1}(roi_match(i),:)) * squeeze( sum( sum( filtered_comp ,1) ,2) ) / sum(sum( abs(component_filter(:,:,1)) )); 
    
    clear('component_filter', 'filtered_comp'); % Cleanup memory
    
    % Estimate the statistical significance of the SCA component with one-sample t-test
    disp('Estimating the statistical significance of the SCA component with one-sample t-test.')
    
    % One-sample t-test in the rightwards (positive) direction 
    % (Note: The two-tailed or negative directions are meaningless, 
    %  since they would allow an imagery SCA component of opposite polarity to be tested.)
    [~,data_scastats.p_value_uncor(i),~,stats] = ttest( area_measure(i,:)' , 0, 'Alpha',cfg.p_thres, 'Tail','right');
    
    % Store the test information on mean, stdev, n (total trials), t-statistic, and confidence intervals
    data_scastats.mean = mean(area_measure, 2)'; % Mean
    data_scastats.stdev(i) = stats.sd; % Standard deviation
    data_scastats.n_trials(i) = length(area_measure(i,:)); % Number of trials = sample size
    data_scastats.t_value(i) = stats.tstat; % t-statistic
    t_thres = tinv((1 - cfg.p_thres / 1 ), data_scastats.n_trials(i)-1 ); % t-treshold for the p-threshold
    sem = data_scastats.stdev(i) / ( sqrt(data_scastats.n_trials(i)) ); % Standard error of the mean
    data_scastats.ci_ratio(:,i) = mean(area_measure(i,:)) + t_thres * sem * [ +1 ; -1 ]; % Confidence interval ratio to mean amplitude
    
    if cfg.showsteps && i<=cfg.ncomps
        
        % Prepare showing asterisks symbol (*) if SCA component was significant
        if data_scastats.p_value_uncor(i)<cfg.p_thres
            significant = '*';
        else
            significant = '';
        end
        
        % Show the distribution of the summed area amplitude in blue color
        figure(distribution)
        subplot(inspect_res,inspect_res,i)
        [count, amplitude] = hist( comps.sign(roi_match(i) ) * ( area_measure(i,:) ));
        bar(amplitude, count, 'b')
        
        % Show the observations contradicting the rejection of the null-hypothesis in red color
        if comps.sign(roi_match(i)) == -1 % if significant in negative direction
            contraindex = amplitude >= 0;
        elseif comps.sign(roi_match(i)) == +1 % if significant in positive direction
            contraindex = amplitude <= 0;
        end
        if sum(contraindex) > 0
            hold on
            bar(amplitude, count .* contraindex ,'r')
        end
        
        % Ensure the amplitude of the tested SCA component always is shown increasing in the rightwards direction on the histogram
        if comps.sign(roi_match(i)) == -1
            set(gca,'Xdir','reverse')
        end
        
        title([comps.label{roi_match(i)},'; p=', num2str(data_scastats.p_value_uncor(i)), significant],'fontsize',13)
        xlabel('Amplitude','fontsize',13)
        ylabel('Number of trials','fontsize',13)
        set(gca,'fontsize',13)
        xlim( max(abs(get(gca,'Xlim'))) * [-1 1] )
        
        clear('significant', 'count', 'amplitude', 'contraindex'); % Cleanup memory
        
    end
    
    if data_scastats.p_value_uncor(i)<cfg.p_thres
        fprintf(['* Significant at p<',num2str(cfg.p_thres),' (without multiple testing correction).\n\n'])
    else
        fprintf(['Not significant at p<',num2str(cfg.p_thres),'.\n\n'])
    end
    
    clear('sem', 'stats', 't_thres'); % Cleanup memory
    
end

clear('area_measure'); % Cleanup memory

if cfg.showsteps && ~isempty(roi_match)
    
    % Add legend to the spatiotemporal SCA filter figure
    figure(spattempfilt)
    subplot(inspect_res, inspect_res, min(cfg.ncomps,length(roi_match))+1)
    plot(comps.time{1}*1000, nan*ones(1,length(comps.time{1})),'b')
    hold on
    plot(comps.time{1}*1000, nan*ones(1,length(comps.time{1})), 'black','linewidth',2)
    axis off
    legend({'Unprocessed trials','Spatiotemporal SCA filter'},'fontsize',13)
    
    clear('inspect_res'); % Cleanup memory
    
end

% Apply the Benjamini-Yekutieli false discovery rate (BY-FDR) correction to estimate the significance of the tested SCA components
fprintf('Applying the Benjamini-Yekutieli false discovery rate correction. \n')
data_scastats.significant_uncor( ~isnan(data_scastats.p_value_uncor) ) = data_scastats.p_value_uncor( ~isnan(data_scastats.p_value_uncor) ) < cfg.p_thres; % Significant uncorrected p-values
data_scastats.significant_byfdr( ~isnan(data_scastats.p_value_uncor) ) = fdr_by(data_scastats.p_value_uncor( ~isnan(data_scastats.p_value_uncor) ), cfg.p_thres); % Significant BY-FDR corrected p-values
data_scastats.significant_byfdr( isnan(data_scastats.p_value_uncor) ) = false;


%% Handle the resulting output from the SCA statistics

% Extract the significant SCA components
fprintf(['\nExtracting ',num2str(sum(data_scastats.significant_byfdr)),' statistically significant SCA components, \n'])
data_scastats.time = comps.time{1,1};
data_scastats.avg = zeros(length(comps.topolabel),length(comps.time{1}));
data_scastats.label = comps.topolabel;
data_scastats.dimord = 'chan_time';
significant_comp_index = find(data_scastats.significant_byfdr);
data_scastats.avg = zeros(length(comps.topolabel), size(comps.trial{1},2)); % Prepare the statistically thresholded output data (channel, time)
for i=1:length(significant_comp_index) % Loop across significant SCA components within the ROI
    
    component_projection = repmat(comps.topo(:,roi_match(significant_comp_index(i))), [1, size(comps.trial{1},2)]) .* repmat(comps.trial{1}(roi_match(significant_comp_index(i)),:), [length(comps.topolabel), 1]); % Project the significant SCA component into channel space
    data_scastats.avg = data_scastats.avg + component_projection;
    
end
clear('significant_comp_index', 'component_projection'); % Cleanup memory

% Project the confidence intervals (-/+% of amplitude) onto the channels using the CI for the significant SCA component with maximum SNIR
disp('and adding estimated confidence intervals. ')
if sum(data_scastats.significant_byfdr)>0
    
    significant_comps = find(data_scastats.significant_byfdr); 
    [~,max_snir_index] = max( data_scastats.t_value(significant_comps) ); % Find the SCA component with maximum signal-to-noise-and-interference-ratio based on the t-value
    
    data_scastats.ci_upper = data_scastats.ci_ratio(1, significant_comps(max_snir_index)) / data_scastats.mean(significant_comps(max_snir_index)) * data_scastats.avg;
    data_scastats.ci_lower = data_scastats.ci_ratio(2, significant_comps(max_snir_index)) / data_scastats.mean(significant_comps(max_snir_index)) * data_scastats.avg;
    
    clear('max_snir_index', 'significant_comps')
    
else
    
    data_scastats.ci_upper = data_scastats.avg;
    data_scastats.ci_lower = data_scastats.avg;
    
end

if cfg.showsteps && sum(data_scastats.significant_byfdr)>0
    
    % Show the Significant SCA components with confidence intervals
    
    figure('name','Significant SCA components with confidence intervals','color','w')
    [~,peak_channel_index] = max(max(abs(data_scastats.avg),[],2),[],1);
    [~,time_sample_index] = max(abs(data_scastats.avg(peak_channel_index,:)));
    plot(data_scastats.time*1000, data_scastats.avg(peak_channel_index,:), 'k')
    patch([comps.time{1}, comps.time{1}(end:-1:1)]*1000, [data_scastats.ci_upper(peak_channel_index,:), data_scastats.ci_lower(peak_channel_index,end:-1:1)], 'black', 'FaceAlpha',.2, 'EdgeColor', 'none')
    box off
    xlim([data_scastats.time(1) data_scastats.time(end)]*1000)
    xlabel('Time (ms)','fontsize',13)
    ylabel(['Amplitude (peak channel: ',data_scastats.label{peak_channel_index},')'],'fontsize',13)
    set(gca,'fontsize',13)
    title('Significant SCA components with confidence intervals','interpreter','none','fontsize',13)
    if data_scastats.avg(peak_channel_index, time_sample_index)<0
        set(gca,'Ydir','reverse')
    end
    legend({'Average', [num2str(100*(1-cfg.p_thres)),'% CI']}, 'location','NorthWestOutside','fontsize',13)
    
    clear('peak_channel_index', 'time_sample_index')
    
end

% Show the results of the SCA statistics

if cfg.showresult && sum(data_scastats.significant_byfdr)>0
    
    if ~isempty(cfg.layout)
        
        % Show the individual statistically thresholded SCA components
        fprintf('\nStarting ft_databrowser to show the individual statistically thresholded SCA components.\n')
        cfg_inspect = struct;
        cfg_inspect.layout = cfg.layout;
        cfg_inspect.channel = roi_match(find(data_scastats.significant_byfdr));
        cfg_inspect.viewmode = 'component';
        cfg_inspect.zlim = 'maxabs';
        evalc('ft_databrowser(cfg_inspect, comps)');
        title('Separated statistically thresholded SCA components','fontsize',13)
        
        % Show the sum of the statistically thresholded SCA components
        fprintf('\nStarting ft_multiplotER to show the sum of statistically thresholded SCA components.\n')
        figure('name','Sum of statistically thresholded SCA components','color','w')
        cfg_inspect = struct;
        cfg_inspect.layout = cfg.layout;
        evalc('ft_multiplotER(cfg_inspect, data_scastats)');
        title('Sum of statistically thresholded SCA components','fontsize',13)
        
        clear('cfg_inspect'); % Cleanup memory
        
    else
        warning('No channel layout was provided! Skipping plotting of separated statistically thresholded SCA components.')
    end
    
end

if sum(data_scastats.significant_byfdr)>0
    fprintf(['\n* Found significant SCA components at p<',num2str(cfg.p_thres),' after multiple testing correction.\n\n'])
else
    fprintf(['\nNo significant SCA components were found at p<',num2str(cfg.p_thres),' after multiple testing correction.\n\n'])
end


%% Benjamini-Yekutieli false discovery rate function

function [h] = fdr_by(p, q)

% FDR false discovery rate
%
% Use as
%   h = fdr_by(p, q)
%
% The input argument p is a vector or matrix with (uncorrected) p-values, the input argument
% q is a scalar that reflects the critical alpha-threshold for the inferential decision. The
% output argument h is a boolean matrix (same size as p) denoting for each sample whether 
% the null hypothesis can be rejected. 
%
% This implements
%   Genovese CR, Lazar NA, Nichols T.
%   Thresholding of statistical maps in functional neuroimaging using the false discovery rate.
%   Neuroimage. 2002 Apr;15(4):870-8.
%
% There are two types of FDR correction (Benjamini-Hochberg & Benjamini-Yekutieli), of
% which the second is currently implemented.
% 
% Correction by NTH 20210421: p-values must be less than the corrected
% threshold (not equal to or less than). 

% Copyright (C) 2005-2015, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$
% 

% convert the input into a row vector
dim = size(p);
p = reshape(p, 1, numel(p));

% sort the observed uncorrected probabilities
[ps, indx] = sort(p);

% count the number of voxels
V = length(p);

% compute the threshold probability for each voxel
pi = ((1:V)/V)  * q / c(V);

% if any(ps<=pi)
if any(ps<pi)
  h = ps<=(max(ps(ps<pi)));
else
  h = false(size(ps));
end

% undo the sorting
[dum, unsort] = sort(indx);
h = h(unsort);

% convert the output back into the original format
h = reshape(h, dim);

function s = c(V)
% See Genovese, Lazar and Holmes (2002) page 872, second column, first paragraph
if V<1000
  % compute it exactly
  s = sum(1./(1:V));
else
  % approximate it
  s = log(V) + 0.57721566490153286060651209008240243104215933593992359880576723488486772677766467093694706329174674951463144724980708248096050401448654283622417399764492353625350033374293733773767394279259525824709491600873520394816567;
end
