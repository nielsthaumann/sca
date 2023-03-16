function data_hsac = runsca_hsac(cfg, comps, data, standard)

% Spike density components with half-split average consistency.
% Projects SCA components (output of the RUNSCA function) back into channel space,
% if the SCA components are consistently detected in half-split averages.
% The output data can be visualized and processed with FieldTrip functions.
%
% Use as
%
%   data_hsac = runsca_hsac(cfg, comps, data)
%
%   or
%
%   data_hsac = runsca_hsac(cfg, comps, data, standard)
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
% Region of interest settings
%
% SCA components are excluded if they do not...
% cfg.roi.channels       = peak in any of the defined channels (e.g., {'F3','Fz','F4'})
%                          (see FT_CHANNELSELECTION) (default is channels with same polarity as the peak amplitude)
% cfg.roi.polarity       = show a negative peak = -1, positive peak = +1, or peak of any polarity = 0 (default is same polarity as the peak amplitude)
% cfg.roi.latency        = peak in the latency range [min max] in seconds (e.g., [0.100 0.200]) (default is range from the peak amplitude to zero-crossing or 0 or end time boundaries)
%
% Half-split average consistency settings
%
% cfg.hsac.mc_iterations = number of Monte Carlo simulated half-split averages (default = 100)
% cfg.hsac.p_match       = significance threshold for the correlation between an SCA component and half-split average (default = .001)
% cfg.hsac.p_comp        = minimum probability that the SCA component signficantly correlates with half-split averages (HSAs)
%                          (.50 = the SCA component is significant in half of the HSAs, 1.00 = in all the HSAs) (default = .70)
% Visualization settings
%
% cfg.inspect            = visually inspect the SCA HSAC results, 'yes' or 'no' (default)
% cfg.layout             = is needed for inspecting the topographies. Use the output of FT_PREPARE_LAYOUT as input
%
% Preprocessing
% 
% cfg.demean             = 'no' or 'yes', whether to apply baseline correction to the average data (and standard) (default = 'no')
% cfg.baselinewindow     = [begin end] in seconds, the default is the complete trial (default is all time samples)
% 
% 
% This function runs in the Matlab environment and requires the FieldTrip toolbox to be installed (see https://github.com/fieldtrip).
%
%
% Beta version 20210611.
%
% When applying this function the following publication must be cited:
%
% Bruzzone, S E P; Haumann, N T; Kliuchko, M; Vuust, P, Brattico, E;
% "Applying Spike-density Component Analysis for high-accuracy auditory event-related potentials in children",
% Clinical Neurophysiology (2021), https://doi.org/10.1016/j.clinph.2021.05.007
%

warning('It is recommended applying RUNSCA_STATS instead of RUNSCA_HSAC (see Haumann et al. (2023), "Mismatch negativity as a marker of music perception in individual cochlear implant users: A spike density component analysis study", Clinical Neurophysiology, https://doi.org/10.1016/j.clinph.2023.01.015).')

%% Prepare settings

% Verify required inputs are provided
if nargin<3
    error('Please provide the required inputs (cfg, comps, data, ...). Type help runsca_hsac for more information.')
end

% Verify second input is SCA components
if ~isfield(comps,'sigma')
    error('Please provide SCA components as the second input. Type help runsca_hsac for more information.')
end

% Verify third input is epoched data
if isfield(data,'trial')
    if ~(length(data.trial)>1)
        error('Please epoched/segmented trials as the third input. Type help runsca_hsac for more information.')
    end
elseif ~isfield(data,'trial')
    error('Please epoched/segmented trials as the third input. Type help runsca_hsac for more information.')
end

% Verify fourth input is epoched data
if nargin==4
    if isfield(standard,'trial')
        if ~(length(standard.trial)>1)
            error('Please epoched/segmented standard trials as the fourth input. Type help runsca_hsac for more information.')
        end
    elseif ~isfield(standard,'trial')
        error('Please epoched/segmented standard trials as the fourth input. Type help runsca_hsac for more information.')
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
    
    % Correct the order of channel indices for matching between the comps and data channels
    comps_channels = [];
    for i=1:length(channel_labels)
        comps_channels(i) = find(ismember(comps.topolabel,channel_labels{i}));
    end
    data_channels = [];
    for i=1:length(channel_labels)
        data_channels(i) = find(ismember(data.label,channel_labels{i}));
    end
    
elseif nargin==4
    
    channel_labels = intersect(channel_labels, standard.label); % Shared channel labels, also contrained by standard channels
    if isempty(channel_labels)
        error('The comps, data, and standard must have shared channels.')
    end
    if ( length(comps.topolabel) ~= length(data.label) ) || ( length(comps.topolabel) ~= length(standard.label) ) || ( length(data.label) ~= length(standard.label)  )
        warning(['The template, comps, and standard have different channels. Using only the ',num2str(length(channel_labels)),' shared channels.'])
    end
    comps_channels = [];
    for i=1:length(channel_labels)
        comps_channels(i) = find(ismember(comps.topolabel,channel_labels{i}));
    end
    data_channels = [];
    for i=1:length(channel_labels)
        data_channels(i) = find(ismember(data.label,channel_labels{i}));
    end
    standard_channels = [];
    for i=1:length(channel_labels)
        standard_channels(i) = find(ismember(standard.label,channel_labels{i}));
    end
    
end

% Ensure to process only the shared comps and data (and standard) time samples
time = intersect(round(comps.time{1,1}*1000), round(data.time{1,1}*1000))/1000; % Shared time points
if nargin==3
    
    if isempty(time)
        error('The comps and data must have shared time points.')
    end
    if length(comps.time{1,1}) ~= length(data.time{1,1})
        warning(['The comps and data have different time points. Using only the ',num2str(length(time)),' shared time points.'])
    end
    
    % Correct the time sample indices for matching between the comps and data time points
    comps_samples = find(ismember(round(comps.time{1,1}*1000), round(time*1000)));
    data_samples = find(ismember(round(data.time{1,1}*1000), round(time*1000)));
    
elseif nargin==4
    
    time = intersect(round(time*1000), round(standard.time{1,1}*1000))/1000; % Shared time points, also contrained by standard time points
    if isempty(time)
        error('The comps, data, and standard must have shared time points.')
    end
    if ( length(comps.time{1,1}) ~= length(data.time{1,1}) ) || ( length(comps.time{1,1}) ~= length(standard.time{1,1}) ) || ( length(data.time{1,1}) ~= length(standard.time{1,1})  )
        warning(['The template, comps, and standard have different time points. Using only the ',num2str(length(time)),' shared time points.'])
    end
    % Correct the time sample indices for matching between the comps and data time points
    comps_samples = find(ismember(round(comps.time{1,1}*1000), round(time*1000)));
    data_samples = find(ismember(round(data.time{1,1}*1000), round(time*1000)));
    standard_samples = find(ismember(round(standard.time{1,1}*1000), round(time*1000)));
    
end
percent_consecutive_time = 100 - ( sum( diff(time) >=  1.5*diff(comps.time{1,1}(1:2)) )/length(time)*100 ); % Check whether any time steps are non-consecutive (>=1.5 sampling interval)
if percent_consecutive_time ~= 100
    warning([num2str(100 - percent_consecutive_time),'% of the shared time points are non-consecutive. There is up to ',num2str( max(abs(diff(time) - diff(comps.time{1,1}(1:2)))) * 1000),' ms jitter in the time steps within the applied time range. (This issue might be related to diverging preprossesing procedures for the inputs.)'])
end


% If no configuration is provided, apply the default settings

cfg.demean              = ft_getopt(cfg, 'demean', 'no');
cfg.baselinewindow      = ft_getopt(cfg, 'baselinewindow', [min(time) max(time)]);

avg_data = nanmean(cat(3,data.trial{:}),3);
if strcmp(cfg.demean,'yes')
    disp('Applying baseline correction to the average data.')
    avg_data = avg_data - repmat(mean(avg_data(:, data.time{1,1} >= cfg.baselinewindow(1) & data.time{1,1} <= cfg.baselinewindow(2) ),2), [1, size(avg_data,2)]);
end

% Calculate trial-level difference waves
if nargin==4
    
    avg_standard = nanmean(cat(3,standard.trial{:}),3);
    if strcmp(cfg.demean,'yes')
        disp('Applying baseline correction to the average standard.')
        avg_standard = avg_standard - repmat(mean(avg_standard(:, standard.time{1,1} >= cfg.baselinewindow(1) & standard.time{1,1} <= cfg.baselinewindow(2) ),2), [1, size(avg_data,2)]);
    end
    disp('Subtracting the average standard from each deviant trial.')
    for i=1:length(data.trial)
        data.trial{i}(data_channels,data_samples) = data.trial{i}(data_channels,data_samples) - avg_standard(standard_channels,standard_samples);
    end
    avg_data(data_channels,data_samples) = avg_data(data_channels,data_samples) - avg_standard(standard_channels,standard_samples);
    
end

if ~isfield(cfg,'roi') % If no ROI settings are defined, create a ROI structure
    cfg.roi = struct;
end

if isfield(cfg.roi,'channels') % Check if exist
    valid_data_channels_labels = channel_labels(ismember(channel_labels,ft_channelselection(cfg.roi.channels, channel_labels))); % Valid shared channels within ROI to check
    valid_data_channels = [];
    for i=1:length(valid_data_channels_labels)
        valid_data_channels(i) = find(ismember(data.label,valid_data_channels_labels{i}));
    end
else
    valid_data_channels = data_channels; % If not exist, check all shared channels for now
end

if isfield(cfg.roi,'polarity') % Check if exist
    polarity = cfg.roi.polarity;
else
    polarity = 0; % If not, check any direction for now
end

if isfield(cfg.roi,'latency') % Check if exist
    valid_data_samples = find(data.time{1,1} >= cfg.roi.latency(1) & data.time{1,1} <= cfg.roi.latency(2)); % If exist, check the valid time samples within the ROI latency
else
    valid_data_samples = data_samples; % If not exist, check all shared time samples for now
end

if polarity~=0 % If ROI polarity direction is defined, find the peak in the defined direction
    [~, peak_chan_id] = max(max(polarity*(avg_data(valid_data_channels,valid_data_samples)),[],2)); % Find the peak channel
    [~, peak_sample_id] = max(max(polarity*(avg_data(valid_data_channels(peak_chan_id),valid_data_samples)),[],1)); % Find the peak time sample
else
    [~, peak_chan_id] = max(max(abs(avg_data(valid_data_channels,valid_data_samples)),[],2)); % Find the peak channel
    [~, peak_sample_id] = max(max(abs(avg_data(valid_data_channels(peak_chan_id),valid_data_samples)),[],1)); % Find the peak time sample
end
peak_chan_value = avg_data(valid_data_channels(peak_chan_id),valid_data_samples(peak_sample_id));

zero_cross = find(diff(sign(avg_data(valid_data_channels(peak_chan_id),:)))~=0); % Find the zero-crossing time samples
zero_cross = unique([find(data.time{1,1}>=0,1,'first') zero_cross length(data.time{1,1})]); % Add the begin and end time samples
latency_samples = [ zero_cross( find( data.time{1,1}(zero_cross) <= (data.time{1,1}(valid_data_samples(peak_sample_id))) ,1,'last') ), ...
    zero_cross( find( data.time{1,1}(zero_cross) >= (data.time{1,1}(valid_data_samples(peak_sample_id))) ,1,'first') )]; % Find time samples around the peak latency until zero-crossing or time boundaries
latency = data.time{1,1}(latency_samples); % Default ROI latency range [begin end] time points (around the peak latency reaching zero-crossing or time boundaries)
same_sign_channels = channel_labels(sign(avg_data(valid_data_channels,valid_data_samples(peak_sample_id))) == sign(peak_chan_value)); % Default ROI channels (any channels with same polarity as the peak amplitude channel)

if ~isfield(cfg.roi,'channels') % Check if not exist
    cfg.roi.channels = same_sign_channels; % If not, insert default
end
if ~isfield(cfg.roi,'polarity') % Check if not exist
    cfg.roi.polarity = sign(peak_chan_value); % If not, insert default
end
if ~isfield(cfg.roi,'latency') % Check if not exist
    cfg.roi.latency = latency; % If not, insert default
end
if ~isfield(cfg,'hsac') % If no HSAC settings are defined, create a HSAC structure
    cfg.hsac = struct;
end
if ~isfield(cfg.hsac,'mc_iterations') % Check if not exist
    cfg.hsac.mc_iterations = 100; % If not, insert default
end
if ~isfield(cfg.hsac,'p_match') % Check if not exist
    cfg.hsac.p_match = .001; % If not, insert default
end
if ~isfield(cfg.hsac,'p_comp') % Check if not exist
    cfg.hsac.p_comp = .70; % If not, insert default
end
cfg.inspect             = ft_getopt(cfg, 'inspect', 'no');
cfg.layout              = ft_getopt(cfg, 'layout', []);


%% Visualize the data topographies and waveforms

if strcmp(cfg.inspect,'yes')
    
    fig_hsac = figure('name','SCA HSAC','color','w');
    
    % Original topography
    if ~isempty(cfg.layout)
        subplot(2,2,1)
        cfg_inspect = {};
        cfg_inspect.channel = channel_labels;
        cfg_inspect.xlim = [data.time{1,1}(find(data.time{1,1} >= cfg.roi.latency(1),1,'first')) data.time{1,1}(find(data.time{1,1} <= cfg.roi.latency(2),1,'last'))];
        cfg_inspect.colorbar = 'no';
        cfg_inspect.comment = ' ';
        cfg_inspect.baseline = 'no';
        cfg_inspect.layout = cfg.layout;
        cfg_inspect.zlim = 'maxabs';
        avg_data_temp = {};
        avg_data_temp.avg = avg_data(data_channels,:);
        avg_data_temp.time = data.time{1,1};
        avg_data_temp.label = data.label(data_channels);
        avg_data_temp.dimord = 'chan_time';
        evalc('ft_topoplotER(cfg_inspect, avg_data_temp)');
        title('Original data')
    else
        warning('No channel layout was provided! Skipping topography plotting.')
    end
    % Original waveform
    subplot(2,2,3)
    plot(data.time{1,1}*1000, avg_data(data_channels,:)) % butterfly plot
    hold on
    plot(data.time{1,1}*1000, zeros(length(avg_data),1),'black') % baseline
    xlim([time(1) time(end)]*1000)
    ylims = max(max(abs(avg_data(data_channels,data_samples))))*[-1 +1];
    if ylims(1)~=0
        ylim( ylims )
    end
    xlabel('ms')
    title('Original data')
    
end


%% Region of interest analysis

tic; % Reset time counter

% Find components with matching peak channels and polarity
cfg.roi.channels = ft_channelselection(cfg.roi.channels, channel_labels);
if cfg.roi.polarity == 0
    [~,peak_channel_id_pos] = max(comps.topo,[],1);
    [~,peak_channel_id_neg] = min(comps.topo,[],1);
    peak_channels_pos = comps.topolabel(peak_channel_id_pos);
    peak_channels_neg = comps.topolabel(peak_channel_id_neg);
    roi_match_peak_channel = find( ismember( peak_channels_pos, cfg.roi.channels) | ismember(peak_channels_neg, cfg.roi.channels) );
else
    [~,peak_channel_id] = max( cfg.roi.polarity*comps.topo ,[],1);
    peak_channels = comps.topolabel(peak_channel_id);
    roi_match_peak_channel = find(ismember(peak_channels, cfg.roi.channels));
end

% Find components with matching peak latencies
roi_match_peak_latency = find( comps.latency >= cfg.roi.latency(1)*1000 & comps.latency <= cfg.roi.latency(2)*1000 );

roi_match = intersect(roi_match_peak_channel, roi_match_peak_latency); % Matching peak channel, polarity, and latency

ptpc_roi = toc/length(comps.label); % Save processing time (s.) per component for ROI match

fprintf('ROI channels: Peak in')
for i=1:length(cfg.roi.channels)
    fprintf([' ', cfg.roi.channels{i}])
end
fprintf('.\n')
if cfg.roi.polarity==-1
    disp('ROI polarity: Negative peak.')
elseif cfg.roi.polarity==+1
    disp('ROI polarity: Positive peak.')
else
    disp('ROI polarity: Any negative or positive peak.')
end
disp(['ROI latency : ',num2str(cfg.roi.latency(1)*1000),' - ',num2str(cfg.roi.latency(2)*1000),' ms peak latency.'])

disp(['Excluded ',num2str(size(comps.trial{1,1},1)-length(roi_match)),' SCA components outside the region of interest.'])
disp(['Keeping ',num2str(length(roi_match)), ' SCA components for half-split average consistency analysis.'])


%% Monte Carlo simulation of averages for uniformly distributed half-split trial combinations

disp(['Monte Carlo simulation of ',num2str(cfg.hsac.mc_iterations),' half-split averages with uniformly distributed trial combinations.'])

tic; % Reset time counter

% Create set of uniformly distributed half-split trial combinations
data.hs_set = [];
for i=1:cfg.hsac.mc_iterations
    data.hs_set(i,:) = randperm(length(data.trial),round(length(data.trial)/2));
end

% Calculate averages for each half-split trial combination
for i=1:cfg.hsac.mc_iterations
    data.hs_avg{i} = zeros(size(data.trial{data.hs_set(i,1)}));
    for j=1:size(data.hs_set,2)
        data.hs_avg{i} = data.hs_avg{i} + (data.trial{data.hs_set(i,j)}/size(data.hs_set,2));
    end
end

% Calculate the average
avg = mean(cat(3,data.trial{:}),3);

if strcmp(cfg.inspect,'yes')
    
    figure('name','Half-split averages','color','w')
    
    % Show distribution of trial occurences across Monte Carlo simulations
    subplot(1,3,1)
    [N,~] = hist(data.hs_set(:), length(data.trial));
    bar(1:length(data.trial),N/cfg.hsac.mc_iterations*100)
    xlim([1 length(data.trial)])
    ylim([0 100])
    xlabel('trial #')
    ylabel('% across HSAs')
    title('Trial distributions across HSAs','interpreter','none')
    
    % Show Monte Carlo simulated half-split trial combinations
    channels = ismember(data.label,cfg.roi.channels);
    avg_sim = zeros(size(data.hs_avg{1}));
    for i=1:cfg.hsac.mc_iterations
        avg_sim = avg_sim + (data.hs_avg{i}/cfg.hsac.mc_iterations);
    end
    subplot(1,3,2)
    hold on
    for i=1:length(data.trial)
        plot(data.time{1}*1000, mean(data.trial{i}(channels,:),1),'r')
        if i==1
            plot(data.time{1}*1000, mean(avg(channels,:),1),'linewidth',3,'color','g')
            legend({'Trials','Average'},'location','northoutside')
        end
    end
    plot(data.time{1}*1000, mean(avg(channels,:),1),'linewidth',3,'color','g')
    y_lims = get(gca,'Ylim'); ylim(max(abs(y_lims))*[-1 1])
    xlim([time(1) time(end)]*1000)
    xlabel('ms')
    ylabel('Average across ROI channels')
    title('Trials')
    subplot(1,3,3)
    hold on
    for i=1:cfg.hsac.mc_iterations
        plot(data.time{1}*1000, mean(data.hs_avg{i}(channels,:),1),'black')
        if i==1
            plot(data.time{1}*1000, mean(avg(channels,:)),'linewidth',4,'color','g')
            legend({'HSAs','Average'},'location','northoutside')
        end
    end
    plot(data.time{1}*1000, mean(avg(channels,:)),'linewidth',4,'color','g')
    y_lims = get(gca,'Ylim'); ylim(max(abs(y_lims))*[-1 1])
    xlim([time(1) time(end)]*1000)
    xlabel('ms')
    ylabel('Average across ROI channels')
    title('HSAs')
    
end


%% Half-split average consistency

% Isolate SCA components which are consistently detected across half-split averages
disp(['Significance threshold for the correlation between an SCA component and HSA, p_match<.',num2str(cfg.hsac.p_match)])
disp(['Minimum probability that the SCA component signficantly correlates with HSAs, p_comp>=',num2str(cfg.hsac.p_comp)])
fprintf('Detecting SCA components which are consistent across half-split averages... ')
comp_r_wave(1:size(comps.trial{1,1},1), 1:cfg.hsac.mc_iterations) = 0; % Initially disqualify correlations for SCA components outside the region of interest
comp_p_wave(1:size(comps.trial{1,1},1), 1:cfg.hsac.mc_iterations) = 1; % Initially disqualify significance for SCA components outside the region of interest
comp_space = []; % Prepare component space waveforms

for i=1:length(roi_match) % Loop over SCA components with ROI match
    
    for j=1:cfg.hsac.mc_iterations % Loop over Monte Carlo simulations of half-split averages
        
        % Project the SCA component back to channel space
        component_projection = repmat(comps.topo(comps_channels,roi_match(i)), [1, length(comps_samples)]) .* repmat(comps.trial{1,1}(roi_match(i),comps_samples)', [1, length(comps_channels)])';
        
        % Transform the multi-channel waveforms into component space and remove all other SCA components except the tested SCA component
        comp_space(roi_match(i),j,:) = mean( ( data.hs_avg{j}(data_channels,data_samples) - ( avg(data_channels,data_samples) - component_projection ) ).*repmat(comps.topo(comps_channels,roi_match(i)), [1, length(comps_samples)]) ,1);
        
        % Calculate the Pearson correlation, r
        [R,P] = corrcoef(comps.trial{1,1}(roi_match(i),comps_samples)', squeeze(comp_space(roi_match(i),j,:)));
        comp_r_wave(roi_match(i),j) = R(1,2);
        comp_p_wave(roi_match(i),j) = P(1,2);
        
        % Disqualify anti-correlations
        if comp_r_wave(roi_match(i),j) < 0
            comp_r_wave(roi_match(i),j) = 0;
            comp_p_wave(roi_match(i),j) = 1;
        end
        
    end
    
end
fprintf('Done.\n')

ptpc_hsac = toc/length(roi_match); % Show processing time (s.) per component for HSAC
disp(['(Average processing time per SCA component is ',num2str(ptpc_hsac),' seconds.)'])

% Show component consistency across half-split averages
if strcmp(cfg.inspect,'yes')
    figure('name','SCA component consistency','color','w')
    subplot(2,1,1)
    imagesc(1-(double(comp_p_wave < cfg.hsac.p_match))',[0 1])
    set(gca,'Ydir','normal')
    xlim([0 size(comp_p_wave,1)])
    xlabel('Component #')
    ylabel('Sign. correlation with HSA #')
    subplot(2,1,2), bar( sum(comp_p_wave < cfg.hsac.p_match, 2) )
    xlim([0 size(comp_p_wave,1)])
    hold on
    line([1 size(comp_p_wave,1)],[cfg.hsac.p_comp cfg.hsac.p_comp]*cfg.hsac.mc_iterations)
    bar(sum(comp_p_wave < cfg.hsac.p_match,2).*(sum(comp_p_wave < cfg.hsac.p_match,2) < cfg.hsac.p_comp * cfg.hsac.mc_iterations),'r')
    xlabel('Component #')
    ylabel('% sign. correlations with HSAs')
    ylim([0 100])
    
end

consistent_comps = find( sum(comp_p_wave < cfg.hsac.p_match,2)/cfg.hsac.mc_iterations >= cfg.hsac.p_comp ); % Find consistent components


%% Extract SCA components which are consistently matched with half-split averages

disp(['Projecting ',num2str(length(consistent_comps)),' consistent SCA components into channel space. '])

% Prepare the output
data_hsac = {};
data_hsac.time = comps.time{1,1};
data_hsac.fsample = 1/diff(comps.time{1,1}(1:2));
data_hsac.avg = zeros( length(comps.topolabel) , length(comps.time{1,1}) );
data_hsac.label = comps.topolabel;
data_hsac.dimord = comps.dimord;
data_hsac.extracted_sca_comps = consistent_comps;

for i=1:length(consistent_comps) % loop across consistent SCA components within ROI
    
    component_projection = repmat(comps.topo(:,consistent_comps(i)),1,size(data_hsac.avg,2)).*repmat(comps.trial{1,1}(consistent_comps(i),:)',1,size(data_hsac.avg,1))'; % Project the component back into channel space
    
    data_hsac.avg = data_hsac.avg + component_projection;
    
end


%% Show fits between half-split averages and consistent SCA components

if strcmp(cfg.inspect,'yes')
    
    figure('name','consistent comps','color','w')
    
    n_plots = ceil(sqrt(length(consistent_comps)));
    
    for i=1:length(consistent_comps) % loop across consistent SCA components within ROI
        
        subplot(n_plots,n_plots,i)
        hold on
        plot(comps.time{1,1}(comps_samples)*1000, squeeze(comp_space(consistent_comps(i),:,:)),'black')
        plot(comps.time{1,1}(comps_samples)*1000, mean( repmat(comps.trial{1,1}(consistent_comps(i),comps_samples)', [1, length(comps_channels)])' .* repmat(comps.topo(comps_channels,consistent_comps(i)).^2, [1, length(comps.trial{1,1}(consistent_comps(i),comps_samples))]), 1), 'g','linewidth',4)
        if (sum(comp_p_wave(consistent_comps(i),:) >= cfg.hsac.p_match)>0)
            plot(comps.time{1,1}*1000, squeeze(comp_space(consistent_comps(i),comp_p_wave(consistent_comps(i),:) >= cfg.hsac.p_match,:)),'r');
        end
        
        title(['SCA component #',num2str(consistent_comps(i)),' and HSAs'])
        xlabel('ms')
        ylabel('Component space (a.u.)')
    end
    
end

%% Visualize the extracted SCA components topography and waveform

if strcmp(cfg.inspect,'yes')
    
    figure(fig_hsac)
    
    % Extracted SCA components topography
    if ~isempty(cfg.layout)
        subplot(2,2,2)
        cfg_inspect = {};
        cfg_inspect.channel = channel_labels;
        cfg_inspect.xlim = [data_hsac.time(comps_samples(1)) data_hsac.time(comps_samples(end))];
        cfg_inspect.colorbar = 'no';
        cfg_inspect.comment = ' ';
        cfg_inspect.baseline = 'no';
        cfg_inspect.layout = cfg.layout;
        cfg_inspect.zlim = 'maxabs';
        evalc('ft_topoplotER(cfg_inspect, data_hsac)');
        title('Extracted SCA components')
    else
        warning('No channel layout was provided! Skipping topography plotting.')
    end
    
    % Extracted SCA components waveform
    subplot(2,2,4)
    plot(data_hsac.time*1000, data_hsac.avg(comps_channels,:)) % butterfly plot
    hold on
    plot(data_hsac.time*1000, zeros(length(data_hsac.time),1),'black') % baseline
    xlim([ max([data.time{1,1}(1), comps.time{1,1}(1)]) min([data.time{1,1}(end), comps.time{1,1}(end)]) ]*1000)
    if ylims(1)~=0
        ylim( ylims )
    end
    xlabel('ms')
    title('Extracted SCA components')
    
end

end
