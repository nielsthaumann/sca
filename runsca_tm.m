function data_tm = runsca_tm(cfg, comps, template)

% Spike density components with template match.
% Projects SCA components (output of the RUNSCA function) back into channel space,
% if the SCA components matches the topography and morphology of a template
% (e.g., grand-average data).
% The output data can be visualized and processed with FieldTrip functions.
%
% Use as
%
%   data_tm = runsca_tm(cfg, comps, template)
%
% where cfg is a configuration structure,
% comps are SCA components obtained with RUNSCA,
% and template is average evoked response MEG or EEG waveforms
% obtained with FT_TIMELOCKANALYSIS or FT_TIMELOCKGRANDAVERAGE,
%
%
% The configuration or part of the configuration can simply be empty (e.g., cfg = []),
% in which case default settings are applied (see below).
%
% Settings
%
% cfg.latency    = template latency range in seconds [begin end] (e.g., [0.100 0.200] ) (default = all time points around the peak latency until zero-crossing or 0 or end time boundaries)
% cfg.inspect    = visually inspect the template matched topographies and waveforms ('yes' or 'no' (default))
% cfg.layout     = is needed for inspecting the topographies. Use the output of FT_PREPARE_LAYOUT as input
%
% This function runs in the Matlab environment and requires the FieldTrip toolbox to be installed (see https://github.com/fieldtrip).
%
% For more information, see:
%
% Haumann, N T; Petersen, B; Friis Andersen, A S; Faulkner, K S; Brattico, E; Vuust, P;
% "Mismatch negativity as a marker of music perception in individual cochlear implant users: 
% A spike density component analysis study",
% Clinical Neurophysiology (2023), https://doi.org/10.1016/j.clinph.2023.01.015
% 
% Bruzzone, S E P; Haumann, N T; Kliuchko, M; Vuust, P, Brattico, E;
% "Applying Spike-density Component Analysis for high-accuracy auditory event-related potentials in children",
% Clinical Neurophysiology (2021), https://doi.org/10.1016/j.clinph.2021.05.007
%
% Haumann, N T; Hansen, B; Huotilainen, M; Vuust, P; Brattico, E;
% "Applying Stochastic Spike train theory for high-accuracy human MEG/EEG"
% Journal of Neuroscience Methods (2020), doi: https://doi.org/10.1016/j.jneumeth.2020.108743
%

%% Prepare settings

warning('It is recommended applying RUNSCA_STATS instead of RUNSCA_TM (see Haumann et al. (2023), "Mismatch negativity as a marker of music perception in individual cochlear implant users: A spike density component analysis study", Clinical Neurophysiology, https://doi.org/10.1016/j.clinph.2023.01.015).')

% Verify required inputs are provided
if nargin<3
    error('Please provide the required inputs (cfg, comps, template). Type help runsca_tm for more information.')
end

% Verify second input is SCA components
if ~isfield(comps,'sigma')
    error('Please provide SCA components as the second input. Type help runsca_tm for more information.')
end

% Verify FieldTrip is installed
if ~exist('ft_getopt','file')
    error('Could not find the function ''ft_getopt''. Please ensure that the FieldTrip Toolbox is installed, and related functions are added to the paths with ''ft_defaults''.')
end

% Ensure the sampling rate of comps and template is the same
if 1/diff(comps.time{1,1}(1:2)) ~= 1/diff(template.time(1:2))
    error(['The sampling rate of the template (',num2str(1/diff(template.time(1:2))),' Hz) and comps (',num2str(1/diff(comps.time{1,1}(1:2))),' Hz) must be the same.'])
end

% Ensure to process only the shared comp and template channels
channel_labels = intersect(comps.topolabel, template.label); % Shared channel labels
if isempty(channel_labels)
    error('The template and comps must have shared channels for the matching.')
end
if length(comps.topolabel) ~= length(template.label)
    warning(['The template and comps have different channels. Using only the ',num2str(length(intersect(comps.topolabel, template.label))),' shared channels.'])
end

% Correct the order of channel indices for matching between the template and comps channels
comps_channels = [];
for i=1:length(channel_labels)
    comps_channels(i) = find(ismember(comps.topolabel,channel_labels{i}));
end
template_channels = [];
for i=1:length(channel_labels)
    template_channels(i) = find(ismember(template.label,channel_labels{i}));
end

% Ensure to process only the shared comps and template time samples
time = intersect(round(comps.time{1,1}*1000), round(template.time*1000))/1000; % Shared time points
if isempty(time)
    error('The comps and template must have shared time points.')
end
if length(comps.time{1,1}) ~= length(template.time)
    warning(['The comps and template have different time points. Using only the ',num2str(length(time)),' shared time points.'])
end

% Correct the time sample indices for matching between the comps and data time points
comps_samples = find(ismember(round(comps.time{1,1}*1000), round(time*1000)));
template_samples = find(ismember(round(template.time*1000), round(time*1000)));
percent_consecutive_time = 100 - ( sum( diff(time) >=  1.5*diff(comps.time{1,1}(1:2)) )/length(time)*100 ); % Check whether any time steps are non-consecutive (>=1.5 sampling interval)
if percent_consecutive_time ~= 100
    warning([num2str(100 - percent_consecutive_time),'% of the shared time points are non-consecutive. There is up to ',num2str( max(abs(diff(time) - diff(comps.time{1,1}(1:2)))) * 1000),' ms jitter in the time steps within the applied time range. (This issue might be related to diverging preprossesing procedures for the inputs.)'])
end


% If no configuration is provided, apply the default settings

cfg.latency               = ft_getopt(cfg, 'latency', [0 time(end)]);
cfg.inspect               = ft_getopt(cfg, 'inspect', 'no');
cfg.layout                = ft_getopt(cfg, 'layout', []);


%% Find valid time samples for the template waveform until zero-crossing

sample_range = find( template.time >= cfg.latency(1) & template.time <= cfg.latency(2) ); % Time samples range
[~,peak_id] = max(max(abs( template.avg(:, sample_range)),[],1)); % Find the peak time within the time samples range
peak_sample = sample_range(peak_id);
peak_time = template.time(peak_sample);
[~,peak_channel] = max(abs(template.avg(:, peak_sample))); % Find the template peak channel
template_latency_samples = find(template.time >= cfg.latency(1),1,'first'):find(template.time <= cfg.latency(2),1,'last'); % Find the template latency samples
zero_cross = find(diff(sign(template.avg(peak_channel, template_latency_samples )))~=0); % Find the template zero-crossing time samples
zero_cross = template_latency_samples(zero_cross);
zero_cross = unique([template_latency_samples(1) zero_cross template_latency_samples(end)]); % Add the start and end time samples
valid_template_samples = zero_cross( find(template.time(zero_cross) <= peak_time,1,'last') ) : zero_cross( find(template.time(zero_cross) >= peak_time,1,'first') ); % Find all time sample around the peak latency until zero-crossing or latency range time boundaries
valid_time = intersect(round(time*1000),round(template.time(valid_template_samples)*1000))/1000;
valid_comps_samples = find(ismember(round(comps.time{1,1}*1000),round(valid_time*1000)));

if length(valid_time) < 2
    error('Found less than 2 valid time samples to compare between the template and comps. Consider changing the latency constraints. Type help runsca_tm for more information.')
end


%% Visualize the data and template topographies and waveforms

if strcmp(cfg.inspect,'yes')
    
    fig_match = figure('name','SCA template match','color','w');
    
    % Template topography
    if ~isempty(cfg.layout)
        subplot(2,3,1)
        cfg_inspect = {};
        cfg_inspect.channel = channel_labels;
        cfg_inspect.xlim = [template.time(valid_template_samples(1)) template.time(valid_template_samples(end))];
        cfg_inspect.colorbar = 'no';
        cfg_inspect.comment = ' ';
        cfg_inspect.baseline = 'no';
        cfg_inspect.layout = cfg.layout;
        cfg_inspect.zlim = 'maxabs';
        evalc('ft_topoplotER(cfg_inspect, template)');
        title('Template')
    else
        warning('No channel layout was provided! Skipping topography plotting.')
    end
    
    % Template waveform
    subplot(2,3,4)
    plot(template.time(template_samples)*1000, template.avg(template_channels,template_samples)) % butterfly plot
    hold on
    plot(template.time*1000, zeros(length(template.time),1),'black') % baseline
    xlim([ max([template.time(1), comps.time{1,1}(1)]) min([template.time(end), comps.time{1,1}(end)]) ]*1000)
    ylim( max(max(abs(template.avg(template_channels,template_samples))))*[-1 +1] )
    xlabel('ms')
    title('Template')
    
    % Project all SCA components back to channel space
    data_all_sca = {};
    data_all_sca.time = comps.time{1,1};
    data_all_sca.fsample = 1/diff(comps.time{1,1}(1:2));
    data_all_sca.avg = zeros( length(comps.topolabel) , length(comps.time{1,1}) );
    data_all_sca.label = comps.topolabel;
    data_all_sca.dimord = comps.dimord;
    for i=1:size(comps.trial{1,1},1) % Add each SCA component
        component_projection = repmat(comps.topo(:,i),1,size(data_all_sca.avg,2)).*repmat(comps.trial{1,1}(i,:)',1,size(data_all_sca.avg,1))'; % Project the component back into channel space
        data_all_sca.avg = data_all_sca.avg + component_projection;
    end
    
    % All SCA components topography
    if ~isempty(cfg.layout)
        subplot(2,3,2)
        cfg_inspect = {};
        cfg_inspect.channel = channel_labels;
        cfg_inspect.xlim = [data_all_sca.time(valid_comps_samples(1)) data_all_sca.time(valid_comps_samples(end))];
        cfg_inspect.colorbar = 'no';
        cfg_inspect.comment = ' ';
        cfg_inspect.baseline = 'no';
        cfg_inspect.layout = cfg.layout;
        cfg_inspect.zlim = 'maxabs';
        evalc('ft_topoplotER(cfg_inspect, data_all_sca)');
        title('All SCA components')
    else
        warning('No channel layout was provided! Skipping topography plotting.')
    end
    
    % All SCA components waveform
    subplot(2,3,5)
    plot(data_all_sca.time*1000, data_all_sca.avg(comps_channels,:)) % butterfly plot
    hold on
    plot(data_all_sca.time*1000, zeros(length(data_all_sca.time),1),'black') % baseline
    xlim([ max([template.time(1), comps.time{1,1}(1)]) min([template.time(end), comps.time{1,1}(end)]) ]*1000)
    xlabel('ms')
    ylims = max(max(abs(data_all_sca.avg(comps_channels,comps_samples))))*[-1 +1];
    if ylims(1)~=0
        ylim( ylims )
    end
    title('All SCA components')
    
end


%% Find SCA components that matches the template

disp(['Estimating match between the template and ',num2str(size(comps.trial{1,1},1)),' SCA components with Pearson correlations '])
disp(['based on ',num2str(length(channel_labels)),' shared channels, '])
disp([num2str(template.time(peak_sample)*1000),' ms peak in the template, and ',num2str(valid_time(1)*1000),' - ',num2str(valid_time(end)*1000),' ms shared latency time range.' ])

comp_r_topo = []; comp_p_topo = []; comp_r_wave = []; comp_p_wave = []; % Initialize variables

% Estimate the template topography for mathing with the components
template_topo = mean(template.avg(template_channels,valid_template_samples),2);

for i=1:size(comps.topo,2) % Loop over SCA components
    
    % Topography correlation
    [R,P] = corrcoef( comps.topo(comps_channels,i)*mean(comps.trial{1,1}(i,valid_comps_samples),2) , template_topo); % Topography correlation
    comp_r_topo(i) = R(1,2);
    comp_p_topo(i)= P(1,2);
    
    % Waveform correlation
    r_channel = []; p_channel = [];
    for j=1:size(comps.topo,1) % Loop over channels
        
        [R,P] = corrcoef(comps.trial{1,1}(i,valid_comps_samples)'*comps.topo(comps_channels(j),1), template.avg(template_channels(j), valid_template_samples)'); % Waveform correlation for each channel
        r_channel(i) = R(1,2);
        p_channel(i)= P(1,2);
        
    end
    comp_r_wave(i) = mean(r_channel); % Mean correlation r of waveform across channels
    comp_p_wave(i) = mean(p_channel); % Mean correlation p of waveform across channels
    
end

corr_sign = (1-2*( sign(comp_r_topo)<1 | sign(comp_r_topo)<1) ); % Find correlation sign for each component (must be positive for both topogaphy and waveform correlations)
match_r = corr_sign.*abs(comp_r_topo).*abs(comp_r_wave); % Calculate signed pseudo-r-squared estimates combining the topography and waveform matches for each component
exclude = corr_sign' == -1 | comps.latency < (template.time(valid_template_samples(1))*1000) | comps.latency > (template.time(valid_template_samples(end))*1000); % Exclude components with negative topography or waveform correlation or peaking outside the latency range
match_r(exclude) = 0; % Disqualify components with negative correlations to the template or peaking outside latency range
[values,match_comps] = sort(match_r,'descend'); % Arrange the pseudo-r-squared for the components from highest to lowest values and obtain the component numbers
match_comps = match_comps(values>0); % Remove the excluded component numbers


%% Extract matching components

disp('Extracting SCA components matching the template, from highest to lowest match,')
disp('as long as the average waveform correlation across channels increases.')

r_change = [];
mean_r_change = [];

data_tm = {};
data_tm.time = comps.time{1,1};
data_tm.fsample = 1/diff(comps.time{1,1}(1:2));
data_tm.avg = zeros( length(comps.topolabel) , length(comps.time{1,1}) );
data_tm.label = comps.topolabel;
data_tm.dimord = comps.dimord;
data_tm.extracted_sca_comps = [];

for i=1:length(match_comps) % Loop over components, from highest to lowest match
    
    component_projection = repmat(comps.topo(:,match_comps(i)),1,size(data_tm.avg,2)).*repmat(comps.trial{1,1}(match_comps(i),:)',1,size(data_tm.avg,1))'; % Project the component back into channel space
    
    % Check whether inclusion of the component increases the covariance
    % between the component and template
    for j=1:length(channel_labels) % Loop over shared channels
        if nansum(nansum(abs(data_tm.avg))) > 0
            R = corrcoef(data_tm.avg(comps_channels(j),valid_comps_samples)'+component_projection(comps_channels(j),valid_comps_samples)',template.avg(template_channels(j),valid_template_samples)') - corrcoef(data_tm.avg(j,valid_comps_samples)',template.avg(j,valid_template_samples)'); % Calculate r-change for each channel
            r_change(j) = R(1,2);
        else
            R = corrcoef(data_tm.avg(comps_channels(j),valid_comps_samples)'+component_projection(comps_channels(j),valid_comps_samples)',template.avg(template_channels(j),valid_template_samples)'); % Calculate r-change for each channel
            r_change(j) = R(1,2);
        end
    end
    
    mean_r_change(i) = mean(r_change);
    
    if mean_r_change(i) > 0 % If covarance increases by adding the component, add the component
        data_tm.avg = data_tm.avg + component_projection;
        data_tm.extracted_sca_comps = vertcat(data_tm.extracted_sca_comps,match_comps(i));
        
    else % else stop adding more components
        break
    end
    
end

disp(['Projected ',num2str(length(data_tm.extracted_sca_comps)),' matching SCA components into channel space.'])
disp(['The estimated match between the SCA components and template is ',num2str(round(sum(mean_r_change(1:end-1))*100)),'%.'])


%% Visualize the extracted SCA components topography and waveform

if strcmp(cfg.inspect,'yes')
    
    figure(fig_match)
    
    % Extracted SCA components topography
    if ~isempty(cfg.layout)
        subplot(2,3,3)
        cfg_inspect = {};
        cfg_inspect.channel = channel_labels;
        cfg_inspect.xlim = [data_tm.time(valid_comps_samples(1)) data_tm.time(valid_comps_samples(end))];
        cfg_inspect.colorbar = 'no';
        cfg_inspect.comment = ' ';
        cfg_inspect.baseline = 'no';
        cfg_inspect.layout = cfg.layout;
        cfg_inspect.zlim = 'maxabs';
        evalc('ft_topoplotER(cfg_inspect, data_tm)');
        title('Extracted SCA components')
    else
        warning('No channel layout was provided! Skipping topography plotting.')
    end
    
    % Extracted SCA components waveform
    subplot(2,3,6)
    plot(data_tm.time*1000, data_tm.avg(comps_channels,:)) % butterfly plot
    hold on
    plot(data_tm.time*1000, zeros(length(data_tm.time),1),'black') % baseline
    xlim([ max([template.time(1), comps.time{1,1}(1)]) min([template.time(end), comps.time{1,1}(end)]) ]*1000)
    xlabel('ms')
    if ylims(1)~=0
        ylim( ylims )
    end
    title('Extracted SCA components')
    
end

end
