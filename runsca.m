function comp = runsca(cfg, data)

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

%% Prepare settings

if ~exist('ft_getopt','file')
    error('Could not find the function ''ft_getopt''. Please ensure that the FieldTrip Toolbox is installed, and related functions are added to the paths with ''ft_defaults''.')
end

% If no configuration is provided, apply the default settings

cfg.channel             = ft_getopt(cfg, 'channel',             'all');
cfg.search_time         = ft_getopt(cfg, 'search_time',         [data.time(1) data.time(end)]);
cfg.baseline_correct    = ft_getopt(cfg, 'baseline_correct',    'mean');
cfg.baseline_window     = ft_getopt(cfg, 'baseline_window',     [data.time(1) 0]);
cfg.model               = ft_getopt(cfg, 'model',               'gauss');
cfg.amplitude_threshold = ft_getopt(cfg, 'amplitude_threshold', 0);
cfg.sort_by             = ft_getopt(cfg, 'sort_by',             'latency');
cfg.spatial_correlation = ft_getopt(cfg, 'spatial_correlation', 1);
cfg.show_components     = ft_getopt(cfg, 'show_components',     'no');
cfg.show_residuals      = ft_getopt(cfg, 'show_residuals',      'yes');
cfg.demomode            = ft_getopt(cfg, 'demomode',            'off');


%% Channel selection

chanid = ismember(data.label, ft_channelselection(cfg.channel, data.label)); 
fprintf('Selecting %d channels...\n', sum(chanid))
data.avg = data.avg(chanid,:);
data.label = data.label(chanid);
if sum(chanid)==0
  error('No channels were selected');
end

% Avoid issues related to spaces in either data labels or layout labels
if isfield(cfg,'layout')
    for i=1:length(data.label)
        data.label{i}(regexp(data.label{i},' ')) = [];
    end
    for i=1:length(cfg.layout.label)
        cfg.layout.label{i}(regexp(cfg.layout.label{i},' ')) = [];
    end
    if sum(ismember(cfg.layout.label,data.label)) < length(data.label)
        error('Channel labels in the data are missing in the topography layout.')
    end
end


%% Scaling and unit settings

if ~isfield(cfg,'scaling_factor')
    
    peak_amplitude = max(max(abs(data.avg)));
    
    if peak_amplitude > 10^-1
        if strcmp(upper(cfg.channel),'EEG')
            cfg.scaling_unit = '\muV';
        elseif strcmp(upper(cfg.channel),'MEG') || strcmp(upper(cfg.channel),'MEGMAG') || strcmp(upper(cfg.channel),'MEGGRAD')
            cfg.scaling_unit = 'fT';
        else
            cfg.scaling_unit = 'unknown unit';
        end
        cfg.scaling_factor = 1;
        
    elseif peak_amplitude <= 10^-1 && peak_amplitude > 10^-7
        cfg.scaling_unit = '\muV';
        cfg.scaling_factor = 10^6;
    elseif peak_amplitude <= 10^-10 && peak_amplitude > 10^-16 && ( strcmp(upper(cfg.channel),'MEG') || strcmp(upper(cfg.channel),'MEGMAG') )
        cfg.scaling_unit = 'fT';
        cfg.scaling_factor = 10^15;
    elseif peak_amplitude <= 10^-10 && peak_amplitude > 10^-16 && ( strcmp(upper(cfg.channel),'MEGGRAD') )
        cfg.scaling_unit = 'fT/cm';
        cfg.scaling_factor = 10^13;
    else
        cfg.scaling_factor = 1;
    end
    
    warning(['ASSUMING THAT THE MEG/EEG DATA IS MEASURED IN ',cfg.scaling_unit,' scaled at ', num2str(1/cfg.scaling_factor),' !'])
    pause(3)
    
end
if ~isfield(cfg,'scaling_unit')
    cfg.scaling_unit = 'unknown unit';
end


%% Baseline correction

if strcmp(cfg.baseline_correct,'mean')
    fprintf('\nCorrecting baseline based on mean across %f to %f seconds.\n',cfg.baseline_window(1),cfg.baseline_window(2))
    data.avg = data.avg - repmat(mean(data.avg(:,data.time >= cfg.baseline_window(1) & data.time <= cfg.baseline_window(2)),2),1,size(data.avg,2));

elseif strcmp(cfg.baseline_correct,'median')
    fprintf('\nCorrecting baseline based on median across %f to %f seconds.\n',cfg.baseline_window(1),cfg.baseline_window(2))
    data.avg = data.avg - repmat(median(data.avg(:,data.time >= cfg.baseline_window(1) & data.time <= cfg.baseline_window(2)),2),1,size(data.avg,2));
   
end


%% Estimate spike density components

if (strcmp(cfg.model,'gauss') || strcmp(cfg.model,'sine')) && ~license('test','Curve_Fitting_Toolbox')
    warning('''Runsca'' requires the Curve Fitting Toolbox to be installed for Gaussian or sine models. ')
end
if strcmp(cfg.model,'gamma') && ~license('test','Statistics_Toolbox')
    warning('For experimental comparison with gamma models, the Statistics and Machine Learning Toolbox needs to be installed. ')
end

fprintf('\nRunning spike density component analysis...\n\n')

if strfind(cfg.scaling_unit,'\mu')==1
    mu_index = find(ismember(cfg.scaling_unit,'\mu'));
    cfg.scaling_unit(mu_index(end)) = char(181);
    cfg.scaling_unit(mu_index(1:2)) = [];
end

% (Visualization for demo mode)
if strcmp(cfg.demomode,'on')
    demo_figure = [];
    if ishandle(demo_figure)
        clf(demo_figure); % Clear the previous demo figure
    else
        demo_figure = figure('Color','White');
        screensize = get( groot, 'Screensize' ); % Obtain screen size information
        x0=screensize(1);y0=screensize(2);width=screensize(3);height=screensize(4); % Increase the figure size to fill the whole screen
        set(gcf,'units','points','position',[x0,y0,width,height]);     
    end
end

% Repeat until residuals increase or all data is modeled
tic;
time_samples = find(data.time >= cfg.search_time(1) & data.time <= cfg.search_time(2));
data.sampleinfo = [1 length(time_samples)];
data.time = data.time(time_samples);
data.avg = data.avg(:,time_samples);
residuals = data;
if isfield(residuals,'trial')
    residuals = rmfield(residuals,'time');
    residuals = rmfield(residuals,'trial');
    residuals.time = data.time{1,1};
    residuals.avg = data.trial{1,1};
end
initial_residuals = nansum(nansum(abs(residuals.avg)));
residual_increase = false;
i = 1;
component = [];
peak_amplitude = inf;

while residual_increase==false && nansum(nansum(abs(residuals.avg))) > 0 && peak_amplitude > cfg.amplitude_threshold
    
    % Calculate the explained variance (r-squared) across channels
    if i==1
        component.rsq(i) = 0; 
    else
        component.rsq(i-1) = nanmean(1 - nansum(residuals.avg.^2,2)./nansum((data.avg).^2,2)) - sum(component.rsq);
    end
    
    % (Visualization for demo mode)
    if strcmp(cfg.demomode,'on')
        if i==1
            subplot(5,4,1)
            plot(residuals.time*1000, residuals.avg*cfg.scaling_factor)
            title(['Step ',num2str(i),'. Data (',num2str(sum(component.rsq)*100),'% variance explained)'])
        elseif i<=5
            subplot(5,4,1+((i-1)*4))
            plot(residuals.time*1000, residuals.avg*cfg.scaling_factor)
            title(['Step ',num2str(i),'. Residuals (',num2str(sum(component.rsq)*100),'% variance explained)'])
        else
            subplot(5,4,1+((5-1)*4))
            cla
            plot(residuals.time*1000, residuals.avg*cfg.scaling_factor)
            title(['Step ',num2str(i),'. Residuals (',num2str(sum(component.rsq)*100),'% variance explained)'])
        end
        xlabel('ms')
        ylabel(cfg.scaling_unit)
        if i==1
            limits = axis;
        end
        xlim([min(residuals.time*1000) max(residuals.time*1000)])
        ylim(limits(3:4))
        if i<=5
           pause(1) 
        end
    end
    
    
    % Find the peak amplitude
    
    [peak_amplitude,peak_sample] = max(max(abs(residuals.avg),[],1),[],2);
    [~,peak_channel] = max(abs(residuals.avg(:,peak_sample)),[],1);
    
    fprintf(['Step ',num2str(i),'. Explained variance = ',num2str(sum(component.rsq)*100),' percent. Residual peak amplitude ',num2str(nansum(residuals.avg(peak_channel,peak_sample))*cfg.scaling_factor),' ',cfg.scaling_unit,'\n'])
    
    % Find the ends of the sorrounding valleys or nearest baseline crossing (bc)

    decrease_samples = sign([0 diff(residuals.avg(peak_channel,:))])==-sign(residuals.avg(peak_channel,peak_sample));
    index = find(decrease_samples(1:peak_sample),1,'last');
    if isempty(index)
        valley_sample(1) = 1;
    else
        valley_sample(1) = index;
    end
    index = find(~decrease_samples(peak_sample+1:end),1,'first');
    if isempty(index)
        valley_sample(2) = size(residuals.avg,2);
    else
        valley_sample(2) = index + peak_sample -1;
    end
    bc_samples = sign(residuals.avg(peak_channel,:))==-sign(residuals.avg(peak_channel,peak_sample));
    index = find(bc_samples(1:peak_sample-1),1,'last');
    if isempty(index)
        bc_sample(1) = 1;
    else
        bc_sample(1) = index + 1;
    end
    index = find(bc_samples(peak_sample+1:end),1,'first');
    if isempty(index)
        bc_sample(2) = size(residuals.avg,2);
    else
        bc_sample(2) = index + peak_sample - 1;
    end
    component_sample(1) = max([valley_sample(1),bc_sample(1)]);
    component_sample(2) = min([valley_sample(2),bc_sample(2)]);
    component_samples = component_sample(1):component_sample(2); 
    
    residuals_time = residuals.time(component_samples);
    residuals_curve = residuals.avg(peak_channel,component_samples)*cfg.scaling_factor;
    
    % (Visualization for demo mode)
    if strcmp(cfg.demomode,'on')
        if i<=5
            subplot(5,4,2+(i-1)*4)
        else
            subplot(5,4,2+(5-1)*4)
            cla
        end
        rectangle_time(1) = residuals.time(1)*1000;
        rectangle_time(2) = residuals.time(min(component_samples))*1000;
        rectangle_time(3) = residuals.time(max(component_samples))*1000;
        rectangle_time(4) = residuals.time(end)*1000;
        hold on
        if rectangle_time(2)-rectangle_time(1) > 0
            rectangle('Position',[rectangle_time(1) limits(3) rectangle_time(2)-rectangle_time(1) limits(4)-limits(3)],'FaceColor',[.7 .7 .7])
        end
        if rectangle_time(4)-rectangle_time(3) > 0
            rectangle('Position',[rectangle_time(3) limits(3) rectangle_time(4)-rectangle_time(3) limits(4)-limits(3)],'FaceColor',[.7 .7 .7])
        end
        plot(residuals.time*1000,residuals.avg(peak_channel,:)*cfg.scaling_factor,'blue')
        ylim(limits(3:4))
        xlim([min(residuals.time*1000) max(residuals.time*1000)])
        xlabel('ms')
        ylabel(cfg.scaling_unit)
        if strcmp(cfg.model,'gauss')
            title('Gaussian fit to peak signal ~ SNIR > 1')
        elseif strcmp(cfg.model,'sine')
            title('Sine arc fit to peak signal ~ SNIR > 1')
        elseif strcmp(cfg.model,'gamma')
            title('Gamma fit to peak signal ~ SNIR > 1')
        end
        if i<=5
           pause(1) 
        end
    end
    
    
    % Find the best model fit or if fails use raw curve (NB: scaling is important for precision)

    try 

        if strcmp(cfg.model,'gauss') % Fit gaussian distribution parameters
            f = fit(residuals_time.',residuals_curve.','gauss1');
            coeffvals = coeffvalues(f); % (1) alpha, (2) mu, (3) sigma
            a1 = coeffvals(1); b1 = coeffvals(2); c1 = coeffvals(3);
            component.waveform(i,:) = a1*exp(-((residuals.time-b1)/c1).^2)/cfg.scaling_factor;
            component.latency(i) = round(b1*1000); % Store component latency in ms
            component.channel{i} = data.label{peak_channel}; % Store component peak channel
            component.amplitude(i) = a1; % Store component amplitude
            component.sigma(i) = c1; % Store component sigma (spread or timing uncertainty)
            component.shape_param(i) = nan; % Store component shape parameter
            component.fit_success(i) = true; 

        elseif strcmp(cfg.model,'sine') % Fit sine arc parameters
            f = fit(residuals_time',residuals_curve.','sin1');
            coeffvals = coeffvalues(f); % (1) amplitude, (2) f(Hz), (3) phase shift
            a1 = coeffvals(1); b1 = coeffvals(2); c1 = coeffvals(3);
            component.waveform(i,:) = a1*sin(b1*residuals.time+c1)/cfg.scaling_factor;
            cut_samples = find(sign(component.waveform(i,:)) ~= sign(component.waveform(i,peak_sample)));
            cut_sample(1) = max(cut_samples(cut_samples < peak_sample))+1;
            cut_sample(2) = min(cut_samples(cut_samples > peak_sample))-1;
            component.waveform(i,~ismember(1:length(residuals.time),cut_sample(1):cut_sample(2))) = 0;
            [~,latency_sample] = max(abs(component.waveform));
            component.latency(i) = round(residuals.time(latency_sample)*1000); % Store component latency in ms
            component.channel{i} = data.label{peak_channel}; % Store component peak channel
            component.amplitude(i) = a1; % Store component amplitude
            component.sigma(i) = nan; % Store component sigma (spread or timing uncertainty)
            component.shape_param(i) = nan; % Store component shape parameter
            component.fit_success(i) = true; 
            
        elseif strcmp(cfg.model,'gamma') % Fit gamma distribution parameters
            
            % Initial alpha = adjustment * peak amplitude, initial shape is set to regular value of 10, 
            % and since shape * scaling parameter equals the mean, initial scaling
            % parameter is set to 0.1*peak latency.
            initial_pars = [];
            initial_pars(1) = 0.1*abs(residuals.avg(peak_channel,peak_sample))*cfg.scaling_factor; 
            initial_pars(2) = 10; 
            initial_pars(3) = 0.1*residuals.time(peak_sample); 
            
            f = nlinfit(residuals_time',residuals_curve',@gamma_function,initial_pars);
            
            coeffvals = f; % (1) alpha, (2) k (shape), (3) theta (scaling)
            a1 = coeffvals(1); b1 = coeffvals(2); c1 = coeffvals(3);
            component.waveform(i,:) = a1 * gampdf(residuals.time,b1,c1)/cfg.scaling_factor;
            if sum(isnan(component.waveform(i,:)))>0 % Check for issue in the gamma fit causing NaNs in modelled waveform
                error('Some of the modelled data does not contain real numbers.')
            end
            [~,latency_sample] = max(abs(component.waveform));
            component.latency(i) = round(residuals.time(latency_sample)*1000); % Store component latency in ms
            component.channel{i} = data.label{peak_channel}; % Store component peak channel
            [component.amplitude(i)] = component.waveform(i,latency_sample)*cfg.scaling_factor; % Store component amplitude
            component.sigma(i) = nan; % Store component sigma (spread or timing uncertainty)
            component.shape_param(i) = b1; % Store component shape parameter
            component.fit_success(i) = true; 

        end

    catch err % if fails use baseline corrected raw curve
        warning(err.message)
        disp('Modeling failed. Applying raw curve')
        component.waveform(i,:) = [zeros(1,component_sample(1)-1), residuals_curve - min(sign(residuals.avg(peak_channel,peak_sample))*[residuals_curve(1),residuals_curve(end)]), zeros(1,length(residuals.time)-length(component_samples)-component_sample(1)+1)]; % Extract waveform from raw curve
        component.waveform(i,:) = component.waveform(i,:)/cfg.scaling_factor; % Revert scaling
        component.latency(i) = round(data.time(peak_sample)*1000); % Store component latency in ms
        component.channel{i} = data.label{peak_channel}; % Store component peak channel
        component.amplitude(i) = residuals.avg(peak_channel,peak_sample)*cfg.scaling_factor; % Store component amplitude
        component.sigma(i) = nan; % Store component sigma (spread or timing uncertainty)
        component.shape_param(i) = nan; % Store component shape parameter
        component.fit_success(i) = false; 
    end
     
    component.sign(i) = +1;
    if sign(component.waveform(i,peak_sample))==-1
        component.waveform(i,:) = -component.waveform(i,:);
        component.sign(i) = -1;
    end
    
    % (Visualization for demo mode)
    if strcmp(cfg.demomode,'on')
        plot(residuals.time*1000,component.sign(i)*component.waveform(i,:)*cfg.scaling_factor,'black','Linewidth',2)
        if i<=5
           pause(1) 
        end
    end
    
    
    % Estimate the component projection on the channels
    
    % (Visualization for demo mode)
    if strcmp(cfg.demomode,'on')
        if i<=5
            subplot(5,4,3+(i-1)*4)
        else
            subplot(5,4,3+(5-1)*4)
            cla
        end
        title(['Component ',num2str(i),' regression on channels'])
        axis off
    end
    
    for j=1:size(residuals.avg,1)
        exlcudenans = isnan(residuals.avg(j,:)); % Exclude any NaNs
        component.topo(i,j) = mldivide(component.waveform(i,~exlcudenans)',residuals.avg(j,~exlcudenans)'); % Find the slope with least squares regression
        
        % (Visualization for demo mode)
        if strcmp(cfg.demomode,'on') && i==1
            hold on
            chn_index = ismember(cfg.layout.label,residuals.label{j});
            scatter(cfg.layout.pos(chn_index,1),cfg.layout.pos(chn_index,2),[],component.topo(i,j)'*limits(4),'Filled')
            xlim([-1 1])
            ylim([-.5 .5])
            caxis([limits(3) limits(4)])
            if i<=5
                pause(.005)
            end
        end
        
    end
    
    if strcmp(cfg.demomode,'on')
        cla
        
        chn_index = [];
        for k=1:length(residuals.label)
            chn_index(k) = find(ismember(cfg.layout.label,residuals.label{k}));
        end
        
        interpolation = 1/300;
        [xi,yi] = meshgrid(min(cfg.layout.pos(:,1)):interpolation:max(cfg.layout.pos(:,1)), min(cfg.layout.pos(:,2)):interpolation:max(cfg.layout.pos(:,2)));
        zi = griddata(cfg.layout.pos(chn_index,1),cfg.layout.pos(chn_index,2),component.topo(i,:)'*limits(4),xi,yi);
        topo = pcolor(xi,yi,zi);
        set(topo,'EdgeColor','none')
        xlim([-1 1])
        ylim([-.5 .5])
        caxis([limits(3) limits(4)])
        axis off
        clear('xi','yi','zi','topo');
        title(['Component ',num2str(i),' regression on channels'])
    end
    
    
    % Remove the component estimate from the remaining residuals
    
    component_projection = repmat(component.topo(i,:)',1,size(residuals.avg,2)).*repmat(component.waveform(i,:)',1,size(residuals.avg,1))';
    
    % Check whether removal of the component reduces the residuals
    if nansum(nansum(abs(residuals.avg - component_projection))) < nansum(nansum(abs(residuals.avg))) && sum(component.waveform(i,:))~=0
        residuals.avg = residuals.avg - component_projection;
        
        % (Visualization for demo mode)
        if strcmp(cfg.demomode,'on')
            if i<=5
                pause(1) 
            end
            if i<=5
                subplot(5,4,4+(i-1)*4)
            else
                subplot(5,4,4+(5-1)*4)
                cla
            end
            plot(residuals.time*1000,(component_projection)*cfg.scaling_factor)
            ylim(limits(3:4))
            xlabel('ms')
            ylabel(cfg.scaling_unit)
            title(['Component ',num2str(i),' projection on channels'])
            if i<=5
               pause(1) 
            end
        end
        
    else
        component.waveform(i,:) = [];
        component.latency(i) = [];
        component.channel(i) = [];
        component.amplitude(i) = [];
        component.sigma(i) = [];
        component.shape_param(i) = []; 
        component.topo(i,:) = [];
        component.fit_success(i) = [];
        component.sign(i) = [];
        residual_increase = true;
        
        % (Visualization for demo mode)
        if strcmp(cfg.demomode,'on')
            for j=2:4
                subplot(5,4,j+((5-1)*4))
                cla
                axis off
            end
        end
        
    end
    
    i = i + 1;
    
end

if isempty(component)
    if nansum(nansum(data.avg))==0
        warning('The input data does not contain any signal!')
    end
    error('Modelling failed!')
end

if nansum(nansum(abs(residuals.avg))) == 0 || peak_amplitude <= cfg.amplitude_threshold
    component.rsq(i) = nanmean(1 - nansum(residuals.avg.^2,2)./nansum((data.avg - repmat(nanmean(data.avg,2),1,size(data.avg,2))).^2,2)) - sum(component.rsq);
end

[~,peak_sample] = max(max(abs(residuals.avg),[],1),[],2);
[~,peak_channel] = max(abs(residuals.avg(:,peak_sample)),[],1);

fprintf('Done.\n') 

processing_time_s = toc

fprintf(['\nModel fit success rate is ',num2str(mean(component.fit_success)*100),' percent.\n'])

summed_var = sum(component.rsq(component.fit_success))*100;
var_per_comp = mean(component.rsq(component.fit_success))*100;

fprintf(['\n',num2str(summed_var),' percent variance is explained by '])
if strcmp(cfg.model,'gauss')
    fprintf('gaussian distributions')
elseif strcmp(cfg.model,'gamma')
    fprintf('gamma distributions')
elseif strcmp(cfg.model,'sine')
    fprintf('sine arcs')
end
fprintf('.\n')

fprintf(['\nResidual peak amplitude is ',num2str(nansum(residuals.avg(peak_channel,peak_sample))*cfg.scaling_factor),' ',cfg.scaling_unit,'.\n'])


%% Combine spatially identical components

if cfg.spatial_correlation < 1

    fprintf(['Combining components with spatial correlation > ',num2str(cfg.spatial_correlation),'...'])

    tic;
    correlations = [];
    added_components = [];
    for i=1:size(component.topo,1)-1
       for j=i+1:size(component.topo,1)
           if ~ismember(j,added_components)
                correlation = corr(component.topo(i,:)',component.topo(j,:)');
                if abs(correlation) > cfg.spatial_correlation
                    % Recalibrate weights of first component by weighted mean
                    % (NB: This will add some modeling error)
                    alpha_comp1 = max(abs(component.waveform(i,:))); 
                    alpha_comp2 = max(abs(component.waveform(j,:)));
                    sum_alpha = alpha_comp1 + alpha_comp2;
                    component.topo(i,:) = ((component.topo(i,:)*alpha_comp1) + (sign(correlation)*component.topo(j,:)*alpha_comp2)) / sum_alpha;
                    component.waveform(i,:) = component.waveform(i,:) + sign(correlation)*component.waveform(j,:); % Add spatially identical componentto first component
                    added_components = [added_components j];
                end
           end
       end
    end
    % Remove components that were added to combinations
    component.waveform(added_components,:) = [];
    component.topo(added_components,:) = [];
    component.rsq(added_components) = [];
    component.latency(added_components) = [];
    component.channel(added_components) = [];
    component.amplitude(added_components) = [];
    component.sigma(added_components) = [];
    component.shape_param(added_components) = []; 
    component.fit_success(added_components) = [];
    component.sign(added_components) = [];
    
    fprintf(' Done. \n')
    toc
end


%% Sort the components

if strcmp(cfg.sort_by,'amplitude') && sum(component.fit_success)>0
    fprintf('\nSorting components by peak amplitude.\n')
    [component.amplitude,amplitude_order] = sort(max(abs(component.waveform),[],2).*max(abs(component.topo),[],2),'descend');
    component.waveform = component.waveform(amplitude_order,:);
    component.topo = component.topo(amplitude_order,:);
    component.latency = component.latency(amplitude_order);
    component.channel = component.channel(amplitude_order);
    component.sigma = component.sigma(amplitude_order);
    component.shape_param = component.shape_param(amplitude_order); 
    component.sign = component.sign(amplitude_order);
    component.fit_success = component.fit_success(amplitude_order);
    
elseif strcmp(cfg.sort_by,'variance') && sum(component.fit_success)>0
    fprintf('\nSorting components by explained variance across channels.\n')
    [~,variance_order] = sort(component.rsq,'descend');
    component.waveform = component.waveform(variance_order,:);
    component.topo = component.topo(variance_order,:);
    component.rsq = component.rsq(variance_order);
    component.amplitude = max(abs(component.waveform),[],2).*max(abs(component.topo),[],2);
    component.latency = component.latency(variance_order);
    component.channel = component.channel(variance_order);
    component.sigma = component.sigma(variance_order);
    component.shape_param = component.shape_param(variance_order); 
    component.sign = component.sign(variance_order);
    component.fit_success = component.fit_success(variance_order);
    
elseif strcmp(cfg.sort_by,'latency') && sum(component.fit_success)>0
    fprintf('\nSorting components by latency.\n')
    [~,latency_order] = sort(component.latency,'ascend');
    component.waveform = component.waveform(latency_order,:);
    component.topo = component.topo(latency_order,:);
    component.rsq = component.rsq(latency_order);
    component.amplitude = max(abs(component.waveform),[],2).*max(abs(component.topo),[],2);
    component.latency = component.latency(latency_order);
    component.channel = component.channel(latency_order);
    component.sigma = component.sigma(latency_order);
    component.shape_param = component.shape_param(latency_order); 
    component.sign = component.sign(latency_order);
    component.fit_success = component.fit_success(latency_order);
    
else
    component.amplitude = max(abs(component.waveform),[],2).*max(abs(component.topo),[],2);
    
end


%% Store the components in FieldTrip data format

comp = residuals; % Copy residuals
comp = rmfield(comp,'avg');
comp.time = {};
comp.time{1,1} = residuals.time;
comp.trial{1,1} = component.waveform;
comp.topolabel = comp.label;
comp.label = {};
for i=1:size(component.waveform,1)
   cat_levels = 1; % Find the necessary categorization levels (1 = latency, 2 = + channel, 3 = + amplitude, 4 = + letters)
   if sum(ismember(component.latency,component.latency(i)))>1
       cat_levels = 2;
       if sum(ismember(component.channel(ismember(component.latency,component.latency(i))),component.channel(i)))>1
           cat_levels = 3;
       end
   end
   if cat_levels==1
       comp.label{i,1} = [sprintf(num2str(component.latency(i))),'(/',component.channel{i},'/',sprintf('%.1f',component.sign(i)*component.amplitude(i)*cfg.scaling_factor),')'];
   elseif cat_levels==2
       comp.label{i,1} = [sprintf(num2str(component.latency(i))),'/',component.channel{i},'(/',sprintf('%.1f',component.sign(i)*component.amplitude(i)*cfg.scaling_factor),')'];
   elseif cat_levels==3
       comp.label{i,1} = [sprintf(num2str(component.latency(i))),'/',component.channel{i},'/',sprintf('%.1f',component.sign(i)*component.amplitude(i)*cfg.scaling_factor)];
   end
   
end
% Check whether numeric suffix should be added after otherwise identical component labels
for i=1:size(component.waveform,1)
    if sum(ismember(comp.label,comp.label{i,1}))>1
        duplicates = find(ismember(comp.label,comp.label{i,1}));
        for j=1:length(duplicates)
            comp.label{j,1} = [comp.label{j,1},'(/',num2str(j),')'];
        end
    end
end
comp.topo = component.topo';
comp.topo(isnan(comp.topo)) = 0;
comp.topo(abs(comp.topo)==inf) = 0;
comp.unmixing = pinv(comp.topo);
comp.latency = component.latency';
comp.channel = component.channel';
comp.amplitude = component.amplitude;
comp.sign = component.sign';
comp.sigma = component.sigma';
comp.shape = component.shape_param';
% Calculates the differential entropy on ms spike bins. 
% NB: Assumes that spikes are distributed across at least two milliseconds
% and appropriate sampling rate is applied. 
% Measures the spike timing uncertainty, as average spike timing information in bits per ms.
if strcmp(cfg.model,'gauss')
    component.entropy = 1/2*log(2*pi*exp(1)*((component.sigma*1000).^2)); 
    comp.entropy = component.entropy';
end
comp.rsq = component.rsq';
comp.fit_success = component.fit_success';


%% Show the components

if strcmp(cfg.show_components,'yes')
    db_cfg = []; db_cfg.layout = cfg.layout; db_cfg.channel = 1:10; db_cfg.viewmode = 'component'; db_cfg.zlim = 'maxabs';
    evalc('ft_databrowser(db_cfg,comp)');
    figure 
    plot(data.time*1000, component.waveform'.*repmat(max(abs(component.topo),[],2),1,size(component.waveform,2))'*cfg.scaling_factor,'linewidth',3)
    set(gca,'fontsize',13)
    if strcmp(cfg.show_residuals,'yes')
        hold on
        plot(data.time*1000, max(abs(residuals.avg),[],1)*cfg.scaling_factor,'color','black','LineStyle',':','linewidth',3)
        legend([comp.label;'Residuals'],'Location','NortheastOutside','fontsize',10)
    else
        legend(comp.label,'Location','NortheastOutside','fontsize',10)
    end
    xlabel('ms','fontsize',13)
    ylabel(cfg.scaling_unit,'fontsize',13)
    title('Max. proj. on channels','interpreter','none','fontsize',13)
end


end