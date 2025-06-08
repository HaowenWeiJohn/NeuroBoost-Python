function visualize_eeg_evoked(EEG, varargin)
    % VISUALIZE_EEG_EVOKED - Create publication-quality visualization of TMS-evoked potentials
    %
    % This function generates a comprehensive plot of EEG evoked responses following
    % TMS stimulation. It applies baseline correction, averages across trials, and
    % creates a formatted visualization suitable for analysis and publication.
    %
    % SYNTAX:
    %   visualize_eeg_evoked(EEG)
    %   visualize_eeg_evoked(EEG, 'Parameter', Value, ...)
    %
    % REQUIRED INPUT:
    %   EEG              - EEGLAB EEG structure containing epoched data
    %                      Must have fields: .data, .times, .chanlocs
    %                      Data format: [channels × timepoints × trials]
    %
    % OPTIONAL PARAMETERS (Name-Value Pairs):
    %   'BaselineWindow'   - [start_ms, end_ms] Time window for baseline correction
    %                        Default: [EEG.times(1), 0] (epoch start to stimulus onset)
    %                        Example: [-1000, -2] removes mean from 1000ms to 2ms pre-stimulus
    %                        Use [] to skip baseline correction if already applied
    %   'TimeLimits'       - [start_ms, end_ms] X-axis display range for the plot
    %                        Default: [EEG.times(1), EEG.times(end)] (entire epoch)
    %                        Example: [-100, 400] shows 100ms pre to 400ms post-stimulus
    %   'AmplitudeLimits'  - [min_uV, max_uV] Y-axis amplitude range for the plot
    %                        Default: [-50, 50] (±50 µV range)
    %                        Example: [-20, 30] sets amplitude scale from -20 to +30 µV
    %   'PlotTitle'        - String for plot title
    %                        Default: 'TMS-Evoked Potentials'
    %   'ChannelNames'     - Cell array of channel names to plot, or single string
    %                        Default: {} (plot all channels)
    %                        Example: {'Cz', 'FCz', 'Pz'} plots only these channels
    %   'SaveFig'          - Boolean indicating whether to save the figure
    %                        Default: false
    %   'OutputPath'       - Path to save the figure (required if SaveFig is true)
    %                        Default: pwd (current directory)
    %   'Filename'         - Base name for saved figure files (required if SaveFig is true)
    %                        Default: 'eeg_evoked_plot'
    %   'TitleInterpreter' - String specifying title interpreter behavior
    %                        Default: '' (use MATLAB default)
    %                        Options: 'none', 'tex', 'latex'
    %
    % OUTPUTS:
    %   Creates a figure window with the evoked potential visualization
    %   Selected channels are plotted as separate lines on the same axes
    %
    % EXAMPLES:
    %   % Minimal usage (all defaults)
    %   visualize_eeg_evoked(EEG);
    %
    %   % Specify only what you need
    %   visualize_eeg_evoked(EEG, 'TimeLimits', [-100, 400], 'AmplitudeLimits', [-20, 30]);
    %
    %   % Plot specific channels with custom title
    %   visualize_eeg_evoked(EEG, 'ChannelNames', {'Cz', 'FCz', 'Pz'}, ...
    %                       'PlotTitle', 'Midline Channels', 'SaveFig', true);
    %
    %   % Full customization
    %   visualize_eeg_evoked(EEG, 'BaselineWindow', [-1000, -2], ...
    %                       'TimeLimits', [-100, 400], 'AmplitudeLimits', [-30, 40], ...
    %                       'PlotTitle', 'TMS-evoked potentials - Condition_A', ...
    %                       'ChannelNames', {'C5', 'C3', 'C1', 'Cz'}, ...
    %                       'SaveFig', true, 'OutputPath', '/path/to/output', ...
    %                       'Filename', 'tep_analysis', 'TitleInterpreter', 'none');
    %
    %   % Backward compatibility (old positional syntax still works)
    %   visualize_eeg_evoked(EEG, [-1000, -2], [-100, 400], [-20, 30], 'My Title');
    %
    % NOTES:
    %   - Function assumes data is already epoched around stimulus events
    %   - Baseline correction is applied to a copy, original data unchanged
    %   - Channel names are matched against EEG.chanlocs.labels (case-insensitive)
    %   - Time zero should correspond to TMS pulse onset in the epoch structure
    %   - Use TitleInterpreter 'none' for titles with underscores to prevent subscript formatting
    %   - Vertical dashed line marks TMS onset (time = 0)
    %
    % See also: POP_RMBASE, MEAN, PLOT, INPUTPARSER

    % Add robust input validation first
    if ~isstruct(EEG) || ~isfield(EEG, 'data') || ~isfield(EEG, 'times') || ~isfield(EEG, 'chanlocs')
        error('EEG must be a valid EEGLAB structure with .data, .times, and .chanlocs fields');
    end

    % Check for empty or invalid data
    if isempty(EEG.data) || size(EEG.data, 1) == 0 || size(EEG.data, 2) == 0
        error('EEG.data cannot be empty');
    end

    % Check for valid time vector
    if isempty(EEG.times) || length(EEG.times) ~= size(EEG.data, 2)
        error('EEG.times must match the number of time points in EEG.data');
    end

    % Check for valid channel locations
    if length(EEG.chanlocs) ~= size(EEG.data, 1)
        error('Number of channel locations must match number of data channels');
    end

    % More nuanced check for epoched vs continuous data
    if size(EEG.data, 3) == 1 && length(EEG.times) > 1000
        warning('EEG data appears to be continuous (single trial with many timepoints). Expected epoched data for averaging.');
    elseif size(EEG.data, 3) == 1
        warning('EEG data contains only one trial. No averaging will be performed.');
    end

    % Check if using old positional syntax for backward compatibility
    if nargin >= 2 && ~ischar(varargin{1}) && ~isstring(varargin{1})
        % Old positional syntax detected
        fprintf('Using legacy positional argument syntax. Consider switching to name-value pairs for better flexibility.\n');
        
        % Map old positional arguments to new parameter names
        baselineWindow = [];
        timeLimits = [];
        amplitudeLimits = [];
        plotTitle = 'TMS-Evoked Potentials';
        channel_names = {};
        save_fig = false;
        output_path = pwd;
        filename = 'eeg_evoked_plot';
        titleInterpreter = '';
        
        if nargin >= 2, baselineWindow = varargin{1}; end
        if nargin >= 3, timeLimits = varargin{2}; end
        if nargin >= 4, amplitudeLimits = varargin{3}; end
        if nargin >= 5, plotTitle = varargin{4}; end
        
        % Handle the complex backward compatibility logic for remaining arguments
        if nargin >= 6 && (islogical(varargin{5}) || (isnumeric(varargin{5}) && (varargin{5} == 0 || varargin{5} == 1)))
            % Old syntax: 6th argument is save_fig (boolean)
            if nargin >= 6, save_fig = varargin{5}; end
            if nargin >= 7, output_path = varargin{6}; end
            if nargin >= 8, filename = varargin{7}; end
            if nargin >= 9, titleInterpreter = varargin{8}; end
        else
            % New syntax: 6th argument is channel_names
            if nargin >= 6, channel_names = varargin{5}; end
            if nargin >= 7, save_fig = varargin{6}; end
            if nargin >= 8, output_path = varargin{7}; end
            if nargin >= 9, filename = varargin{8}; end
            if nargin >= 10, titleInterpreter = varargin{9}; end
        end
        
    else
        % New name-value pair syntax
        p = inputParser;
        
        % Define default values
        % Smart default for baseline window: use epoch start to closest-to-zero time
        if any(EEG.times <= 0)
            % If epoch includes pre-stimulus time, use epoch start to time closest to 0
            baseline_end = EEG.times(find(EEG.times <= 0, 1, 'last'));
        else
            % If epoch starts after stimulus, use first 10% of epoch for baseline
            baseline_end = EEG.times(1) + 0.1 * (EEG.times(end) - EEG.times(1));
        end
        defaultBaselineWindow = [EEG.times(1), baseline_end];
        
        defaultTimeLimits = [EEG.times(1), EEG.times(end)];
        defaultAmplitudeLimits = [-50, 50];
        defaultPlotTitle = 'TMS-Evoked Potentials';
        defaultChannelNames = {};
        defaultSaveFig = false;
        defaultOutputPath = pwd;
        defaultFilename = 'eeg_evoked_plot';
        defaultTitleInterpreter = '';
        
        % Add parameters to parser with improved validation
        addParameter(p, 'BaselineWindow', defaultBaselineWindow, @(x) isnumeric(x) && (isempty(x) || length(x) == 2));
        addParameter(p, 'TimeLimits', defaultTimeLimits, @(x) isnumeric(x) && (isempty(x) || (length(x) == 2 && x(1) < x(2))));
        addParameter(p, 'AmplitudeLimits', defaultAmplitudeLimits, @(x) isnumeric(x) && (isempty(x) || (length(x) == 2 && x(1) < x(2))));
        addParameter(p, 'PlotTitle', defaultPlotTitle, @(x) ischar(x) || isstring(x));
        addParameter(p, 'ChannelNames', defaultChannelNames, @(x) iscell(x) || ischar(x) || isstring(x) || isempty(x));
        addParameter(p, 'SaveFig', defaultSaveFig, @(x) islogical(x) || isnumeric(x));
        addParameter(p, 'OutputPath', defaultOutputPath, @(x) ischar(x) || isstring(x));
        addParameter(p, 'Filename', defaultFilename, @(x) ischar(x) || isstring(x));
        addParameter(p, 'TitleInterpreter', defaultTitleInterpreter, @(x) ischar(x) || isstring(x));
        
        % Parse inputs
        parse(p, varargin{:});
        
        % Extract parsed values
        baselineWindow = p.Results.BaselineWindow;
        timeLimits = p.Results.TimeLimits;
        amplitudeLimits = p.Results.AmplitudeLimits;
        plotTitle = p.Results.PlotTitle;
        channel_names = p.Results.ChannelNames;
        save_fig = p.Results.SaveFig;
        output_path = p.Results.OutputPath;
        filename = p.Results.Filename;
        titleInterpreter = p.Results.TitleInterpreter;
        
        % Additional validation after parsing
        if ~isempty(baselineWindow)
            if baselineWindow(1) < EEG.times(1) || baselineWindow(2) > EEG.times(end)
                warning('BaselineWindow [%.1f, %.1f] extends beyond data range [%.1f, %.1f]. Clipping to data range.', ...
                    baselineWindow(1), baselineWindow(2), EEG.times(1), EEG.times(end));
                baselineWindow(1) = max(baselineWindow(1), EEG.times(1));
                baselineWindow(2) = min(baselineWindow(2), EEG.times(end));
            end
        end
        
        if ~isempty(timeLimits)
            if timeLimits(1) < EEG.times(1) || timeLimits(2) > EEG.times(end)
                warning('TimeLimits [%.1f, %.1f] extends beyond data range [%.1f, %.1f]. Clipping to data range.', ...
                    timeLimits(1), timeLimits(2), EEG.times(1), EEG.times(end));
                timeLimits(1) = max(timeLimits(1), EEG.times(1));
                timeLimits(2) = min(timeLimits(2), EEG.times(end));
            end
        else
            timeLimits = [EEG.times(1), EEG.times(end)];
        end
        
        if isempty(amplitudeLimits)
            % Auto-calculate amplitude limits based on data
            temp_data = mean(EEG.data, 3);
            data_range = [min(temp_data(:)), max(temp_data(:))];
            margin = 0.1 * (data_range(2) - data_range(1));
            amplitudeLimits = [data_range(1) - margin, data_range(2) + margin];
        end
    end

    % Convert logical save_fig to boolean if needed
    save_fig = logical(save_fig);

    % Validate and process channel_names parameter
    if isempty(channel_names)
        channel_indices = 1:EEG.nbchan;  % Use all channels
        fprintf('Plotting all %d channels\n', EEG.nbchan);
    else
        % Validate channel_names input
        if ischar(channel_names) || isstring(channel_names)
            channel_names = {char(channel_names)};  % Convert single string to cell array
        elseif ~iscell(channel_names)
            error('ChannelNames must be a cell array of strings or a single string');
        end
        
        % Find indices of requested channels (more efficient)
        channel_indices = [];
        missing_channels = {};
        
        % Pre-compute lowercase labels for faster comparison
        eeg_labels_lower = lower({EEG.chanlocs.labels});
        
        for i = 1:length(channel_names)
            idx = find(strcmpi(eeg_labels_lower, lower(channel_names{i})), 1);
            if ~isempty(idx)
                channel_indices(end+1) = idx;  % More efficient than concatenation
            else
                missing_channels{end+1} = channel_names{i};
            end
        end
        
        % Report missing channels
        if ~isempty(missing_channels)
            warning('The following channels were not found: %s', strjoin(missing_channels, ', '));
        end
        
        % Check if any valid channels were found
        if isempty(channel_indices)
            error('None of the specified channels were found in the EEG data');
        end
        
        fprintf('Plotting %d of %d requested channels\n', length(channel_indices), length(channel_names));
    end

    % Validate save_fig parameters
    if save_fig
        if ~ischar(output_path) && ~isstring(output_path)
            error('OutputPath must be a valid string or char array when SaveFig is true');
        end
        if ~ischar(filename) && ~isstring(filename)
            error('Filename must be a valid string or char array when SaveFig is true');
        end
        % Create output directory if it doesn't exist
        if ~exist(output_path, 'dir')
            mkdir(output_path);
        end
    end

    % Apply baseline correction to a working copy (preserves original data)
    EEG_epoch = EEG;
    if ~isempty(baselineWindow)
        EEG_epoch = pop_rmbase(EEG_epoch, baselineWindow);
    end
    
    % Calculate grand average evoked potential across all trials
    % Result: [channels × timepoints] matrix of averaged responses
    evoked_data = mean(EEG_epoch.data, 3);
    
    % Select only the requested channels
    evoked_data_selected = evoked_data(channel_indices, :);
    
    % Create the evoked potential visualization
    figure;
    plot(EEG_epoch.times, evoked_data_selected');  % Transpose for proper orientation
    hold on;
    
    % Set axis limits and labels
    xlim(timeLimits);
    ylim(amplitudeLimits);
    
    % Set title with optional interpreter specification
    if isempty(titleInterpreter)
        % Use MATLAB default behavior (don't specify interpreter)
        title(plotTitle, 'FontSize', 11, 'FontWeight', 'bold');
    else
        % Use specified interpreter
        title(plotTitle, 'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', titleInterpreter);
    end
    
    xlabel('Time (ms)', 'FontSize', 12);
    ylabel('Amplitude (µV)', 'FontSize', 12);
    
    % Add legend if specific channels were selected (and not too many)
    if ~isempty(channel_names) && length(channel_indices) <= 10
        legend({EEG.chanlocs(channel_indices).labels}, 'Location', 'best', 'FontSize', 8);
    end
    
    % Add vertical line at stimulus onset (time = 0)
    if timeLimits(1) <= 0 && timeLimits(2) >= 0
        plot([0 0], amplitudeLimits, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
    
    % Enhance plot appearance
    grid on;
    grid minor;
    % Set tick label font size only (not affecting title, xlabel, ylabel)
    set(gca, 'FontSize', 10, 'TitleFontSizeMultiplier', 0.8);  % Preserve smaller title size
    
    hold off;

    % Add optional figure saving with error handling
    if save_fig
        try
            % Remove extension from filename if present to avoid duplication
            [~, base_filename, ~] = fileparts(filename);
            
            exportgraphics(gcf, fullfile(output_path, [base_filename '.png']), 'Resolution', 300);
            saveas(gcf, fullfile(output_path, [base_filename '.fig']));
            fprintf('Figure saved successfully to: %s\n', output_path);
        catch ME
            warning('MATLAB:SaveFigureFailed', 'Failed to save figure: %s', ME.message);
        end
    end
end