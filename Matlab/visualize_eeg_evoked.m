function visualize_eeg_evoked(EEG, baselineWindow, timeLimits, amplitudeLimits, plotTitle, save_fig, output_path, filename, titleInterpreter)
    % VISUALIZE_EEG_EVOKED - Create publication-quality visualization of TMS-evoked potentials
    %
    % This function generates a comprehensive plot of EEG evoked responses following
    % TMS stimulation. It applies baseline correction, averages across trials, and
    % creates a formatted visualization suitable for analysis and publication.
    %
    % SYNTAX:
    %   visualize_eeg_evoked(EEG, baselineWindow, timeLimits, amplitudeLimits, plotTitle)
    %   visualize_eeg_evoked(EEG, baselineWindow, timeLimits, amplitudeLimits, plotTitle, save_fig, output_path, filename)
    %   visualize_eeg_evoked(EEG, baselineWindow, timeLimits, amplitudeLimits, plotTitle, save_fig, output_path, filename, titleInterpreter)
    %
    % INPUTS:
    %   EEG              - EEGLAB EEG structure containing epoched data
    %                      Must have fields: .data, .times, .chanlocs
    %                      Data format: [channels × timepoints × trials]
    %   baselineWindow   - [start_ms, end_ms] Time window for baseline correction
    %                      Example: [-1000, -2] removes mean from 1000ms to 2ms pre-stimulus
    %                      Use [] to skip baseline correction if already applied
    %   timeLimits       - [start_ms, end_ms] X-axis display range for the plot
    %                      Example: [-100, 400] shows 100ms pre to 400ms post-stimulus
    %                      Use [] to display the entire epoch duration
    %   amplitudeLimits  - [min_uV, max_uV] Y-axis amplitude range for the plot
    %                      Example: [-20, 30] sets amplitude scale from -20 to +30 µV
    %                      Choose based on your expected TEP amplitude range
    %   plotTitle        - String for plot title (default: 'TMS-Evoked Potentials')
    %                      Should describe the processing stage or condition
    %   save_fig         - Boolean indicating whether to save the figure (default: false)
    %   output_path      - Path to save the figure (required if save_fig is true)
    %   filename         - Base name for saved figure files (required if save_fig is true)
    %   titleInterpreter - String specifying title interpreter behavior (default: not specified)
    %                      Options: 'none' (literal text, no interpretation), 
    %                               'tex' (TeX markup), 'latex' (LaTeX markup)
    %                      If not specified, uses MATLAB's default behavior
    %
    % OUTPUTS:
    %   Creates a figure window with the evoked potential visualization
    %   All channels are plotted as separate lines on the same axes
    %
    % EXAMPLE:
    %   % Standard TEP visualization after TESA processing
    %   visualize_eeg_evoked(EEG, [-1000 -2], [-100 400], [-20 30], ...
    %                       'Final TESA-processed TMS-evoked potentials');
    %
    %   % With figure saving and title containing underscores
    %   visualize_eeg_evoked(EEG, [-1000 -2], [-100 400], [-20 30], ...
    %                       'TMS-evoked potentials - Center_110_part1', ...
    %                       true, '/path/to/output', 'tep_figure', 'none');
    %
    % NOTES:
    %   - Function assumes data is already epoched around stimulus events
    %   - Baseline correction is applied to a copy, original data unchanged
    %   - All channels are plotted; consider channel selection for cleaner visualization
    %   - Time zero should correspond to TMS pulse onset in the epoch structure
    %   - Use titleInterpreter 'none' for titles with underscores to prevent subscript formatting
    %
    % See also: POP_RMBASE, MEAN, PLOT

    % Add robust input validation first
    if ~isstruct(EEG) || ~isfield(EEG, 'data') || ~isfield(EEG, 'times')
        error('EEG must be a valid EEGLAB structure with .data and .times fields');
    end

    if size(EEG.data, 3) < 2
        warning('EEG data appears to be continuous. Expected epoched data for averaging.');
    end

    % Input validation and default parameter handling
    if nargin < 9 || isempty(titleInterpreter)
        titleInterpreter = '';  % Empty means don't specify interpreter (use MATLAB default)
    end
    if nargin < 8 || isempty(filename)
        filename = 'eeg_evoked_plot';
    end
    if nargin < 7 || isempty(output_path)
        output_path = pwd; % Default to current directory
    end
    if nargin < 6 || isempty(save_fig)
        save_fig = false;
    end
    if nargin < 5 || isempty(plotTitle)
        plotTitle = 'TMS-Evoked Potentials';
    end
    if nargin < 4 || isempty(amplitudeLimits)
        amplitudeLimits = [-50 50];  % Default ±50 µV range
    end
    if nargin < 3 || isempty(timeLimits)
        timeLimits = [EEG.times(1), EEG.times(end)]; % Use the entire epoch duration
    end
    if nargin < 2 || isempty(baselineWindow)
        baselineWindow = [EEG.times(1), 0]; % Default: epoch start to stimulus onset
    end

    % Validate save_fig parameters
    if save_fig
        if ~ischar(output_path) && ~isstring(output_path)
            error('output_path must be a valid string or char array when save_fig is true');
        end
        if ~ischar(filename) && ~isstring(filename)
            error('filename must be a valid string or char array when save_fig is true');
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
    
    % Create the evoked potential visualization
    figure;
    plot(EEG_epoch.times, evoked_data');  % Transpose for proper orientation
    hold on;
    
    % Set axis limits and labels
    xlim(timeLimits);
    ylim(amplitudeLimits);
    
    % Set title with optional interpreter specification
    if isempty(titleInterpreter)
        % Use MATLAB default behavior (don't specify interpreter)
        title(plotTitle, 'FontSize', 14, 'FontWeight', 'bold');
    else
        % Use specified interpreter
        title(plotTitle, 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', titleInterpreter);
    end
    
    xlabel('Time (ms)', 'FontSize', 12);
    ylabel('Amplitude (µV)', 'FontSize', 12);
    
    % Add vertical line at stimulus onset (time = 0)
    if timeLimits(1) <= 0 && timeLimits(2) >= 0
        plot([0 0], amplitudeLimits, 'k--', 'LineWidth', 1.5);
    end
    
    % Enhance plot appearance
    grid on;
    grid minor;
    set(gca, 'FontSize', 10);
    
    hold off;

    % Add optional figure saving with error handling
    if save_fig
        try
            exportgraphics(gcf, fullfile(output_path, [filename '.png']), 'Resolution', 300);
            saveas(gcf, fullfile(output_path, [filename '.fig']));
            fprintf('Figure saved successfully to: %s\n', output_path);
        catch ME
            warning('MATLAB:SaveFigureFailed', 'Failed to save figure: %s', ME.message);
        end
    end
end