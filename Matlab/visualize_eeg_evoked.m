function visualize_eeg_evoked(EEG, baselineWindow, timeLimits, amplitudeLimits, plotTitle)
    % visualizeEEGEvoked - Visualizes the EEG evoked potentials after basic preprocessing.
    %
    % Syntax:
    %   visualizeEEGEvoked(EEG, baselineWindow, timeLimits, amplitudeLimits, plotTitle)
    %
    % Inputs:
    %   EEG              - EEG structure after preprocessing
    %   baselineWindow   - Two-element vector [start, end] for baseline removal (e.g., [-750 0])
    %   timeLimits       - Two-element vector [start, end] for x-axis limits (e.g., [-100 400])
    %   amplitudeLimits  - Two-element vector [min, max] for y-axis limits (e.g., [-50 50])
    %   plotTitle        - Title for the plot (e.g., 'After Basic Preprocessing')
    %
    % Example:
    %   visualizeEEGEvoked(EEG, [-750 0], [-100 400], [-50 50], 'Evoked Potential Visualization');

    if nargin < 5
        plotTitle = 'Evoked';
    end
    if nargin < 4
        amplitudeLimits = [-50 50];
    end
    if nargin < 3 || isempty(timeLimits)
        timeLimits = [EEG.times(1), EEG.times(end)]; % Use the entire epoch
    end
    if nargin < 2 || isempty(baselineWindow)
        baselineWindow = [EEG.times(1), 0]; % Use from the first sample to 0 ms as baseline
    end

    % Baseline correction
    EEG_epoch = EEG;
    EEG_epoch = pop_rmbase(EEG_epoch, baselineWindow);
    
    % Calculate the average evoked potential across trials
    evoked_dummy = mean(EEG_epoch.data, 3);
    
    % Plotting
    figure;
    plot(EEG_epoch.times, evoked_dummy(:,:));
    hold on;
    xlim(timeLimits);
    ylim(amplitudeLimits);
    title(plotTitle);
    xlabel('Time (ms)');
    ylabel('Amplitude (µV)');
    hold off;
end