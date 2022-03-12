function Entropy = IT_Single_Neuron_Entropy(spikes, Fs, binSize, nStates, varargin)
% Entropy = IT_Single_Neuron_Entropy(spikes, Fs, timeBinSize, stateBinEdges...)
%
% Calculates the entropy of a single neuron. Spike counts are binned w.r.t
% time, and the probability of the neuron in one of nStates is calculated
% relative to the maximum firing rate of the neuron.
%
% A warning is given if one of your states has a count of less than 10. If
% one state has a count of 0, a NaN entropy result will be produced.
%
%
% Inputs            spikes              index of spike times (in s)
%
%                   Fs                  sampling period of spikes (e.g.,
%                                       .001 for 1000-Hz sampling rate)
%
%                   binSize             size of time-based bin to
%                                       discretize data (in s)
%
%                   nStates             number of possible states (integer
%                                       value, states are determined by
%                                       equal-width binning)
%
%                   Optional            'Name', Value
%
%                   PlotOutput          creates a plot of the absolute and
%                                       relative probability distributions
%                                       used to calculate entropy
%                                       (true/false, default is false)
%
%
% Output            Entropy             structure containing user-selected
%                                       parameters, probability
%                                       distributions, and single-neuron
%                                       entropy (in bits)
%
% Contributed by Ryan Glanz (ryan-glanz@uiowa.edu)
% Last updated 9.6.2020 by RG
%

%% Parameters
params = inputParser;
params.addRequired('spikes', @(x) isnumeric(x));
params.addRequired('Fs', @(x) isnumeric(x));
params.addRequired('binSize', @(x) isnumeric(x));
params.addRequired('nStates', @(x) isnumeric(x));
params.addParameter('PlotOutput', false, @islogical);
params.parse(spikes, Fs, binSize, nStates, varargin{:});

plotOutput = params.Results.PlotOutput;

%% Bin data w.r.t. time
logSpikes = ismember(0:Fs:spikes(end), spikes); % Logical array of spike times
timeToTruncate = mod(length(logSpikes), binSize * (1/Fs)) * Fs;
warningMessage  = strcat("Truncating ", num2str(timeToTruncate),...
    " s from data.");
warning(warningMessage)

logSpikes = logSpikes(1:end-(timeToTruncate * (1/Fs))); % Truncate logical vector for binning
binSpikes = reshape(logSpikes, [], binSize * (1/Fs));
binSpikes = sum(binSpikes, 2);  % Sum of time bins

%% State determination
maxBinSpikes = max(binSpikes);
nBinSpikes = binSpikes / maxBinSpikes;    % Normalize firing rate values
stateBinEdges = linspace(0, 1, nStates + 1);    % State bins
stateCount = histcounts(nBinSpikes, stateBinEdges);
if any(stateCount <= 10)
    warning("At least one state bin has a count of <= 10")
end

%% Assign variables and calculate entropy
Entropy.Params.Fs = Fs;
Entropy.Params.binSize = binSize;
Entropy.Params.nStates = nStates;
Entropy.Params.timeTruncated = timeToTruncate;
Entropy.MaxFiringRate = maxBinSpikes / binSize;
Entropy.StateBinEdges = maxBinSpikes / binSize * stateBinEdges;
Entropy.ProbDistribution.Absolute = stateCount;
Entropy.ProbDistribution.Relative = stateCount / sum(stateCount);

Entropy.H = sum(arrayfun(@(x) Entropy.ProbDistribution.Relative(x)*...
    log2(1/Entropy.ProbDistribution.Relative(x)),...
    1:length(Entropy.ProbDistribution.Relative)));

%% Plot
if plotOutput
    figure
    subplot(1, 2, 1)
    bar(Entropy.ProbDistribution.Absolute)
    title('Absolute Probability Distribution')
    ylabel('Counts')
    xlabel('State (sps, right-edge)')
    xticks(1:length(Entropy.StateBinEdges) - 1)
    xticklabels(Entropy.StateBinEdges(2:end))
    
    subplot(1, 2, 2)
    bar(Entropy.ProbDistribution.Relative)
    title('Relative Probability Distribution')
    ylabel('Percent')
    xlabel('State (percentile, right-edge)')
    xticks(1:length(stateBinEdges) - 1)
    xticklabels(stateBinEdges(2:end))
end

end
