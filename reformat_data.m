%% Script to put Rat Scavenge input data into correct Format for running through model
% Data taken from the following link:  https://www.dropbox.com/s/94dhsgnx2cfs3jx/hc_data_raw.mat?dl=0
% Script following Kording Lab: https://github.com/KordingLab/Neural_Decoding/blob/master/Examples_hippocampus/Example_format_data_hc.ipynb
% Dataset described as: "Multi-unit recordings from the rat hippocampus made during open field foraging" 

% Importing the raw data
rawData = load("hc_data_raw.mat");

% Parsing input data and defining variables
posStartTimes = rawData.pos_times(1); % Start times at which the positions were recorded 
posEndTimes = rawData.pos_times(end);

% Bin data
timeBins = .2;  % time in sec of time bins
formattedSpikeData = binSpikes(rawData.spike_times, posStartTimes, posEndTimes, timeBins);
formattedPosData = binData(rawData.pos, rawData.pos_times, timeBins, posStartTimes, posEndTimes);

% Output data
outDir = pwd;
outFile = fullfile(outDir, 'reformatted_scavenge_data.mat');
save(outFile, 'formattedPosData', 'formattedSpikeData');


%% Functions

function binnedSpikes = binSpikes(spikeTimes, startWin, endWin, binNum)
    % Put the spikes for the neurons into a specified number of bins

    % spikeTimes: matrix containing the spike times for a given neuron
    % startWin  : start time for spike binning
    % endWin    : end time for spike binning
    % binNum    : bin size 
    
    binEdge = startWin:binNum:endWin;  % Find the start and end points of each bins
    binCount = length(binEdge) - 1;    % Count the number of bins
    neuronNum = length(spikeTimes);    % Count total number of neurons
   
    binnedSpikes = zeros(binCount, neuronNum);
    
    % Create an array containing number of spikes in each bin for a given neuron
    for neuron = 1:neuronNum
        counts = histcounts(spikeTimes{neuron}, binEdge);
        binnedSpikes(:, neuron) = counts';
    end

end

function binnedOutput = binData(positionData, times, binWidth, startWin, endWin)
    binEdge = startWin:binWidth:endWin;  % Find the start and end points of each bins
    binCount = length(binEdge) - 1;    % Count the number of bins
    
    features = size(positionData, 2) % Count dimensions of data
    binnedOutput = zeros(binCount, features);

    for bin = 1:binCount
        i = times >=binEdge(bin) & times < binEdge(bin + 1);
        for j = 1:features
            binnedOutput(bin, j) = mean(positionData(i, j), 'omitnan');
        end
    end
end   
