% Conservative optimized version of StOptimization.m.
%
% Goals:
% 1) Keep the original THD calculation logic based on MATLAB thd().
% 2) Sweep multiple p values and DAC quantization steps R.
% 3) Compare quantized S(t) with quantized Se(t).
% 4) Automatically record the best R values and export result tables.
% 5) Use clearer variable names and comments.
%
% Note:
% Set thdSamplingMode = "original" to preserve the original code behavior
% using fs = 99999. Set it to "effective" to test whether using the actual
% expanded waveform sampling rate changes the result.

clc;
clear;
close all;

%% User settings
pList = [3 7 15 31];       % For S(t), OSR = 4*(p+1)
RList = 4:128;             % DAC quantization steps

targetSignalFrequency = 5e5;
numCycles = 20;
holdSamples = 100;         % Original code repeats each DAC point 100 times
numHarmonics = 100;

LPFBandwidthRatio = 2;     % LPF cutoff frequency = 2 * signal frequency
LPFOrder = 2;

targetTHDPercent = 0.10;
targetErrorPercent = 0.05;

thdSamplingMode = "original";   % "original" or "effective"
originalTHDSamplingRate = 99999;

%% Main sweep
allResults = table();
summary = table();

figure(1);
plotLayout = tiledlayout(2, 1);
title(plotLayout, 'Conservative S(t) and Se(t) Parameter Sweep');

for pS = pList
    OSR = 4*(pS + 1);
    pSe = OSR/4;                  % Se(t) uses p = OSR/4
    dacClock = targetSignalFrequency * OSR;
    effectiveSamplingRate = dacClock * holdSamples;

    if thdSamplingMode == "effective"
        thdSamplingRate = effectiveSamplingRate;
    else
        thdSamplingRate = originalTHDSamplingRate;
    end

    % Preserve the original endpoint-included time vector.
    time = 0 : 1/dacClock : numCycles/targetSignalFrequency;

    waveformS = generateSWaveform(time, targetSignalFrequency, pS, OSR);
    waveformSe = generateSeWaveform(time, targetSignalFrequency, pSe, OSR);

    [metricsS, resultsS] = evaluateWaveformAcrossR( ...
        waveformS, RList, dacClock, holdSamples, targetSignalFrequency, ...
        LPFBandwidthRatio, LPFOrder, thdSamplingRate, numHarmonics, ...
        "S", pS, OSR, pSe);

    [metricsSe, resultsSe] = evaluateWaveformAcrossR( ...
        waveformSe, RList, dacClock, holdSamples, targetSignalFrequency, ...
        LPFBandwidthRatio, LPFOrder, thdSamplingRate, numHarmonics, ...
        "Se", pS, OSR, pSe);

    allResults = [allResults; resultsS; resultsSe]; %#ok<AGROW>

    summary = [summary; summarizeMetrics(metricsS, "S", pS, OSR, pSe, dacClock, ...
        targetTHDPercent, targetErrorPercent)]; %#ok<AGROW>
    summary = [summary; summarizeMetrics(metricsSe, "Se", pS, OSR, pSe, dacClock, ...
        targetTHDPercent, targetErrorPercent)]; %#ok<AGROW>

    nexttile(1);
    hold on;
    plot(RList, metricsS.THDPercent, '-', 'LineWidth', 1.1, ...
        'DisplayName', sprintf('S(t), p=%d, OSR=%d', pS, OSR));
    plot(RList, metricsSe.THDPercent, '--', 'LineWidth', 1.1, ...
        'DisplayName', sprintf('Se(t), OSR=%d', OSR));
    ylabel('THD (%)');
    grid on;

    nexttile(2);
    hold on;
    plot(RList, metricsS.ErrorPercent, '-', 'LineWidth', 1.1, ...
        'DisplayName', sprintf('S(t), p=%d, OSR=%d', pS, OSR));
    plot(RList, metricsSe.ErrorPercent, '--', 'LineWidth', 1.1, ...
        'DisplayName', sprintf('Se(t), OSR=%d', OSR));
    xlabel('DAC Quantization Steps R');
    ylabel('Weighted Error (%)');
    grid on;
end

nexttile(1);
yline(targetTHDPercent, ':k', 'Target THD');
title('THD from Original thd() Logic');
legend('Location', 'northeastoutside');

nexttile(2);
yline(targetErrorPercent, ':k', 'Target Error');
title('Weighted Harmonic Error from Original Logic');
legend('Location', 'northeastoutside');

disp('Summary of recommended R values:');
disp(summary);

writetable(allResults, 'StOptimization_conservative_all_results.csv');
writetable(summary, 'StOptimization_conservative_summary.csv');
save('StOptimization_conservative_results.mat', 'allResults', 'summary');

disp('Saved: StOptimization_conservative_all_results.csv');
disp('Saved: StOptimization_conservative_summary.csv');
disp('Saved: StOptimization_conservative_results.mat');

%% Local functions
function waveform = generateSWaveform(time, f, p, OSR)
    phase = 2*pi*f*time;
    waveform = square(phase);

    for i = 1:p
        phaseShift = i*2*pi/OSR;
        weight = cos(phaseShift);
        waveform = waveform ...
            + weight*square(phase + phaseShift) ...
            + weight*square(phase - phaseShift);
    end
end

function waveform = generateSeWaveform(time, f, pSe, OSR)
    phase = 2*pi*f*time;
    waveform = zeros(size(time));

    for i = 1:pSe
        phaseShift = i*2*pi/OSR - pi/OSR;
        weight = cos(phaseShift);
        waveform = waveform ...
            + weight*square(phase + phaseShift) ...
            + weight*square(phase - phaseShift);
    end
end

function [metrics, resultTable] = evaluateWaveformAcrossR( ...
    waveform, RList, dacClock, holdSamples, signalFrequency, ...
    LPFBandwidthRatio, LPFOrder, thdSamplingRate, numHarmonics, ...
    waveformName, pS, OSR, pSe)

    peakToPeak = max(waveform) - min(waveform);
    effectiveSamplingRate = dacClock * holdSamples;
    cutoffFrequency = LPFBandwidthRatio * signalFrequency;
    [filterB, filterA] = butter(LPFOrder, cutoffFrequency/(effectiveSamplingRate/2), 'low');

    THDPercent = zeros(numel(RList), 1);
    ErrorPercent = zeros(numel(RList), 1);
    RequiredBits = zeros(numel(RList), 1);

    for rIndex = 1:numel(RList)
        R = RList(rIndex);
        RequiredBits(rIndex) = ceil(log2(R + 1));

        scaledWaveform = R * (waveform/peakToPeak + 0.5);
        quantizedWaveform = round(scaledWaveform);

        expandedWaveform = repelem(quantizedWaveform, holdSamples);
        filteredWaveform = filter(filterB, filterA, expandedWaveform);

        numSamples = length(filteredWaveform);
        midStartPoint = floor(numSamples * 0.25) + 1;
        midEndPoint = floor(numSamples * 0.75);
        midSectionSignal = filteredWaveform(midStartPoint:midEndPoint);

        harmonicPercent = zeros(1, numHarmonics);
        harmonicError = zeros(1, numHarmonics);

        % This preserves the original code's incremental thd() logic.
        for harmonicOrder = 3:numHarmonics
            thdCurrent = thd(midSectionSignal, thdSamplingRate, harmonicOrder);
            thdCurrentPercent = 100 * (10^(thdCurrent/20));

            thdPrevious = thd(midSectionSignal, thdSamplingRate, harmonicOrder - 1);
            thdPreviousPercent = 100 * (10^(thdPrevious/20));

            differenceSquared = (thdCurrentPercent/100)^2 - (thdPreviousPercent/100)^2;
            harmonicPercent(harmonicOrder) = sqrt(max(differenceSquared, 0)) * 100;
            harmonicError(harmonicOrder) = harmonicPercent(harmonicOrder) / harmonicOrder;
        end

        totalTHD = thd(midSectionSignal, thdSamplingRate, numHarmonics);
        THDPercent(rIndex) = 100 * (10^(totalTHD/20));
        ErrorPercent(rIndex) = sum(harmonicError(3:numHarmonics));
    end

    metrics = table(RList(:), RequiredBits, THDPercent, ErrorPercent, ...
        'VariableNames', {'R', 'RequiredBits', 'THDPercent', 'ErrorPercent'});

    waveformColumn = repmat(waveformName, numel(RList), 1);
    pSColumn = repmat(pS, numel(RList), 1);
    pSeColumn = repmat(pSe, numel(RList), 1);
    osrColumn = repmat(OSR, numel(RList), 1);
    clockColumn = repmat(dacClock/1e6, numel(RList), 1);

    resultTable = table(waveformColumn, pSColumn, pSeColumn, osrColumn, clockColumn, ...
        RList(:), RequiredBits, THDPercent, ErrorPercent, ...
        'VariableNames', {'Waveform','p_for_S','p_for_Se','OSR','DAC_clock_MHz', ...
        'R','Required_bits','THD_percent','Weighted_error_percent'});
end

function oneSummary = summarizeMetrics(metrics, waveformName, pS, OSR, pSe, dacClock, ...
    targetTHDPercent, targetErrorPercent)

    valid = find(metrics.THDPercent <= targetTHDPercent ...
        & metrics.ErrorPercent <= targetErrorPercent, 1, 'first');

    if isempty(valid)
        firstRMeetingTargets = NaN;
        bitsAtTarget = NaN;
        thdAtTarget = NaN;
        errorAtTarget = NaN;
    else
        firstRMeetingTargets = metrics.R(valid);
        bitsAtTarget = metrics.RequiredBits(valid);
        thdAtTarget = metrics.THDPercent(valid);
        errorAtTarget = metrics.ErrorPercent(valid);
    end

    [minTHD, minTHDIndex] = min(metrics.THDPercent);
    [minError, minErrorIndex] = min(metrics.ErrorPercent);

    oneSummary = table(waveformName, pS, pSe, OSR, dacClock/1e6, ...
        firstRMeetingTargets, bitsAtTarget, thdAtTarget, errorAtTarget, ...
        metrics.R(minTHDIndex), minTHD, metrics.R(minErrorIndex), minError, ...
        'VariableNames', {'Waveform','p_for_S','p_for_Se','OSR','DAC_clock_MHz', ...
        'First_R_meeting_targets','Required_bits','THD_at_target_R_percent', ...
        'Error_at_target_R_percent','R_for_min_THD','Min_THD_percent', ...
        'R_for_min_error','Min_error_percent'});
end