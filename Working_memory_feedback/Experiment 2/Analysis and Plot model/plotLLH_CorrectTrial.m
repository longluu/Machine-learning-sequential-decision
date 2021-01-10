%%%%%%%%%%%%%%%%%%%%%%% Plot LLH and fit params from correct trials %%%%%%%%%%%%%%%%%%%%%%%

subjID = {'ll', 'pw', 'eh', 'bh', 'ln', 'at', 'dh'}; % 

%% No resample
paramModel = NaN(length(subjID), 9);
negLLH = NaN(length(subjID), 1);
path_fitResult = 'C:\Users\kwsl455\Machine-learning-sequential-decision\Working_memory_feedback\Experiment 2\Model fit\Fit result\Version 2 (free pC)\NoResample\';
for kk = 1 : length(subjID)
    fileName = [path_fitResult 'FitResult-' subjID{kk} '-extracted.txt'];
    fileID = fopen(fileName);
    paramAll = textscan(fileID,'%f %f %f %f %f %f %f %f %f %f','CommentStyle','//');
    result_mat = cell2mat(paramAll);
    paramModel(kk, :) = result_mat(end, 2:end);
    negLLH(kk, :) = result_mat(end, 1);
    fclose(fileID);
end

noiseSensoryExp1 = paramModel(:, 1:2);
noiseMemoryExp1 = paramModel(:, 6);
noiseAll = [noiseSensoryExp1 noiseMemoryExp1];

% Plot the parameters with bars
figure;
hold on
subplot(2, 2, 1);
bar(noiseAll)
ylim([0 19])
xlabel ('Subject')
ylabel('Noise SD (deg)')
box off

subplot(2, 2, 2);
bar(negLLH)
ylim([0 6000])
xlabel ('Subject')
ylabel('-LLH')
box off

%% No resample, m_m sampled from theta
paramModel = NaN(length(subjID), 9);
negLLH = NaN(length(subjID), 1);
path_fitResult = 'C:\Users\kwsl455\Machine-learning-sequential-decision\Working_memory_feedback\Experiment 2\Model fit\Fit result\Version 2 (free pC)\NoResample_MmFromTheta\';
for kk = 1 : length(subjID)
    fileName = [path_fitResult 'FitResult-' subjID{kk} '-extracted.txt'];
    fileID = fopen(fileName);
    paramAll = textscan(fileID,'%f %f %f %f %f %f %f %f %f %f','CommentStyle','//');
    result_mat = cell2mat(paramAll);
    paramModel(kk, :) = result_mat(end, 2:end);
    negLLH(kk, :) = result_mat(end, 1);
    fclose(fileID);
end

noiseSensoryExp1 = paramModel(:, 1:2);
noiseMemoryExp1 = paramModel(:, 6);
noiseAll = [noiseSensoryExp1 noiseMemoryExp1];

% Plot the parameters with bars
subplot(2, 2, 3);
bar(noiseAll)
ylim([0 19])
xlabel ('Subject')
ylabel('Noise SD (deg)')
box off

subplot(2, 2, 4);
bar(negLLH)
ylim([0 6000])
xlabel ('Subject')
ylabel('-LLH')
box off