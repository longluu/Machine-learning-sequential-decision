%%%%%%%%%%%%%%%%%%%%%%% Plot the mean estimate vs fit prior and noise %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Fit params from correct - resample %%%%%%%%%%%%%%%%%%%%%%%

%% Compute the subject LLH
subjectIDAll = {'ll', 'pw', 'eh', 'bh', 'ln', 'at', 'dh'}; % 
estimate_Data_cw = NaN(length(subjectIDAll), 2);
fit_prior = NaN(length(subjectIDAll), 1);
fit_memory_noise = NaN(length(subjectIDAll), 1);
fit_sensory_noise = NaN(length(subjectIDAll), 2);

paramsAllSubject = [5.1714    6.5149           0.0000     18.2531    -9.3607   0.9940    2.0902    0.9978    0.6635;
                    8.5063    8.7485           0.0000     33.2603   -21.0845   0.6393    1.8928    0.9098    0.4518;
                    8.7284    8.9988           0.0000     13.8798   -12.3746   0.6037    2.7094    0.9189    0.5136;
                    6.1961    8.3406           0.0000     19.5003   -12.7209   0.7809    2.6041    0.4987    0.4168;
                    5.1595    8.9406           0.0000     32.3917   -17.5501   5.2297    1.5830    0.3252    0.3550;
                    10.0963   15.0721           0.0000     57.4518   -39.0410   0.0101    4.0136    0.4339    0.5875;
                    7.7473    9.7957           0.0000     23.1824   -18.2993   0.6654    3.9069    0.8918    0.4908];

for nn = 1 : length(subjectIDAll)
    subjectID = subjectIDAll{nn};
    if strcmp(subjectID, 'average')
        experimentNumber = 1:5;
    else
        experimentNumber = 1;    
    end
    experimentType = 'MainExperiment';
    experiment = 'Original';
    
    session = 1;
    dataAll = [];
    fontSize = 20;
    lineWidth = 2;

    for ii = 1 : length(experimentNumber)
        if strcmp(experiment, 'PilotData')
            dataFullPath = fullfile('Data', subjectID, experimentType, experiment, ['Session' num2str(session)], ...
                                        ['ConditionDecisionMaking-' num2str(experimentNumber(ii))]);
        elseif strcmp(experimentType, 'MotorNoise')
            dataFullPath = fullfile('Data', subjectID, experimentType, ['Session' num2str(session)], ...
                                        ['MotorNoise-' num2str(experimentNumber(ii))]);
        elseif strcmp(experimentType, 'PerceptNoise')
            dataFullPath = fullfile('Data', subjectID, experimentType, ['Session' num2str(session)], ...
                                        ['PerceptNoise-' num2str(experimentNumber(ii))]);                                
        elseif strcmp(experimentType, 'TrainArray')
            dataFullPath = fullfile('Data', subjectID, experimentType, ['Session' num2str(session)], ...
                                        ['TrainArray-' num2str(experimentNumber(ii))]);  
        else
            dataFullPath = fullfile('Data', subjectID, experimentType, [experiment num2str(session)], ...
                                        [experiment '-' num2str(experimentNumber(ii))]);          
        end
        load(dataFullPath);
        if ii == 1
            dataAll = dataTotal;
            if strcmp(experimentType, 'PerceptNoise') && isfield(params, 'staircaseTrack')
                staircaseTrack = params.staircaseTrack;
            end
        else
            dataAll = [dataAll; dataTotal];
            if strcmp(experimentType, 'PerceptNoise') && isfield(params, 'staircaseTrack')
                staircaseTrack = [staircaseTrack; params.staircaseTrack];
            end
        end
    end
    biasMotor = zeros(1,15);  
    if strcmp(experimentType, 'MotorNoise')
        a = dataAll(:,7);
        dataAll(isnan(a),:) = [];    
    end

    % Extract data
    %  Column 1: angle differences (theta1 - theta2)
    %  Column 2: bar presentation time
    %  Column 3: SOA
    %  Column 4: bar Stimulus noise
    %  Column 5: bar reference angle
    %  Column 6: CW/CCW OR Red/Green
    %  Column 7: estimated angle
    %  Column 8: given CW/CCW
    angleDiff = (unique(dataAll(:,1)))';
    angleEstimate = dataAll(:, 7);
    stdNoiseLevel = unique(dataAll(:,4));
    
    angleTrue = dataAll(:,5) - dataAll(:,1);
    angleTrue(angleTrue>180 | angleTrue<0) = angleTrue(angleTrue>180 | angleTrue<0)...
                                                - 180 * sign(angleTrue(angleTrue>180 | angleTrue<0));
    indexFixing = find(abs(angleTrue-angleEstimate) > 70);
    angleEstimate(indexFixing) = angleEstimate(indexFixing) + sign(angleTrue(indexFixing)-angleEstimate(indexFixing))*180;    
    angleDiffEst = dataAll(:, 5) - angleEstimate;
    indexAdjust1 = abs(angleDiffEst)>90;
    indexAdjust2 = sign(angleDiffEst) ~= sign(dataAll(:, 1));
    indexAdjust = indexAdjust1 & indexAdjust2;
    angleDiffEst(indexAdjust) = angleDiffEst(indexAdjust) - 180 * sign(angleDiffEst(indexAdjust));
    biasAll = angleTrue - angleEstimate;
    binaryFeedbackAll = dataAll(:,8);  
    binaryDecisionAll = dataAll(:,6);   
    indicatorConsistent = sign(angleDiffEst) == binaryFeedbackAll;
    angleDiffEst(~indicatorConsistent) = NaN;
    angleDiffEst(binaryDecisionAll == binaryFeedbackAll) = NaN;
    
    % Get the mean estimate for CW
    estimate_Data_cw(nn, 1) = nanmean(angleDiffEst(dataAll(:, 4) == stdNoiseLevel(1) & binaryFeedbackAll == 1));
    estimate_Data_cw(nn, 2) = nanmean(angleDiffEst(dataAll(:, 4) == stdNoiseLevel(2) & binaryFeedbackAll == 1));
    
    % Get the model params
    fit_prior(nn, :) = abs(paramsAllSubject(nn, 4));
    fit_memory_noise(nn) = paramsAllSubject(nn, 6);
    fit_sensory_noise(nn,:) = paramsAllSubject(nn, 1:2);
end

%% Plot the estimates
mean_est_data_subj = squeeze(nanmean(estimate_Data_cw, 2));
fit_sensory_memory_noise = sqrt(fit_sensory_noise.^2 + fit_memory_noise.^2);

figure
colorName = {'SlateGray', 'DarkMagenta', 'DarkOrange', 'DarkGoldenRod', 'SpringGreen', 'Teal', 'DodgerBlue', 'Navy'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

subplot(1, 4, 1)
hold on
set(gca, 'FontSize', 12)
for nn = 1 : length(subjectIDAll)
    plot(mean_est_data_subj(nn), fit_prior(nn, 1), 'o', 'Color', colorIndex(nn, :))
end
xlabel('Mean estimate - data (deg)')
ylabel('Fit prior (deg)')
r = round(corr(mean_est_data_subj, fit_prior, 'type', 'Pearson'), 2);
title (['Prior, r: ' num2str(r)])
xlim([2 14])

subplot(1, 4, 2)
hold on
set(gca, 'FontSize', 12)
for nn = 1 : length(subjectIDAll)
    plot(estimate_Data_cw(nn, 1), fit_sensory_noise(nn, 1), 'o', 'Color', colorIndex(nn, :))
    plot(estimate_Data_cw(nn, 2), fit_sensory_noise(nn, 2), 'x', 'Color', colorIndex(nn, :))
end
xlabel('Mean estimate - data (deg)')
ylabel('Fit sensory noise (deg)')
r = round(corr(estimate_Data_cw(:), fit_sensory_noise(:), 'type', 'Pearson'), 2);
title (['Sensory noise, r: ' num2str(r)])

subplot(1, 4, 3)
hold on
set(gca, 'FontSize', 12)
for nn = 1 : length(subjectIDAll)
    plot(mean_est_data_subj(nn), fit_memory_noise(nn), 'o', 'Color', colorIndex(nn, :))
end
xlabel('Mean estimate - data (deg)')
ylabel('Fit memory noise (deg)')
r = round(corr(mean_est_data_subj, fit_memory_noise, 'type', 'Pearson'), 2);
title (['Memory noise, r: ' num2str(r)])
xlim([2 14])
ylim([-1 6])

subplot(1, 4, 4)
hold on
set(gca, 'FontSize', 12)
for nn = 1 : length(subjectIDAll)
    plot(estimate_Data_cw(nn, 1), fit_sensory_memory_noise(nn, 1), 'o', 'Color', colorIndex(nn, :))
    plot(estimate_Data_cw(nn, 2), fit_sensory_memory_noise(nn, 2), 'x', 'Color', colorIndex(nn, :))
end
xlabel('Mean estimate - data (deg)')
ylabel('Fit sensory + memory noise (deg)')
r = round(corr(estimate_Data_cw(:), fit_sensory_memory_noise(:), 'type', 'Pearson'), 2);
title (['Sensory + memory noise, r: ' num2str(r)])