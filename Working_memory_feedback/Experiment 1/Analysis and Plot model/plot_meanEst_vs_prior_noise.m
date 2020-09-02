%%%%%%%%%%%%%%%%%%%%%%% Plot the mean estimate vs fit prior and noise %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Fit params from correct - resample %%%%%%%%%%%%%%%%%%%%%%%

%% Compute the subject LLH
subjectIDAll = {'ll', 'an', 'ep', 'jp', 'kc'};
estimate_Data = NaN(length(subjectIDAll), 2, 8);
fit_prior = NaN(length(subjectIDAll), 1);
fit_memory_noise = NaN(length(subjectIDAll), 1);
fit_sensory_noise = NaN(length(subjectIDAll), 2);

paramsAllSubject = [2.6500    6.0895           0.0000     22.2852     1.6506   0.9414    2.0976;
                    3.0023    9.7384           0.0000     34.4053     0.0615   0.9480    3.1069;
                    4.6136   10.4165           0.0000     29.8375     0.1325   0.9940    3.8106;
                    7.7094   11.9114           0.0000     55.7419     0.0083   0.2850    3.8551;
                    5.1033   10.3703           0.0000     46.6421     4.7921   0.8187    3.3313];

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
    stdNoiseLevel = unique(dataAll(:,4));
    angleDiff = (unique(dataAll(:,1)))';
    angleEstimate = dataAll(:, 7);

    [~, ~, ~, estimateData, ~, ~] = dataForFitting(subjectID, 0, 0);

    % Compute meanEstimate
    meanEstimate = NaN(size(estimateData));
    for ii = 1 : size(estimateData, 1)
        for jj = 1 : size(estimateData, 2)
            tempEst = abs(estimateData{ii, jj});
            if sum(~isnan(tempEst)) > 5
                meanEstimate(ii, jj) = nanmean(tempEst);
            end
        end
    end
    estimate_Data(nn, :, :) = meanEstimate;
    
    % Get the model params
    fit_prior(nn) = paramsAllSubject(nn, 4);
    fit_memory_noise(nn) = paramsAllSubject(nn, 5);
    fit_sensory_noise(nn,:) = paramsAllSubject(nn, 1:2);
end

%% Plot the estimates
mean_est_data_subj_sensory = squeeze(nanmean(estimate_Data, 3));
mean_est_data_subj = squeeze(nanmean(mean_est_data_subj_sensory, 2));
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
    plot(mean_est_data_subj(nn), fit_prior(nn), 'o', 'Color', colorIndex(nn, :))
end
xlabel('Mean estimate - data (deg)')
ylabel('Fit prior (deg)')
r = round(corr(mean_est_data_subj, fit_prior, 'type', 'Pearson'), 2);
title (['Prior, r: ' num2str(r)])
xlim([3 14])

subplot(1, 4, 2)
hold on
set(gca, 'FontSize', 12)
for nn = 1 : length(subjectIDAll)
    plot(mean_est_data_subj_sensory(nn, 1), fit_sensory_noise(nn, 1), 'o', 'Color', colorIndex(nn, :))
    plot(mean_est_data_subj_sensory(nn, 2), fit_sensory_noise(nn, 2), 'x', 'Color', colorIndex(nn, :))
end
xlabel('Mean estimate - data (deg)')
ylabel('Fit sensory noise (deg)')
r = round(corr(mean_est_data_subj_sensory(:), fit_sensory_noise(:), 'type', 'Pearson'), 2);
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
xlim([3 14])
ylim([-1 5])

subplot(1, 4, 4)
hold on
set(gca, 'FontSize', 12)
for nn = 1 : length(subjectIDAll)
    plot(mean_est_data_subj_sensory(nn, 1), fit_sensory_memory_noise(nn, 1), 'o', 'Color', colorIndex(nn, :))
    plot(mean_est_data_subj_sensory(nn, 2), fit_sensory_memory_noise(nn, 2), 'x', 'Color', colorIndex(nn, :))
end
xlabel('Mean estimate - data (deg)')
ylabel('Fit sensory + memory noise (deg)')
r = round(corr(mean_est_data_subj_sensory(:), fit_sensory_memory_noise(:), 'type', 'Pearson'), 2);
title (['Sensory + memory noise, r: ' num2str(r)])