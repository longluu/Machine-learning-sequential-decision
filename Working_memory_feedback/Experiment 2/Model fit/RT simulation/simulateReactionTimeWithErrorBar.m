%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulate reaction time of resample model %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract model parameter
subjectIDAll = {'average'};
params_ave = [6.2041    8.3568           0.0000     31.4168   -19.0927   5.4034    2.6349    0.9999    0.5358];
paramsBootstrap = cell(1, length(subjectIDAll));
maxIter = 0;
subjIter = NaN(1, length(subjectIDAll));
for kk = 1 : length(subjectIDAll);
    fileName = ['FitResult-Bootstrap-' subjectIDAll{kk} '.txt'];
    fileID = fopen(fileName);
    paramAll = textscan(fileID,'%f %f %f %f %f %f %f %f %f %f','CommentStyle','//');
    paramsBootstrap{kk} = cell2mat(paramAll(2:end));
    subjIter(kk) = size(paramsBootstrap{kk}, 1);
    if maxIter < subjIter(kk)
        maxIter = subjIter(kk);
    end    
    fclose(fileID);
    plotFig = 0;
    if plotFig
        h = figure;
        figPos = [0.01, 0.2, 0.98, 0.6];
        set(h,'Units','normalized','Position',figPos)
        hold on
        fontSize = 20;
        set(gca, 'FontSize', fontSize)
    end
end
paramsBootstrap = paramsBootstrap{1};
paramsBootstrap = [paramsBootstrap; params_ave];

%% Compute the resampling time for all boostrapped samples
rt_simulate_correct_ave = NaN(size(paramsBootstrap, 1), 2);
rt_simulate_incorrect_ave = NaN(size(paramsBootstrap, 1), 2);

for bb = 1 : size(paramsBootstrap, 1)
    stdSensory = paramsBootstrap(bb, 1:2);
    stdMemory = paramsBootstrap(bb, 6);
    stdResample = sqrt(stdSensory.^2 + stdMemory^2);
    thetaStim = [-12:2:-2 5:5:30];
    dstep = 0.1;
    rangeM = round([thetaStim(1)-7*stdSensory(2) thetaStim(end)+7*stdSensory(2)], -log10(dstep));
    m_original = rangeM(1):dstep:rangeM(2);
    rangeMR = round([min(rangeM)-7*stdResample(2) max(rangeM)+7*stdResample(2)], -log10(dstep));     
    mr = rangeMR(1):dstep:rangeMR(2);

    %% Correct trials
    rt_simulate_correct = NaN(length(stdSensory), length(thetaStim));
    for ii = 1 : length(stdSensory)
        for jj = 1 : length(thetaStim)
            % Sensory distribution p(m|theta)
            indCorrect = sign(m_original) == sign(thetaStim(jj));
            m = m_original(indCorrect);
            pmGtheta = exp(-((m - thetaStim(jj)).^2) / (2*stdSensory(ii)^2));
            pmGtheta = pmGtheta / sum(pmGtheta);

            % Resample distribution p(mr|m)
            pmrGm = zeros(length(m), length(mr));
            pmrGm(1, :) = exp(-((mr - m(1)).^2) / (2*stdResample(ii)^2));
            pmrGm(1, :) = pmrGm(1, :) / sum(pmrGm(1, :));
            for kk = 2 : length(m)
                pmrGm(kk, kk:end) = pmrGm(1, 1:end-kk+1);
            end
            indCorrect = sign(mr) == sign(thetaStim(jj));
            pmrGm(:, indCorrect) = 0;
            pIncorrect = sum(pmrGm, 2);
            rt_simulate_correct(ii, jj) = pmGtheta * pIncorrect;
        end
    end
    rt_simulate_correct_ave(bb, :) = mean(rt_simulate_correct, 2);
    
    %% Incorrect trials
    rt_simulate_incorrect = NaN(length(stdSensory), length(thetaStim));
    for ii = 1 : length(stdSensory)
        for jj = 1 : length(thetaStim)
            % Sensory distribution p(m|theta)
            indIncorrect = sign(m_original) ~= sign(thetaStim(jj));
            m = m_original(indIncorrect);
            pmGtheta = exp(-((m - thetaStim(jj)).^2) / (2*stdSensory(ii)^2));
            pmGtheta = pmGtheta / sum(pmGtheta);

            % Resample distribution p(mr|m)
            pmrGm = zeros(length(m), length(mr));
            pmrGm(1, :) = exp(-((mr - m(1)).^2) / (2*stdResample(ii)^2));
            pmrGm(1, :) = pmrGm(1, :) / sum(pmrGm(1, :));
            for kk = 2 : length(m)
                pmrGm(kk, kk+1:end) = pmrGm(1, 1:end-kk);
            end
    %         if thetaStim(jj) == 10
    %             keyboard
    %         end        
            indCorrect = sign(mr) == sign(thetaStim(jj));
            pmrGm(:, indCorrect) = 0;
            pIncorrect = sum(pmrGm, 2);
            rt_simulate_incorrect(ii, jj) = pmGtheta * pIncorrect;
        end
    end
    rt_simulate_incorrect_ave(bb, :) = mean(rt_simulate_incorrect, 2);
end

%% Compute the 95% CI
lower_CI = NaN(2, length(stdSensory));
upper_CI = NaN(2, length(stdSensory));
for ii = 1 : length(stdSensory)
    ci_correct = prctile(rt_simulate_correct_ave(1:end-1, ii),[2.5 97.5]);
    ci_incorrect = prctile(rt_simulate_incorrect_ave(1:end-1, ii),[2.5 97.5]);
    
    lower_CI(1, ii) = ci_correct(1);
    lower_CI(2, ii) = ci_incorrect(1);
    upper_CI(1, ii) = ci_correct(2);
    upper_CI(2, ii) = ci_incorrect(2);
end

%% Plot the result
color_rgb = [   0    0.4470    0.7410;
                0.8500    0.3250    0.0980;
                0.9290    0.6940    0.1250;
                0.4940    0.1840    0.5560;
                0.4660    0.6740    0.1880;
                0.3010    0.7450    0.9330;
                0.6350    0.0780    0.1840];

figure
rt_ave = [rt_simulate_correct_ave(end, :); rt_simulate_incorrect_ave(end, :)];
errorBarGraph(rt_ave, lower_CI, upper_CI, color_rgb)
xlabel('Trial type')
ylabel('Resampling time')

