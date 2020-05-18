%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulate reaction time of resample model %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model parameter
params_ave = [6.2041    8.3568           0.0000     31.4168   -19.0927   5.4034    2.6349    0.9999    0.5358];

stdSensory = paramsAll(1:2);
stdMemory = paramsAll(6);
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

%% Plot the result
color_rgb = [   0    0.4470    0.7410;
                0.8500    0.3250    0.0980;
                0.9290    0.6940    0.1250;
                0.4940    0.1840    0.5560;
                0.4660    0.6740    0.1880;
                0.3010    0.7450    0.9330;
                0.6350    0.0780    0.1840];
figure
h_legend = NaN(1, 2);
subplot(1, 2, 1)
hold on
for ii = 1 : 2
    h_legend(ii) = plot(thetaStim(1:6), rt_simulate_correct(ii, 1:6), 'o-', 'Color', color_rgb(ii, :));
    plot(thetaStim(7:end), rt_simulate_correct(ii, 7:end), 'o-', 'Color', color_rgb(ii, :))
end
xlabel('Stimulus orientation (deg)')
ylabel('Resampling time')
title('Correct trials')
legend(h_legend, 'Low noise', 'High noise')

subplot(1, 2, 2)
hold on
for ii = 1 : 2
    h_legend(ii) = plot(thetaStim(1:6), rt_simulate_incorrect(ii, 1:6), 'o-', 'Color', color_rgb(ii, :));
    plot(thetaStim(7:end), rt_simulate_incorrect(ii, 7:end), 'o-', 'Color', color_rgb(ii, :))
end
xlabel('Stimulus orientation (deg)')
ylabel('Resampling time')
title('Incorrect trials')

figure
rt_ave_correct = mean(rt_simulate_correct, 2);
rt_ave_incorrect = mean(rt_simulate_incorrect, 2);
bar([1 2], [rt_ave_correct'; rt_ave_incorrect'])
xlabel('Trial type')
ylabel('Resampling time')

