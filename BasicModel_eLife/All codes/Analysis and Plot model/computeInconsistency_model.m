% conditional observer model with memory recall and motor noise
% for orientation stimulus (infinite space approximation) 
% ********** Old version - condition the prior only ********** 
% astocker - lluu
% 06.2017
flagSC = 0; % 1: self-conditioned model
            % 0: standard Bayes
flagDecisionGiven = 0; % 1: stimulus category (cw/ccw) are given (Exp. 3)
                       % 0: observer makes categorical decision (Exp. 1 and 2)
includeIncongruentTrials = 'include';

dstep = 0.1;
subjectID = {'Sc1', 'Sc2', 'll', 'sy', 'cz', 'vs', 'as', 'll', 'xfl', 'aj', 'zw', 'skb'};
paramsAllSubject = [3.4994    5.4317    9.5625      0.0000    25        6.0913    0.2953    0.01;
             4.6025    6.1782    8.9676      0.0000    23.1643   5.8385    0.8455    0.01;
             2.5611    4.8570    7.1153      0.0000    28.6016      0.8805    0.8461   0.01;
             4.0300    5.3505    8.1590      0.0000    41.7649   2.0957    0.0257    0.01;
             3.8692    6.4042   12.3882      0.0000    39.5010  14.0727    0.2322    0.01;
             3.7041    4.6566   10.9916      0.0000    36.4739   0.2433    0.6481    0.01;
             3.9226    5.7843   10.9330      0.0000    39.0534  14.2463    0.1492    0.01;
             2.5611    4.8570    7.1153      0.0000    23.1044    0.8805    0.8461    0.01;
             6.3630    8.4075   14.5464      0.0000    22.5672   5.9826    0.7652    0.01;
             4.6736    6.1762    7.9466      0.0000    15.8459   1.2144    0.3681    0.01;
             4.6585    5.5272    6.8023      0.0000    17.4750   4.1846    0.0934    0.01;
             4.3172    4.5344    8.4083      0.0000    33.9063  14.5556    0.5305    0.01];

thetaStim = -21:0.1:21;  
thetaStim = round(thetaStim, -log10(dstep));
         
figure;
fontSize = 15;
percentInconsistent = NaN(length(subjectID), 3, length(thetaStim));
percentInconsistentTotal = NaN(length(subjectID), 3);
colorName = {'Black', 'Brown', 'RosyBrown', 'SlateGray', 'Crimson', 'DarkMagenta', 'DarkOrange', 'DarkGoldenRod', 'SpringGreen', 'Teal', 'DodgerBlue', 'Navy'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

for ii = 1 : length(subjectID)
    paramsAll = paramsAllSubject(ii, :);
    % lapse
    lapseRate = paramsAll(4);


    % sensory noise
    stdSensory = paramsAll(1:3);

    % memory recall noise
    stdMemory = paramsAll(6);

    % motor noise;
    stdMotor = paramsAll(8);

    % priors
    pC = [0.5, 0.5]'; % [cw ccw]
    pthcw = paramsAll(5);
    pthccw = -paramsAll(5);
    smoothFactor = paramsAll(7);
    
    rangeth = [-42 42];
    th = rangeth(1):dstep:rangeth(2);
    th = round(th, -log10(dstep));
    nth = length(th);

    pthGC = zeros(2,nth);
    pthGC_Discrimination = zeros(2,nth);
    if flagSC
        pthGC(1,:) = TukeyWindow([0 pthcw], 0, smoothFactor, th);
        pthGC(2,:) = TukeyWindow([pthccw 0], 1, smoothFactor, th);
    else
        pth = (TukeyWindow([0 pthcw], 0, smoothFactor, th) + TukeyWindow([pthccw 0], 1, smoothFactor, th))/2;
        pth(th==0) = 0;
        pth(th==0) = max(pth);
        pthGC(1,:) = pth;
        pthGC(2,:) = pth;
    end
    pthGC_Discrimination(1,:) = TukeyWindow([0 pthcw], 0, smoothFactor, th);
    pthGC_Discrimination(2,:) = TukeyWindow([pthccw 0], 1, smoothFactor, th);
    
    for kk=1:length(stdSensory)  
        rangeM = [min(thetaStim)-5*stdSensory(kk) max(thetaStim)+5*stdSensory(kk)];
        if rangeM(2) < rangeth(2)
            rangeM = rangeth;
        end
        nm = 1000;
        m = linspace(rangeM(1), rangeM(2), nm);

        nmm = 1200;
        rangeMM = [min(rangeM)-6*stdMemory max(rangeM)+6*stdMemory];
        if rangeMM(2) < rangeth(2)
            rangeMM = rangeth;
        end        
        mm = linspace(rangeMM(1), rangeMM(2), nmm);

        M = repmat(m',1,nth);
        MM_m = repmat(mm',1,nm);
        MM_th = repmat(mm',1,nth); 
        THm = repmat(th, nm, 1); 
        THmm = repmat(th, nmm, 1); 

        %% Generative (forward)
        % orientation noise p(m|th)
        pmGth = exp(-((M-THm).^2)./(2*stdSensory(kk)^2));
        pmGth = pmGth./(repmat(sum(pmGth,1),nm,1)); 

        %% Inference
        % 1: categorical judgment
        if ~flagDecisionGiven
            PCGm = (pthGC_Discrimination * pmGth') .* repmat(pC,1,nm);
            % fix the issue when sensory noise is too low
            indFirstNonZero = find(PCGm(2,:), 1);
            PCGm(2, 1: indFirstNonZero-1) = PCGm(2, indFirstNonZero);
            indLastNonZero = find(PCGm(1,:), 1, 'last');
            PCGm(1, indLastNonZero+1:end) = PCGm(1, indLastNonZero);
            PCGm = PCGm./(repmat(sum(PCGm,1),2,1));
            % max posterior decision
            PChGm = round(PCGm);
            % marginalization
            PChGtheta = PChGm * pmGth(:, ismember(th, thetaStim));
            PChGtheta_lapse = lapseRate + (1 - 2*lapseRate) * PChGtheta;
            PChGtheta_lapse = PChGtheta_lapse ./ repmat(sum(PChGtheta_lapse, 1), 2, 1);
        else
            PChGtheta_lapse = NaN(2, length(thetaStim));
            PChGtheta_lapse(1, thetaStim > 0) = 1;
            PChGtheta_lapse(1, thetaStim < 0) = 0;
            PChGtheta_lapse(1, thetaStim == 0) = 0.5;
            PChGtheta_lapse(2, :) = 1 - PChGtheta_lapse(1, :);
        end    

        % 2: estimation
        pmmGth = exp(-((MM_th-THmm).^2)./(2*(stdSensory(kk)^2 + stdMemory^2))); % p(mm|th) = N(th, sm^2 + smm^2)
        pmmGth = pmmGth./(repmat(sum(pmmGth,1),nmm,1)); 
        pthGmmChcw = (pmmGth.*repmat(pthGC(1,:),nmm,1))';
        pthGmmChcw = pthGmmChcw./repmat(sum(pthGmmChcw,1),nth,1);
        pthGmmChcw(isnan(pthGmmChcw)) = 0;

        pthGmmChccw = (pmmGth.*repmat(pthGC(2,:),nmm,1))';
        pthGmmChccw = pthGmmChccw./repmat(sum(pthGmmChccw,1),nth,1);
        pthGmmChccw(isnan(pthGmmChccw)) = 0;

        EthChcw = th * pthGmmChcw;
        EthChccw = th * pthGmmChccw;
        % discard repeating/decreasing values (required for interpolation) 
        indKeepCw = 1:length(EthChcw);
        while sum(diff(EthChcw)<=0) >0
            indDiscardCw = [false diff(EthChcw)<=0];
            EthChcw(indDiscardCw) = [];
            indKeepCw(indDiscardCw) = [];
        end
        indKeepCcw = 1:length(EthChccw);
        while sum(diff(EthChccw)<=0) >0
            indDiscardCcw = [diff(EthChccw)<=0 false];
            EthChccw(indDiscardCcw) = [];
            indKeepCcw(indDiscardCcw) = [];
        end

        a = 1./gradient(EthChcw,dstep);
        if ~flagDecisionGiven
            % memory noise
            pmmGm = exp(-((MM_m-repmat(m, nmm, 1)).^2)./(2*stdMemory^2)); 
            pmmGm = pmmGm./(repmat(sum(pmmGm,1),nmm,1));   

            % attention marginalization: compute distribution only over those ms that lead to cw decision!
            pmmGthChcw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(1,:)',1,length(thetaStim)));
            pmmGthChcw = pmmGthChcw ./ repmat(sum(pmmGthChcw,1),nmm,1);
            b = repmat(a',1,length(thetaStim)) .* pmmGthChcw(indKeepCw, :);        
        else
            b = repmat(a',1,length(thetaStim)) .* pmmGth(indKeepCw, ismember(th, thetaStim));   
        end
        pthhGthChcw = interp1(EthChcw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChcw = conv2(pthhGthChcw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChcw(pthhGthChcw < 0) = 0; 

        a = 1./gradient(EthChccw,dstep);
        if ~flagDecisionGiven
            % attention marginalization: compute distribution only over those ms that lead to cw decision!
            pmmGthChccw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(2,:)',1,length(thetaStim)));        
            b = repmat(a',1,length(thetaStim)) .* pmmGthChccw(indKeepCcw, :);        
        else
            b = repmat(a',1,length(thetaStim)) .* pmmGth(indKeepCcw, ismember(th, thetaStim));
        end  
        pthhGthChccw = interp1(EthChccw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChccw = conv2(pthhGthChccw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChccw(pthhGthChccw < 0) = 0; 
        
        PChGtheta_lapse_new = PChGtheta_lapse;
        pthGIndexcongruent = ones(1, size(pthhGthChcw, 2));

        pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
        pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);        
        
        pthhGthChcw = pthhGthChcw.*repmat(PChGtheta_lapse_new(1,:),nth,1);
        pthhGthChccw = pthhGthChccw.*repmat(PChGtheta_lapse_new(2,:),nth,1);
        
        pthhGthChcw_norm = pthhGthChcw / sum(pthhGthChcw(:));
        
        %% Compute the inconsistency fraction
        percentInconsistent(ii, kk, :) = 100 * (sum(pthhGthChcw(th<=0, :), 1) + sum(pthhGthChccw(th>=0, :), 1));
        percentInconsistentTotal(ii, kk) = 100 * sum(sum(pthhGthChcw_norm(th<=0, :)));
        
        %% plot
        %     showrange = [-21 21];
        %     ind = find(thetaStim >= showrange(1) & thetaStim <= showrange(2));
        %     nthshow = length(ind);
        %     
        %     subplot(1, 3, kk);
        %     pthhGthChcw = max(pthhGthChcw(:)) - pthhGthChcw;
        %     xRange = [-22 22];
        %     indX = find(thetaStim >= xRange(1) & thetaStim <= xRange(2));
        %     xMax = length(indX);
        %     yRange = [-40 40];
        %     indY = find(th >= yRange(1) & th <= yRange(2));
        %     yMax = length(indY);
        %     thNew = th(indY);
        %     indYStart = find(thNew == xRange(1));
        %     indYEnd = find(thNew == xRange(2));
        %     imagesc(pthhGthChcw(indY, indX));
        %     hold on;
        %     axis xy;
        %     colormap('gray');
        %     plot([1 xMax],[round(yMax/2) round(yMax/2)],'k:', 'LineWidth', 1);
        %     plot([round(xMax/2) round(xMax/2)],[1 yMax],'k:', 'LineWidth', 1);
        %     plot([1 xMax],[indYStart indYEnd],'w--', 'LineWidth', 1.5);
        %     set(gca, 'ylim', [1 yMax], 'xlim', [1 xMax], ...
        %         'XTick', round(linspace(1,xMax,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
        %         'YTick', round(linspace(1,yMax,5)), 'YTickLabel', num2cell(round(linspace(yRange(1),yRange(2),5))), ...
        %         'FontSize', fontSize)

        subplot(2, 3, kk)
        hold on
        indexPlot = thetaStim >= 0;
        plot(thetaStim(indexPlot), squeeze(percentInconsistent(ii, kk, indexPlot)), 'Color', colorIndex(ii, :))
        xlim([min(thetaStim(indexPlot)) max(thetaStim(indexPlot))])
        xlabel('Stimulus orientation (deg)')
        ylabel('Percent inconsistent (%)')
    end
end

% Plot the legend
subplot(2, 3, 4)
hold on
set(gca, 'FontSize', 30)
nullX = zeros(1, length(subjectID));
nullY = zeros(1, length(subjectID));
h = NaN(1, length(subjectID));
for ii = 1 : length(subjectID)
    h(ii) = scatter(nullX(ii), nullY(ii), 1, colorIndex(ii, :), 'filled');
end
[~, iconHandle]  = legend(h,'Sc1', 'Sc2', 'S1', 'S2', 'S3', 'S4', 'S5', 'S1', 'S6', 'S7', 'S8', 'S9', 'Location', 'NorthWest');
legend('boxoff')
for ii = 1 : length(subjectID)
    iconHandle(ii+length(subjectID)).Children.MarkerSize = 15;
end

% Plot the total inconsistency
hFig = figure;
hAx1 = gca;
[~, hPanel] = errorbar_groups(percentInconsistentTotal', zeros(size(percentInconsistentTotal')), zeros(size(percentInconsistentTotal')), ...
                'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig, 'AxID', hAx1);
set(hPanel, 'EdgeColor', 'none')
set(gca, 'FontSize', fontSize)
xlabel('Subject')
ylabel('Percent inconsistent (%)')

% Plot the memory noise vs. inconsistency
memNoise = paramsAllSubject(3:end, 6);
percentInconsistentAll = mean(percentInconsistentTotal(3:end, :), 2);

figure;
colorName = {'SlateGray', 'Crimson', 'DarkMagenta', 'DarkOrange', 'DarkGoldenRod', 'SlateGray', 'SpringGreen', 'Teal', 'DodgerBlue', 'Navy'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
regressSlope = regress(percentInconsistentAll, [ones(length(memNoise), 1) memNoise]);
Xregress = [min(memNoise) max(memNoise)];
Yregress = Xregress * regressSlope(2) + regressSlope(1);

hold on
set(gca, 'FontSize', 20)
scatter(memNoise, percentInconsistentAll, 13^2, colorIndex, 'filled');
plot(Xregress, Yregress)
xlabel('Memory noise (deg)')
ylabel('Percent inconsistent (predicted)')
title([' r: ', num2str(corr(memNoise, percentInconsistentAll))])


