%%%%%%%% The sampling Bayesian model of conditioned perception %%%%%%%%%%
% ********** Old version - condition the prior only ********** 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
% 


%% Set the parameters of the model
subjectID = 'average';
flagSC = 1; % 1: self-conditioned model
            % 0: standard Bayes
incorrectType = 1; % 1: flip the decision bit
                   % 2: flip the estimates thetaHat
                   % 3: resampling distribution (center: memory sample, width: sensory+memory)
                   % 4: resampling distribution (center: sensory sample, width: memory)
                   % 5: resampling distribution (center: sensory sample, width: sensory+memory)
            
params = [5.1033   10.3703           0.0000     46.6421     0.05   0.8187    3.3313];
stdSensory = params(1:2);
lapseRate = params(3);
stdMemory = params(5);
priorRange = params(4);
stdMotor = params(7);                     
smoothFactor = params(6);
boundaryCutoff = 0;
thetaStim = -21:3:21;
pC = [0.5, 0.5]'; % [cw ccw]
pthcw = priorRange;
pthccw = -priorRange;
nTrialPerCondition = 2000;
fontSize = 20;
fixedConsistencyIndex = 1;
marginalizedVersion = 1;
if fixedConsistencyIndex
    consistencyIndex = 1;
end
includeIncongruentTrials = 0;
rangeCollapse = round(length(thetaStim)/2);

%% Modeling
dstep = 0.1;
rangeth = [-60 60];
th = rangeth(1):dstep:rangeth(2);
nth = length(th);

pthGC = zeros(2,nth);
if flagSC
    pthGC(1,:) = TukeyWindowNew([0 pthcw], smoothFactor, th, boundaryCutoff);
    pthGC(2,:) = TukeyWindowNew([pthccw 0], smoothFactor, th, boundaryCutoff);
else
    pth = (TukeyWindow([0 pthcw], 0, smoothFactor, th) + TukeyWindow([pthccw 0], 1, smoothFactor, th))/2;
    pth(th==0) = 0;
    pth(th==0) = max(pth);
    pthGC(1,:) = pth;
    pthGC(2,:) = pth;    
end

pTheta_CCW = pthGC(2,:);
pTheta_CW = pthGC(1,:);
pTheta_CCW = repmat(pTheta_CCW, nTrialPerCondition, 1);
pTheta_CW = repmat(pTheta_CW, nTrialPerCondition, 1);
pTheta_C_new = zeros(nTrialPerCondition, nth);

%%% SAMPLING version
discriminateSim = NaN(nTrialPerCondition, length(stdSensory), length(thetaStim));
estimateSim = NaN(nTrialPerCondition, length(stdSensory), length(thetaStim));
biasEstimateSim = NaN(nTrialPerCondition, length(stdSensory), length(thetaStim));
indexCorrectSim = NaN(nTrialPerCondition, length(stdSensory), length(thetaStim));
for jj = 1 : length(stdSensory)
    for kk = 1 : length(thetaStim)   
        % Draw a sample and measure
        m = normrnd(thetaStim(kk), stdSensory(jj),[nTrialPerCondition 1]);

        % Discrimination (Inference from p(C|m))
        pM_Theta = normpdf(repmat(m,1,nth), repmat(th,nTrialPerCondition,1), stdSensory(jj));
        pCCW_m = pC(1) .* squeeze(sum((pTheta_CCW .* pM_Theta),2)); 
        pCW_m = pC(2) .* squeeze(sum((pTheta_CW .* pM_Theta),2)); 
        discriminateSim(pCCW_m >= pCW_m, jj, kk) = -1;
        discriminateSim(pCCW_m < pCW_m, jj, kk) = 1;
        tempDiscriminate = squeeze(discriminateSim(:,jj,kk));

        % Correct feedback
        correctFeedback = sign(thetaStim(kk)) * ones(size(m));
        if thetaStim(kk) < 0
            pTheta_C_new = pTheta_CCW;
        elseif thetaStim(kk) > 0
            pTheta_C_new = pTheta_CW;
        else
            inCCW = round(rand(1,nTrialPerCondition));
            pTheta_C_new(inCCW==1, :) = squeeze(pTheta_CCW(inCCW==1,:));
            pTheta_C_new(inCCW==0, :) = squeeze(pTheta_CW(inCCW==0,:));
            correctFeedback(inCCW==1) = -1;
            correctFeedback(inCCW==0) = 1;
        end                   

        % Resample the correct measurement
        Mr = m;
        indCorrect = tempDiscriminate == correctFeedback;
        mCorrect = m(indCorrect);
        if ~isempty(mCorrect)
            pResample = normpdf(repmat(th, length(mCorrect), 1), ...
                                repmat(mCorrect, 1, nth), sqrt(stdMemory^2));            
            if thetaStim(kk) < 0
                pResample(:, th > 0) = 0;
            elseif thetaStim(kk) > 0
                pResample(:,th < 0) = 0;
            else
                inCCW = inCCW(indCorrect);
                pResample(inCCW==1, th > 0) = 0;
                pResample(inCCW==0, th < 0) = 0;
            end

            mResampled = NaN(size(mCorrect));
            for tt = 1 : length(mCorrect)
                mResampled(tt) = randpdf(pResample(tt, :), th, [1 1]);
            end
            Mr(indCorrect) = mResampled;
        end

        % Estimation
        pMr_Theta = normpdf(repmat(Mr,1,nth), repmat(th,nTrialPerCondition,1), sqrt(stdSensory(jj)^2 + stdMemory^2));
        pTheta_Mr_c = pTheta_C_new .* pMr_Theta;
        pMm_c = trapz(th, pTheta_Mr_c, 2);
        pTheta_Mr_Norm = pTheta_Mr_c ./ repmat(pMm_c, 1, nth);
        thetaEstimate = trapz(th, repmat(th, nTrialPerCondition, 1) .* pTheta_Mr_Norm, 2);
        thetaResponse = normrnd(thetaEstimate, stdMotor);
        if includeIncongruentTrials == 0
            indexExcludeIncongruent = sign(thetaResponse) ~= sign(thetaEstimate);
            thetaResponse(indexExcludeIncongruent) = NaN;
            discriminateSim(indexExcludeIncongruent, jj, kk) = NaN;
        end
        estimateSim(:, jj, kk) = thetaResponse;
        biasEstimateSim(:, jj, kk) = thetaStim(kk) - thetaResponse;
        indexCorrectSim(:, jj, kk) = (tempDiscriminate == correctFeedback);
    end
end


%% MARGINALIZED version
if marginalizedVersion == 1
    estimateModelCorrect.Xval = cell(length(stdSensory));
    estimateModelCorrect.Yval = cell(length(stdSensory));  
    estimateTheoryCCW_Correct = NaN(length(stdSensory), length(thetaStim));
    estimateTheoryCW_Correct = NaN(length(stdSensory), length(thetaStim));
    estimateModelIncorrect.Xval = cell(length(stdSensory));
    estimateModelIncorrect.Yval = cell(length(stdSensory));  
    estimateTheoryCCW_Incorrect = NaN(length(stdSensory), length(thetaStim));
    estimateTheoryCW_Incorrect = NaN(length(stdSensory), length(thetaStim));
    

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

        nmr = nm;
        mr = m;
        
        M = repmat(m',1,nth);
        MM_m = repmat(mm',1,nm);
        MM_th = repmat(mm',1,nth); 
        MM_ths = repmat(mm',1,length(thetaStim));
        MR_mm = repmat(mr', 1, nmm);
        MR_th = repmat(mr', 1, nth);
        THm = repmat(th, nm, 1); 
        THmm = repmat(th, nmm, 1);
        THmr = repmat(th, nmr, 1);
        THSmm = repmat(thetaStim, nmm, 1);

        %% Correct trials
        % Generative (forward)
        % orientation noise p(m|th)
        pmGth = exp(-((M-THm).^2)./(2*stdSensory(kk)^2));
        pmGth = pmGth./(repmat(sum(pmGth,1),nm,1)); 

        % Inference
        % 1: categorical judgment
        PCGm = (pthGC * pmGth') .* repmat(pC,1,nm);
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

        % 2: estimation
        % Likelihood function the same as correct decision p(mm|th) = N(th, sm^2 + smm^2)           
        pmrGth = exp(-((MR_th-THmr).^2)./(2*(stdSensory(kk)^2 + stdMemory^2)));
        pmrGth = pmrGth./(repmat(sum(pmrGth,1),nmr,1)); 
        pthGmrChcw = (pmrGth.*repmat(pthGC(1,:),nmr,1))';
        pthGmrChcw = pthGmrChcw./repmat(sum(pthGmrChcw,1),nth,1);
        pthGmrChcw(isnan(pthGmrChcw)) = 0;

        pthGmrChccw = (pmrGth.*repmat(pthGC(2,:),nmr,1))';
        pthGmrChccw = pthGmrChccw./repmat(sum(pthGmrChccw,1),nth,1);
        pthGmrChccw(isnan(pthGmrChccw)) = 0;

        EthChcw = th * pthGmrChcw;
        EthChccw = th * pthGmrChccw;
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

        % Resample m until we have a sample that is consistent with feedback
        % p(mr|m, theta, Chat)
        MR_m = repmat(mr', 1, nm);
        pmrGmth = exp(-((MR_m-repmat(m, nmr, 1)).^2)./(2*(stdMemory^2))); 


        pmrGmthChcw = pmrGmth;
        pmrGmthChcw(mr < 0, :) = 0;
        % put the tail with all 0 to 1 (deal with small memory noise)
        indZero = sum(pmrGmthChcw, 1) == 0;
        pmrGmthChcw(mr > 0, indZero) = 1;
        pmrGmthChcw = pmrGmthChcw./(repmat(sum(pmrGmthChcw,1),nmr,1));
        pmrGmthChcw(mr > 0, indZero) = 1e-50;

        pmrGmthChccw = pmrGmth;
        pmrGmthChccw(mr > 0, :) = 0;
        % put the tail with all 0 to 1 (deal with small memory noise)
        indZero = sum(pmrGmthChccw, 1) == 0;
        pmrGmthChccw(mr < 0, indZero) = 1;        
        pmrGmthChccw = pmrGmthChccw./(repmat(sum(pmrGmthChccw,1),nmr,1));
        pmrGmthChccw(mr < 0, indZero) = 1e-50;

        % Marginalize over m that lead to cw/ccw decision to compute likelihood p(mr|theta, Chat)
        pmGthChcw = pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(1,:)',1,length(thetaStim));
        pmrGthChcw = pmrGmthChcw * pmGthChcw;   
        pmrGthChcw = pmrGthChcw ./ (repmat(sum(pmrGthChcw,1),nmr,1)); 
        pmrGthChcw(isnan(pmrGthChcw)) = 0;

        pmGthChccw = pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(2,:)',1,length(thetaStim));
        pmrGthChccw = pmrGmthChccw * pmGthChccw;
        pmrGthChccw = pmrGthChccw ./ (repmat(sum(pmrGthChccw,1),nmr,1)); 
        pmrGthChccw(isnan(pmrGthChccw)) = 0;

        a = 1./gradient(EthChcw,dstep);
        b = repmat(a',1,length(thetaStim)) .* pmrGthChcw(indKeepCw, :);        
        pthhGthChcw = interp1(EthChcw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChcw = conv2(pthhGthChcw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChcw(pthhGthChcw < 0) = 0; 

        a = 1./gradient(EthChccw,dstep);


        b = repmat(a',1,length(thetaStim)) .* pmrGthChccw(indKeepCcw, :);        
        pthhGthChccw = interp1(EthChccw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChccw = conv2(pthhGthChccw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChccw(pthhGthChccw < 0) = 0; 

        pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
        pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);             

        if includeIncongruentTrials == 0
            % modify psychometric curve p(Chat|theta, Congruent) ~ p(Congruent| Chat, theta) * p(Chat|Theta)
            pCongruentGcwTh = sum(pthhGthChcw(th' >= 0, :));
            pCongruentGccwTh = sum(pthhGthChccw(th' <= 0, :));
            PChGtheta_lapse_new = PChGtheta_lapse .* [pCongruentGcwTh; pCongruentGccwTh];
            PChGtheta_lapse_new = PChGtheta_lapse_new ./ repmat(sum(PChGtheta_lapse_new, 1), 2, 1);

            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th'>= 0, :) = 0;
            pthhGthChcw(th'< 0, :) = 0;
        else
            PChGtheta_lapse_new = PChGtheta_lapse;
        end
        
        % remove incorrect trials
        pthhGthChcw(:, thetaStim < 0) = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
        pthhGthChccw(:, thetaStim > 0) = 0;


        pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            
        mthhGthChcw_correct = th * pthhGthChcw_norm;
        mthhGthChccw_correct = th * pthhGthChccw_norm;
        mthhGthChcw_correct(thetaStim < 0) = NaN;
        mthhGthChccw_correct(thetaStim > 0) = NaN;

        pthhGth_correct = pthhGthChcw.*repmat(PChGtheta_lapse(1,:),nth,1) + pthhGthChccw.*repmat(PChGtheta_lapse(2,:),nth,1);
        pthhAndth_correct = pthhGth_correct;

        estimateModelCorrect.Xval{kk} = th;
        estimateModelCorrect.Yval{kk} = pthhAndth_correct;
        estimateTheoryCCW_Correct(kk,:) = mthhGthChccw_correct;
        estimateTheoryCW_Correct(kk,:) = mthhGthChcw_correct; 
    end
end

% Analyze and plot
estimateResponseSimCCW_correct = NaN(size(estimateSim));
estimateResponseSimCCW_correct(discriminateSim == -1 & indexCorrectSim == 1) = estimateSim(discriminateSim == -1 & indexCorrectSim == 1);
estimateResponseSimCCWAve_correct = squeeze(nanmean(estimateResponseSimCCW_correct, 1));
estimateResponseSimCCWStd_correct = squeeze(nanstd(estimateResponseSimCCW_correct, 1))./(sqrt(squeeze(nansum(~isnan(estimateResponseSimCCW_correct),1))));
estimateResponseSimCW_correct = NaN(size(estimateSim));
estimateResponseSimCW_correct(discriminateSim == 1 & indexCorrectSim == 1) = estimateSim(discriminateSim == 1 & indexCorrectSim == 1);
estimateResponseSimCWAve_correct = squeeze(nanmean(estimateResponseSimCW_correct, 1));
estimateResponseSimCWStd_correct = squeeze(nanstd(estimateResponseSimCW_correct, 1))./(sqrt(squeeze(nansum(~isnan(estimateResponseSimCW_correct),1))));

estimateResponseSimCCW_incorrect = NaN(size(estimateSim));
estimateResponseSimCCW_incorrect(discriminateSim == -1 & indexCorrectSim == 0) = estimateSim(discriminateSim == -1 & indexCorrectSim == 0);
estimateResponseSimCCWAve_incorrect = squeeze(nanmean(estimateResponseSimCCW_incorrect, 1));
estimateResponseSimCCWStd_incorrect = squeeze(nanstd(estimateResponseSimCCW_incorrect, 1))./(sqrt(squeeze(nansum(~isnan(estimateResponseSimCCW_incorrect),1))));
estimateResponseSimCW_incorrect = NaN(size(estimateSim));
estimateResponseSimCW_incorrect(discriminateSim == 1 & indexCorrectSim == 0) = estimateSim(discriminateSim == 1 & indexCorrectSim == 0);
estimateResponseSimCWAve_incorrect = squeeze(nanmean(estimateResponseSimCW_incorrect, 1));
estimateResponseSimCWStd_incorrect = squeeze(nanstd(estimateResponseSimCW_incorrect, 1))./(sqrt(squeeze(nansum(~isnan(estimateResponseSimCW_incorrect),1))));

hAverage = figure;
lineWidth = 4;
figPos = [0.3, 0.2, 0.5, 0.7];
set(hAverage,'Units','normalized','Position',figPos)
hold on
legendName = cell(1,length(stdSensory));
colorName = {'g', 'r', 'b', 'cyan', 'magenta', 'y'};
maxPlot = max(thetaStim);
hLegend = NaN(1, length(stdSensory));
biasMean_correct = NaN(length(stdSensory), rangeCollapse);
biasStd_correct = NaN(length(stdSensory), rangeCollapse);
biasTheory_correct = NaN(length(stdSensory), rangeCollapse);
biasMean_incorrect = NaN(length(stdSensory), rangeCollapse);
biasStd_incorrect = NaN(length(stdSensory), rangeCollapse);
biasTheory_incorrect = NaN(length(stdSensory), rangeCollapse);

for ii = 1 : length(stdSensory)
    % Correct trials
    set(gca,'FontSize',fontSize)
    hShade = shadedErrorBar(thetaStim, estimateResponseSimCCWAve_correct(ii,:), estimateResponseSimCCWStd_correct(ii,:),... 
                         {'Color', colorName{ii}, 'LineWidth', lineWidth},1,0,0);        
    hLegend(ii) = hShade.patch;
    hold on
    shadedErrorBar(thetaStim, estimateResponseSimCWAve_correct(ii, :), estimateResponseSimCWStd_correct(ii, :),... 
                         {'Color', colorName{ii}, 'LineWidth', lineWidth},1,0,0);        
    if  marginalizedVersion
        hLegend(ii) = plot(thetaStim, estimateTheoryCCW_Correct(ii,:), 'Color', colorName{ii}, 'LineWidth', lineWidth);
        plot(thetaStim, estimateTheoryCW_Correct(ii,:), 'Color', colorName{ii}, 'LineWidth', lineWidth);
    end
    legendName{ii} = ['Noise level = ' num2str(stdSensory(ii))];
    xlabel('True angle (degree)')
    ylabel('Angle estimate (degree)')  
    plot([-maxPlot maxPlot], [-maxPlot maxPlot], 'k--', 'LineWidth', 2)
    xlim([-maxPlot maxPlot])
    ylim([-maxPlot maxPlot])
    text(6, -15, '+/-1SEM', 'FontSize', 15)

    % Get the bias
    tempEstimateCW = squeeze(estimateResponseSimCW_correct(:,ii,:));
    tempEstimateCCW = squeeze(estimateResponseSimCCW_correct(:,ii,:));
    tempBias1 = tempEstimateCW(:,rangeCollapse:end) - repmat(thetaStim(rangeCollapse:end),size(tempEstimateCW,1),1);
    tempBias2 = tempEstimateCCW(:,1:rangeCollapse) - repmat(thetaStim(1:rangeCollapse),size(tempEstimateCCW,1),1);
    tempBias = [tempBias1;-tempBias2(:,sort(1:size(tempBias2,2),'descend'))];
    biasMean_correct(ii,:) = squeeze(nanmean(tempBias,1));
    biasStd_correct(ii,:) = squeeze(nanstd(tempBias,1))./sqrt(nansum(~isnan(tempBias),1)); 

    if marginalizedVersion
        tempBiasTheory1 = estimateTheoryCW_Correct(ii,rangeCollapse:end) - thetaStim(rangeCollapse:end);
        tempBiasTheory2 = estimateTheoryCCW_Correct(ii,1:rangeCollapse) - thetaStim(1:rangeCollapse);
        tempBiasTheory = [tempBiasTheory1;-tempBiasTheory2(:,sort(1:size(tempBiasTheory2,2),'descend'))];  
        biasTheory_correct(ii,:) = squeeze(nanmean(tempBiasTheory,1));
    end
end
legend(hLegend, legendName, 'Location', 'NorthWest')

%% Plot the estimate distribution
% Correct trials
figure
counter = 1;
binCenter = linspace(-40, 40, 41);
tempEstimateModelX = estimateModelCorrect.Xval{1};

for kk = 1 : length(estimateModelCorrect.Yval)
    estimateModelY = estimateModelCorrect.Yval{kk};
    estimateModelY = estimateModelY(:, thetaStim >=0);
    estimateSimulation = estimateSim(:, kk, thetaStim >=0);
    indexCorrectSimulation = squeeze(indexCorrectSim(:, kk, thetaStim >=0));
    for jj = 1 : size(estimateModelY, 2)            
        tempEstimateModelY = estimateModelY(:, jj); 
        tempEstimateModelY = tempEstimateModelY ./ trapz(tempEstimateModelX, tempEstimateModelY);
        tempestimateDataX = estimateSimulation(indexCorrectSimulation(:, jj) == 1, 1, jj);
        
        subplot(length(estimateModelCorrect.Yval), size(estimateModelY, 2), counter)
        hold on
        binCount = hist(tempestimateDataX, binCenter);
        binCount = binCount / (sum(binCount) * diff(binCenter(1:2)));
        bar(binCenter, binCount, 1);
 
        plot(tempEstimateModelX, tempEstimateModelY,'r','LineWidth',3)
        plot([0 0], [0 0.1], 'k--')
        xlim([-50 50])
        set(gca,'YTickLabel',[],'YTick',[])
        if kk == 1
            title(num2str(thetaStim(jj+7)))
        end
        box off
        
        counter = counter + 1;
    end
end
tightfig
