%%%%%%%% The full marginalized Bayesian model of conditioned perception %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
% 
clear; 

subjectID = {'skb'};
experimentNumber = 1;
experimentType = 'MainExperiment';
experiment = 'DecisionGiven'; % DecisionGiven   Original
session = 1;
dataAll = [];
fontSize = 20;
if strcmp(subjectID{1}, 'average')
    subjectIDD = '';
else
    subjectIDD = subjectID{1};
end
if strcmp(experiment, 'Original')
    expNumber = 2;
elseif strcmp(experiment, 'DecisionGiven')
    expNumber = 3;
elseif strcmp(experiment, 'ControlReplication')
    expNumber = 1;
end

includeIncongruentTrials = ''; % empty: no incongruent trials; 

flagSC = 1; % 1: Self-conditioned                             
lineWidth = 3;
bootstrapParams = 0;
if strcmp(experiment, 'ControlReplication')
    experimentNum = 'Exp1';
    flagDecisionGiven = 0;
elseif strcmp(experiment, 'Original')
    experimentNum = 'Exp2';
    flagDecisionGiven = 0;
elseif strcmp(experiment, 'DecisionGiven')
    experimentNum = 'Exp3';   
    flagDecisionGiven = 1;   
end
if bootstrapParams
    fileName = ['Bootstrap_' subjectIDD '_' experimentNum '.txt'];
    fileID = fopen(fileName);
    paramAll = textscan(fileID,'%f %f %f %f %f %f %f %f %f','CommentStyle','//');
    fclose(fileID);
    paramsFit = cell2mat(paramAll(:, 2:end));
else
    paramsFit = [4.3172    4.5344    8.4083      0.0000    33.9063  14.5556    0.5305    2.0456];
end
stdSensory = paramsFit(:, 1:3);
priorRange = paramsFit(:, 5);
smoothFactor = paramsFit(:, 7);
stdMemory = paramsFit(:, 6);
if stdMemory == 0
    stdMemory = 0.0001;
end
lapseRate = paramsFit(:, 4);
stdMotor = paramsFit(:, 8);
optimizationAlgorithm = 1;

for ii = 1 : length(experimentNumber)
    if strcmp(experiment, 'PilotData')
        dataFullPath = fullfile('Data', subjectID{1}, experimentType, experiment, ['Session' num2str(session)], ...
                                    ['ConditionDecisionMaking-' num2str(experimentNumber(ii))]);
    elseif strcmp(experimentType, 'MotorNoise')
        dataFullPath = fullfile('Data', subjectID{1}, experimentType, ['Session' num2str(session)], ...
                                    ['MotorNoise-' num2str(experimentNumber(ii))]);
    elseif strcmp(experimentType, 'PerceptNoise')
        dataFullPath = fullfile('Data', subjectID{1}, experimentType, ['Session' num2str(session)], ...
                                    ['PerceptNoise-' num2str(experimentNumber(ii))]);                                
    else
        dataFullPath = fullfile('Data', subjectID{1}, experimentType,[experiment num2str(session)], ...
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

% Extract data
stdNoiseExp = unique(dataAll(:,4));
angleDiff = (unique(dataAll(:,1)))';
dataZeroDiff = dataAll(dataAll(:,1) == 0,:);
nBlocks = ceil(params.nTrialPerCondition/length(params.barAngleReference));
angleEstimate = dataAll(:, nBlocks+6:2*nBlocks+5);

% Average across conditions
angleBar2 = dataAll(:,5) - dataAll(:,1);
indexCorrect = find(abs(angleBar2-angleEstimate) > 90);
angleEstimate(indexCorrect) = angleEstimate(indexCorrect) + sign(angleBar2(indexCorrect)-angleEstimate(indexCorrect))*180;    
angleDiffEst = dataAll(:, 5) - angleEstimate;
biasAll = angleBar2 - angleEstimate;
if strcmp(experiment, 'DecisionGiven')
    binaryDecisionAll = dataAll(:,8);
    indicatorConsistent = sign(angleDiffEst)==dataAll(:,8);
elseif strcmp(experimentType, 'MotorNoise')
    indicatorConsistent = ones(length(biasAll));
    binaryDecisionAll = dataAll(:,6);    
else
    binaryDecisionAll = dataAll(:,6);
    indicatorConsistent = sign(angleDiffEst)==dataAll(:,6);    
end    
if ~isempty(includeIncongruentTrials)
    indicatorConsistent = ones(length(biasAll), 1);
end

angleDiffEstimate = cell(length(stdNoiseExp), length(angleDiff));
biasEstimate = cell(length(stdNoiseExp), length(angleDiff));
binaryDecisionSelf = cell(length(stdNoiseExp), length(angleDiff));
binaryDecisionGiven = binaryDecisionSelf;
for ii = 1 : length(stdNoiseExp)
    for jj = 1 : length(angleDiff)
        tempBias = biasAll(dataAll(:,4)==stdNoiseExp(ii) ...
                           & dataAll(:,1)==angleDiff(jj), :);
        tempIndicatorConsistent = indicatorConsistent(dataAll(:,4)==stdNoiseExp(ii) ...
                                                    & dataAll(:,1)==angleDiff(jj), :);                                
        biasEstimate{ii,jj} = tempBias(:);
        tempAngleEllispe1 = dataAll(:,5);
        tempAngleEllispe1 = tempAngleEllispe1(dataAll(:,4)==stdNoiseExp(ii) ...
                                            & dataAll(:,1)==angleDiff(jj), :);
        tempAngleEstimate = angleEstimate(dataAll(:,4)==stdNoiseExp(ii) ...
                                            & dataAll(:,1)==angleDiff(jj), :); 
        tempAngleEstimate(~tempIndicatorConsistent) = NaN;                                
        angleDiffEstimate{ii,jj} = tempAngleEllispe1-tempAngleEstimate;
        if strcmp(experiment, 'DecisionGiven')
            binaryDecisionSelf{ii,jj} = sign(angleDiffEstimate{ii,jj});
            binaryDecisionGiven{ii,jj} = binaryDecisionAll(dataAll(:,4)==stdNoiseExp(ii) ...
                                            & dataAll(:,1)==angleDiff(jj), :);     
            binaryDecision = binaryDecisionGiven;                            
        else
            binaryDecisionSelf{ii,jj} = binaryDecisionAll(dataAll(:,4)==stdNoiseExp(ii) ...
                                            & dataAll(:,1)==angleDiff(jj), :);
            binaryDecision = binaryDecisionSelf;                            
        end
    end
end

fractionCW = NaN(length(stdNoiseExp),length(params.barAngleDiff));
nTrialsPerCondition = NaN(length(stdNoiseExp),length(params.barAngleDiff));
if strcmp(experiment, 'DecisionGiven')
    binaryDecisionAll = dataAll(:, 8);
else
    binaryDecisionAll = dataAll(:, 6);
end    
for ii = 1 : length(stdNoiseExp)
    for jj = 1 : length(fractionCW)
        indexAll = dataAll(:,1) == params.barAngleDiff(jj) ...
                    & dataAll(:,4) == stdNoiseExp(ii)...
                    & indicatorConsistent;
        indexCW = dataAll(:,1) == params.barAngleDiff(jj) ...
                    & dataAll(:,4) == stdNoiseExp(ii)...
                    & indicatorConsistent ...
                    & binaryDecisionAll == 1;
        nTrialsPerCondition(ii,jj) = nansum((indexAll));
        fractionCW(ii,jj) = nansum(indexCW)/nansum(indexAll);
    end
end


%% Modeling
thetaStim = -21:0.5:21;
pC = [0.5, 0.5]'; % [cw ccw]
pthcw = priorRange;
pthccw = -priorRange;

dstep = 0.1;

rangeth = [-42 42];
th = rangeth(1):dstep:rangeth(2);
nth = length(th);

pthGC = zeros(2,nth);
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

pCw = NaN(length(stdSensory), length(thetaStim));
estimateModel.Xval = cell(length(stdSensory));
estimateModel.Yval = cell(length(stdSensory));  
estimateTheoryCCW = NaN(length(stdSensory), length(thetaStim));
estimateTheoryCW = NaN(length(stdSensory), length(thetaStim));

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
    % orientation noise
    pmGth = exp(-((M-THm).^2)./(2*stdSensory(kk)^2));
    pmGth = pmGth./(repmat(sum(pmGth,1),nm,1)); 

    %% Inference
    % 1: categorical judgment
    if ~flagDecisionGiven
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
    else
        PChGtheta_lapse = NaN(2, length(thetaStim));
        PChGtheta_lapse(1, thetaStim > 0) = 1;
        PChGtheta_lapse(1, thetaStim < 0) = 0;
        PChGtheta_lapse(1, thetaStim == 0) = 0.5;
        PChGtheta_lapse(2, :) = 1 - PChGtheta_lapse(1, :);
    end
        
    % 2: estimation
    if ~flagDecisionGiven
        % memory noise
        pmmGm = exp(-((MM_m-repmat(m, nmm, 1)).^2)./(2*stdMemory^2)); 
        pmmGm = pmmGm./(repmat(sum(pmmGm,1),nmm,1));   

        % attention marginalization: compute distribution only over those ms that lead to cw decision!
        pmmGthChcw = pmmGm * (pmGth .* repmat(PChGm(1,:)',1,nth));
        pmmGthChcw = pmmGthChcw ./ repmat(sum(pmmGthChcw,1),nmm,1);
        pmmGthChccw = pmmGm * (pmGth .* repmat(PChGm(2,:)',1,nth));
        pmmGthChccw = pmmGthChccw ./ repmat(sum(pmmGthChccw,1),nmm,1);

        pthGmmChcw = (pmmGthChcw.*repmat(pthGC(1,:),nmm,1))';
        pthGmmChcw = pthGmmChcw./repmat(sum(pthGmmChcw,1),nth,1);
        pthGmmChcw(isnan(pthGmmChcw)) = 0;      

        pthGmmChccw = (pmmGthChccw.*repmat(pthGC(2,:),nmm,1))';
        pthGmmChccw = pthGmmChccw./repmat(sum(pthGmmChccw,1),nth,1);
        pthGmmChccw(isnan(pthGmmChccw)) = 0;                        
    else
        pmmGth = exp(-((MM_th-THmm).^2)./(2*(stdSensory(kk)^2 + stdMemory^2))); % p(mm|th) = N(th, sm^2 + smm^2)
        pmmGth = pmmGth./(repmat(sum(pmmGth,1),nmm,1)); 

        pthGmmChcw = (pmmGth.*repmat(pthGC(1,:),nmm,1))';
        pthGmmChcw = pthGmmChcw./repmat(sum(pthGmmChcw,1),nth,1);
        pthGmmChcw(isnan(pthGmmChcw)) = 0;

        pthGmmChccw = (pmmGth.*repmat(pthGC(2,:),nmm,1))';
        pthGmmChccw = pthGmmChccw./repmat(sum(pthGmmChccw,1),nth,1);
        pthGmmChccw(isnan(pthGmmChccw)) = 0;            
    end


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
        % attention marginalization: compute distribution only over those ms that lead to cw decision!
        pmmGthChcw = pmmGthChcw(:, ismember(th, thetaStim));   
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
        pmmGthChccw = pmmGthChccw(:, ismember(th, thetaStim));        
        b = repmat(a',1,length(thetaStim)) .* pmmGthChccw(indKeepCcw, :);        
    else
        b = repmat(a',1,length(thetaStim)) .* pmmGth(indKeepCcw, ismember(th, thetaStim));
    end  
    pthhGthChccw = interp1(EthChccw,b,th,'linear','extrap');
    % add motor noise
    pthhGthChccw = conv2(pthhGthChccw,pdf('norm',th,0,stdMotor)','same');
    pthhGthChccw(pthhGthChccw < 0) = 0; 


    if isempty(includeIncongruentTrials)
        pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
        pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            

        % modify psychometric curve p(Chat|theta, Congruent) ~ p(Congruent| Chat, theta) * p(Chat|Theta)
        pCongruentGcwTh = sum(pthhGthChcw(th' >= 0, :));
        pCongruentGccwTh = sum(pthhGthChccw(th' <= 0, :));
        PChGtheta_lapse = PChGtheta_lapse .* [pCongruentGcwTh; pCongruentGccwTh];
        PChGtheta_lapse = PChGtheta_lapse ./ repmat(sum(PChGtheta_lapse, 1), 2, 1);

        % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
        pthhGthChccw(th'>= 0, :) = 0;
        pthhGthChcw(th'< 0, :) = 0;
    end 
    pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
    pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);    
    mthhGthChcw = th * pthhGthChcw;
    mthhGthChccw = th * pthhGthChccw;

    pthhGth = pthhGthChcw.*repmat(PChGtheta_lapse(1,:),nth,1) + pthhGthChccw.*repmat(PChGtheta_lapse(2,:),nth,1);

    pCw(kk, :) = PChGtheta_lapse(1,:);
    estimateModel.Xval{kk} = th;
    estimateModel.Yval{kk} = pthhGth;
    estimateTheoryCCW(kk,:) = mthhGthChccw;
    estimateTheoryCW(kk,:) = mthhGthChcw;
end    
    
    
% Analyze
rangeCollapseTheory = round(length(thetaStim)/2);
biasTheory = NaN(length(stdSensory), rangeCollapseTheory);
for ii = 1 : length(stdSensory)
    tempBiasTheory1 = estimateTheoryCW(ii, rangeCollapseTheory:end) - thetaStim(rangeCollapseTheory:end);
    tempBiasTheory2 = estimateTheoryCCW(ii, 1:rangeCollapseTheory) - thetaStim(1:rangeCollapseTheory);
    tempBiasTheory = [tempBiasTheory1;-tempBiasTheory2(:,sort(1:size(tempBiasTheory2,2),'descend'))];  
    biasTheory(ii,:) = squeeze(nanmean(tempBiasTheory,1));
end

%% Plot the psychometric curves
colorName = {'g', 'r', 'b', 'magenta', 'y'};
figure
hold on
set(gca, 'FontSize', fontSize)
for ii = 1 : size(pCw, 1)
    plot(thetaStim, pCw(ii, :), 'Color', colorName{ii},  'LineWidth', 5)
    plot(angleDiff, fractionCW(ii, :), 'o', 'MarkerSize', 15, 'Color', colorName{ii}, 'MarkerFaceColor', colorName{ii})
end
xlabel('Stimulus orientation (deg)')
ylabel('Fraction "cw"')
xlim([min(angleDiff) max(angleDiff)])

%% Plot the average
hAverage = figure;
figPos = [0.3, 0.2, 0.5, 0.7];
set(hAverage,'Units','normalized','Position',figPos)
hold on
rangeCollapse = round(length(angleDiff)/2);
legendName = cell(1,length(stdSensory));
maxYplot = max(params.barAngleDiff)+10;
maxXplot = max(params.barAngleDiff);
indexColor = 1;
hLegend = NaN(1, length(stdSensory));
boxcarLength = 1;

% Plot true angle
biasMean = NaN(length(stdSensory), rangeCollapse);
biasStd = NaN(length(stdSensory), rangeCollapse);
for ii = 1 : length(stdSensory)
    hold on
    set(gca,'FontSize',fontSize)
    tempAngleDiffEstimate = [angleDiffEstimate{ii,:}];
    tempBinaryDecision = [binaryDecision{ii,:}];
    tempAngleDiffEstimateCW = NaN(size(tempAngleDiffEstimate));
    tempAngleDiffEstimateCW(tempBinaryDecision == 1) = tempAngleDiffEstimate(tempBinaryDecision == 1);
    tempAngleDiffEstimateCWAve = squeeze(nanmean(tempAngleDiffEstimateCW,1));
    tempAngleDiffEstimateCWAve(1:rangeCollapse-1) = NaN;
    tempAngleDiffEstimateCWStd = squeeze(nanstd(tempAngleDiffEstimateCW,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCW),1));
    tempBias1 = tempAngleDiffEstimateCW(:,rangeCollapse:end) - repmat(angleDiff(rangeCollapse:end),size(tempAngleDiffEstimateCW,1),1);
    tempBias1 = tempBias1 - repmat(biasMotor(rangeCollapse:end), size(tempBias1,1),1);
    hShade = shadedErrorBar(angleDiff, tempAngleDiffEstimateCWAve, tempAngleDiffEstimateCWStd,... 
                         {'Color', colorName{ii}, 'LineWidth', lineWidth}, 1, 0, 0);        

    tempAngleDiffEstimateCCW = NaN(size(tempAngleDiffEstimate));
    tempAngleDiffEstimateCCW(tempBinaryDecision == -1) = tempAngleDiffEstimate(tempBinaryDecision == -1);
    tempAngleDiffEstimateCCWAve = squeeze(nanmean(tempAngleDiffEstimateCCW,1));
    tempAngleDiffEstimateCCWAve(rangeCollapse+1:end) = NaN;
    tempAngleDiffEstimateCCWStd = squeeze(nanstd(tempAngleDiffEstimateCCW,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCCW),1));
    tempBias2 = tempAngleDiffEstimateCCW(:,1:rangeCollapse) - repmat(angleDiff(1:rangeCollapse),size(tempAngleDiffEstimateCCW,1),1);
    tempBias2 = tempBias2 - repmat(biasMotor(1:rangeCollapse), size(tempBias2,1),1);
    tempBias = [tempBias1;-tempBias2(:,sort(1:size(tempBias2,2),'descend'))];
    tempBiasMean = nanmean(tempBias,1);
    bias = [ones(1,boxcarLength-1)*tempBiasMean(1) tempBiasMean ones(1,boxcarLength-1)*tempBiasMean(end)];
    smoothBias = conv(bias, ones(1,boxcarLength)/boxcarLength, 'full');
    biasMean(ii,:) = smoothBias(boxcarLength:boxcarLength+rangeCollapse-1);
    biasStd(ii,:) = squeeze(nanstd(tempBias,1))./sqrt(nansum(~isnan(tempBias),1));        
    shadedErrorBar(angleDiff, tempAngleDiffEstimateCCWAve, tempAngleDiffEstimateCCWStd,... 
                         {'Color', colorName{ii}, 'LineWidth', lineWidth}, 1, 0, 0); 
    hLegend(ii) = plot(thetaStim, estimateTheoryCCW(ii,:), 'Color', colorName{ii}, 'LineWidth', lineWidth);
    plot(thetaStim, estimateTheoryCW(ii,:), 'Color', colorName{ii}, 'LineWidth', lineWidth);

    legendName{ii} = ['Stimulus noise = ' num2str(stdSensory(ii))];
    xlabel('\theta_{true} (degree)')
    ylabel('\theta_{estimate} (degree)')  
    title(['Subject ' upper(subjectID{1}) ]);
    plot([angleDiff(1) angleDiff(end)],  [angleDiff(1) angleDiff(end)], '--k')
    xlim([-maxXplot maxXplot])
    ylim([-maxYplot maxYplot])
    box on
    text(6, -15, '+/-1SEM', 'FontSize', 15)
end
legend(hLegend, legendName, 'Location', 'NorthWest')

% Plot the bias
hBias = figure;
hold on
set(gca,'FontSize',fontSize)
for ii = 1 : length(stdSensory)
    hold on
    if strcmp(subjectID{1}, 'average')
        hShade = shadedErrorBar(angleDiff(rangeCollapse:end), biasMean(ii,:), biasStd(ii,:),... 
                             {'Color', colorName{ii}, 'LineWidth', lineWidth},1,0,0); 
    end
    hLegend(ii) = plot(thetaStim(rangeCollapseTheory:end), biasTheory(ii,:), ... 
                    'Color', colorName{ii}, 'LineWidth', lineWidth);
    legendName{ii} = ['noise = ' num2str(stdSensory(ii))];
    xlabel('Absolute true angle (degree)')
    ylabel('Bias (degree)')  
%     title(['Subject ' upper(subjectID) ]);
    xlim([0 maxXplot])
    ylim([-10 10])
    if strcmp(experiment, 'ControlReplication')
        ylim([-10 15])
    else
        ylim([-10 15])
        set(gca, 'YTick', [-10:5:15]) 
    end
end
plot([0 21], [0 0], 'k--', 'LineWidth', lineWidth)
% legend(hLegend, legendName) 

%% Extract the collapsed data
[binaryDecision, percentCW, nTrialsPerCondition, estimateData, angleDiff, stdMotor] = dataForFitting(subjectID{1}, expNumber, includeIncongruentTrials);


%% Plot the estimate distribution
figure
counter = 1;
binCenter = linspace(-40, 40, 29);
tempEstimateModelX = estimateModel.Xval{1};
angleDiffCollapse = angleDiff(angleDiff >= 0);
for kk = 1 : length(estimateModel.Yval)
    estimateModelY = estimateModel.Yval{kk};
    estimateModelY = estimateModelY(:, ismember(thetaStim, angleDiffCollapse));
    for jj = 1 : size(estimateModelY, 2)            
        tempEstimateModelY = estimateModelY(:, jj); 
        tempEstimateModelY = tempEstimateModelY ./ trapz(tempEstimateModelX, tempEstimateModelY);
        tempestimateDataX = estimateData{kk,jj};
        
        subplot(length(estimateModel.Yval), size(estimateModelY, 2), counter)
        hold on
        binCount = hist(tempestimateDataX, binCenter);
        binCount = binCount / (sum(binCount) * diff(binCenter(1:2)));
        bar(binCenter, binCount, 1);
 
        plot(tempEstimateModelX, tempEstimateModelY,'r','LineWidth',3)
        plot([0 0], [0 0.095], 'k--')
        xlim([-42 42])
        ylim([0 0.095])
        set(gca,'YTickLabel',[],'YTick',[])
        if kk == 1
            title(num2str(angleDiffCollapse(jj)))
        end
        box off
        
        counter = counter + 1;
    end
end
