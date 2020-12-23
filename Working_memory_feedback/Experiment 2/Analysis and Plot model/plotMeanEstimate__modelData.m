%%%%%%%%%%%%%%%%%%%%%%% Compute the mean estimate of the models and data %%%%%%%%%%%%%%%%%%%%%%

%% Compute the subject LLH
subjectIDAll = {'ll', 'pw', 'eh', 'bh', 'ln', 'at', 'dh'}; % 
estimateCW_Prior = NaN(length(subjectIDAll), 2);
estimateCCW_Prior = NaN(length(subjectIDAll), 2);
estimateCW_Variance = NaN(length(subjectIDAll), 2);
estimateCCW_Variance = NaN(length(subjectIDAll), 2);
estimateCW_FlipDecisionAddMem = NaN(length(subjectIDAll), 2);
estimateCCW_FlipDecisionAddMem = NaN(length(subjectIDAll), 2);
estimateCW_FlipDecision = NaN(length(subjectIDAll), 2);
estimateCCW_FlipDecision = NaN(length(subjectIDAll), 2);
estimateCW_Resample = NaN(length(subjectIDAll), 2);
estimateCCW_Resample = NaN(length(subjectIDAll), 2);
estimateCW_Surprise = NaN(length(subjectIDAll), 2);
estimateCCW_Surprise = NaN(length(subjectIDAll), 2);

estimateCW_Data = NaN(length(subjectIDAll), 2);
estimateCCW_Data = NaN(length(subjectIDAll), 2);
ciCW_Data = NaN(length(subjectIDAll), 2, 2);
ciCCW_Data = NaN(length(subjectIDAll), 2, 2);

% % No resample for correct trials                
% paramsAllSubject = [4.3425    6.2248           0.0000     19.1377    -9.6045   3.5283    2.0902    0.9927    0.6813;
%                     8.7205    8.8878           0.0000     33.0312   -21.2586   1.2076    1.8928    0.9215    0.4348;
%                     8.3099    8.7268           0.0000     13.9033   -12.6251   0.6127    2.7094    0.9468    0.5099;
%                     6.3379    8.4823           0.0000     19.3091   -12.7690   0.9699    2.6041    0.5558    0.4201;
%                     6.4379    9.9076           0.0000     32.5355   -17.6389   3.0964    1.5830    0.2067    0.4086;
%                     9.7284   15.1847           0.0000     56.3034   -42.2707   0.4378    4.0136    0.9881    0.5523;
%                     7.8510    9.8641           0.0000     22.8922   -18.0641   0.8543    3.9069    0.9646    0.4949];

% No resample for correct trials (fit relapse)               
paramsAllSubject = [2.5841    4.6918           0.0327     19.1233    -9.9073   5.5023    2.0902    0.9984    0.7464;
                    8.6334    8.8191           0.0043     33.0647   -21.3477   0.8755    1.8928    0.9127    0.4313;
                    7.4646    7.8156           0.0241     13.6883   -12.2943   1.0845    2.7094    0.8384    0.5343;
                    5.8649    8.0372           0.0075     19.2718   -12.7879   0.7122    2.6041    0.4691    0.4097;
                    6.2094    9.8275           0.0010     32.6800   -17.6770   3.9740    1.5830    0.1959    0.4021;
                    9.4508   14.9103           0.0007     58.8196   -39.8226   0.1698    4.0136    0.2880    0.5929;
                    6.8495    9.2325           0.0127     22.7921   -18.3783   1.0853    3.9069    0.8675    0.5019];

% % Resample for correct trials                
% paramsAllSubject = [5.1714    6.5149           0.0000     18.2531    -9.3607   0.9940    2.0902    0.9978    0.6635;
%                     8.5063    8.7485           0.0000     33.2603   -21.0845   0.6393    1.8928    0.9098    0.4518;
%                     8.7284    8.9988           0.0000     13.8798   -12.3746   0.6037    2.7094    0.9189    0.5136;
%                     6.1961    8.3406           0.0000     19.5003   -12.7209   0.7809    2.6041    0.4987    0.4168;
%                     5.1595    8.9406           0.0000     32.3917   -17.5501   5.2297    1.5830    0.3252    0.3550;
%                     10.0963   15.0721           0.0000     57.4518   -39.0410   0.0101    4.0136    0.4339    0.5875;
%                     7.7473    9.7957           0.0000     23.1824   -18.2993   0.6654    3.9069    0.8918    0.4908];
                
for nn = 1 : length(subjectIDAll)
    subjectID = subjectIDAll{nn};
    if strcmp(subjectID, 'average')
        experimentNumber = 2:4;
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
    %  Column 1: true angle difference (population)
    %  Column 2: bar presentation time
    %  Column 3: SOA
    %  Column 4: bar noise level
    %  Column 5: bar reference angle
    %  Column 6 to 8: data (categorical decision and estimation)
    %  Column 9: true angle difference (sample) 
    %  Column 10: decision time
    %  Column 11: estimation time    stdNoiseLevel = unique(dataAll(:,4));
    
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
    
    estimateCW_Data(nn, 1) = nanmean(angleDiffEst(dataAll(:, 4) == stdNoiseLevel(1) & binaryFeedbackAll == 1));
    estimateCCW_Data(nn, 1) = abs(nanmean(angleDiffEst(dataAll(:, 4) == stdNoiseLevel(1) & binaryFeedbackAll == -1)));
    estimateCW_Data(nn, 2) = nanmean(angleDiffEst(dataAll(:, 4) == stdNoiseLevel(2) & binaryFeedbackAll == 1));
    estimateCCW_Data(nn, 2) = abs(nanmean(angleDiffEst(dataAll(:, 4) == stdNoiseLevel(2) & binaryFeedbackAll == -1)));
    
    ciCW_Data(nn, 1, :) = bootci(1000, @nanmean, angleDiffEst(dataAll(:, 4) == stdNoiseLevel(1) & binaryFeedbackAll == 1));
    ciCCW_Data(nn, 1, :) = bootci(1000, @nanmean, abs(angleDiffEst(dataAll(:, 4) == stdNoiseLevel(1) & binaryFeedbackAll == -1)));
    ciCW_Data(nn, 2, :) = bootci(1000, @nanmean, angleDiffEst(dataAll(:, 4) == stdNoiseLevel(2) & binaryFeedbackAll == 1));
    ciCCW_Data(nn, 2, :) = bootci(1000, @nanmean, abs(angleDiffEst(dataAll(:, 4) == stdNoiseLevel(2) & binaryFeedbackAll == -1)));

%     estimateCW_Data(nn, 1) = nanmedian(angleDiffEst(dataAll(:, 4) == stdNoiseLevel(1) & binaryFeedbackAll == 1));
%     estimateCCW_Data(nn, 1) = abs(nanmedian(angleDiffEst(dataAll(:, 4) == stdNoiseLevel(1) & binaryFeedbackAll == -1)));
%     estimateCW_Data(nn, 2) = nanmedian(angleDiffEst(dataAll(:, 4) == stdNoiseLevel(2) & binaryFeedbackAll == 1));
%     estimateCCW_Data(nn, 2) = abs(nanmedian(angleDiffEst(dataAll(:, 4) == stdNoiseLevel(2) & binaryFeedbackAll == -1)));
%     ciCW_Data(nn, 1, :) = bootci(1000, @nanmedian, angleDiffEst(dataAll(:, 4) == stdNoiseLevel(1) & binaryFeedbackAll == 1));
%     ciCCW_Data(nn, 1, :) = bootci(1000, @nanmedian, abs(angleDiffEst(dataAll(:, 4) == stdNoiseLevel(1) & binaryFeedbackAll == -1)));
%     ciCW_Data(nn, 2, :) = bootci(1000, @nanmedian, angleDiffEst(dataAll(:, 4) == stdNoiseLevel(2) & binaryFeedbackAll == 1));
%     ciCCW_Data(nn, 2, :) = bootci(1000, @nanmedian, abs(angleDiffEst(dataAll(:, 4) == stdNoiseLevel(2) & binaryFeedbackAll == -1)));
    
    
    %% Compute the mean estimate of models
    flagSC = 1; % 1: self-conditioned model
               % 0: standard Bayes
    includeIncongruentTrials = 0;
    dstep = 0.1;
    paramsAll = paramsAllSubject(nn, :);
    lapseRate = paramsAll(3);

    % stimulus orientation
    thetaStim = [-12:2:0 5:5:30]; % 
    thetaStim = round(thetaStim, -log10(dstep));

    % sensory noise
    stdSensory = paramsAll(1:2);

    % memory recall noise
    stdMemory = paramsAll(6);
    stdMemoryIncorrect = 2 * stdMemory;

    % motor noise;
    stdMotor = paramsAll(7);

    % priors
    smoothFactor = paramsAll(8);

    % LOOP - noise levels
    pCw = paramsAll(9);
    pC = [pCw, 1-pCw]'; % [cw ccw]
    pthcw = paramsAll(4);
    pthccw = paramsAll(5); % paramsAll(4)

    rangeth = [-85 85];
    th = rangeth(1):dstep:rangeth(2);
    th = round(th, -log10(dstep));
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
        % orientation noise
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
        % memory noise
        pmmGm = exp(-((MM_m-repmat(m, nmm, 1)).^2)./(2*stdMemory^2)); 
        pmmGm = pmmGm./(repmat(sum(pmmGm,1),nmm,1));   

        % attention marginalization: compute distribution only over those ms that lead to cw decision!
        pmmGthChcw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(1,:)',1,length(thetaStim)));
        b = repmat(a',1,length(thetaStim)) .* pmmGthChcw(indKeepCw, :);        
        pthhGthChcw = interp1(EthChcw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChcw = conv2(pthhGthChcw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChcw(pthhGthChcw < 0) = 0; 

        a = 1./gradient(EthChccw,dstep);
        % attention marginalization: compute distribution only over those ms that lead to cw decision!
        pmmGthChccw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(2,:)',1,length(thetaStim)));        
        b = repmat(a',1,length(thetaStim)) .* pmmGthChccw(indKeepCcw, :);        
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

        %% Prior only
        pthhGthChcw = repmat(normpdf(th', pthccw/2, stdMotor), 1, length(thetaStim));
        pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);   
        pthhGthChcw = pthhGthChcw  .* repmat(PChGtheta_lapse(1,:),nth,1);

        pthhGthChccw = repmat(normpdf(th', pthcw/2, stdMotor), 1, length(thetaStim)) .* repmat(PChGtheta_lapse(2,:),nth,1); 
        pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1); 
        pthhGthChccw =  pthhGthChccw .* repmat(PChGtheta_lapse(2,:),nth,1); 
           

        if includeIncongruentTrials == 0
            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th'<= 0, :) = 0;
            pthhGthChcw(th'> 0, :) = 0;
        end

        % remove 'correct' trials
        pthhGthChcw(:, thetaStim > 0) = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        pthhGthChccw(:, thetaStim < 0) = 0; 

        pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  

        meanCW = th * pthhGthChcw_norm;
        meanCW(isnan(meanCW)) = [];
        meanCCW = th * pthhGthChccw_norm;
        meanCCW(isnan(meanCCW)) = [];
        weightCW = PChGtheta_lapse(1,thetaStim<=0);
        weightCW = weightCW / sum(weightCW);
        weightCCW = PChGtheta_lapse(2,thetaStim>=0);
        weightCCW = weightCCW / sum(weightCCW);
        
        estimateCCW_Prior(nn, kk) = abs(meanCW * weightCW');
        estimateCW_Prior(nn, kk) = meanCCW * weightCCW';

        %% Flip decision
        pmmGth = exp(-((MM_th-THmm).^2)./(2*(stdSensory(kk)^2 + stdMemory^2))); % p(mm|th) = N(th, sm^2 + smm^2)
        pmmGth = pmmGth./(repmat(sum(pmmGth,1),nmm,1)); 

        pthGmmChcw = (pmmGth.*repmat(pthGC(2,:),nmm,1))';
        pthGmmChcw = pthGmmChcw./repmat(sum(pthGmmChcw,1),nth,1);
        pthGmmChcw(isnan(pthGmmChcw)) = 0;

        pthGmmChccw = (pmmGth.*repmat(pthGC(1,:),nmm,1))';
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
        % memory noise
        pmmGm = exp(-((MM_m-repmat(m, nmm, 1)).^2)./(2*stdMemory^2)); 
        pmmGm = pmmGm./(repmat(sum(pmmGm,1),nmm,1));   

        % attention marginalization: compute distribution only over those ms that lead to cw decision!
        pmmGthChcw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(1,:)',1,length(thetaStim)));
        b = repmat(a',1,length(thetaStim)) .* pmmGthChcw(indKeepCw, :);        

        pthhGthChcw = interp1(EthChcw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChcw = conv2(pthhGthChcw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChcw(pthhGthChcw < 0) = 0; 

        a = 1./gradient(EthChccw,dstep);
        % attention marginalization: compute distribution only over those ms that lead to cw decision!
        pmmGthChccw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(2,:)',1,length(thetaStim)));        
        b = repmat(a',1,length(thetaStim)) .* pmmGthChccw(indKeepCcw, :);        
        pthhGthChccw = interp1(EthChccw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChccw = conv2(pthhGthChccw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChccw(pthhGthChccw < 0) = 0; 
        pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
        pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            

        if includeIncongruentTrials == 0
            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th'<= 0, :) = 0;
            pthhGthChcw(th'> 0, :) = 0;
        end

        % remove 'correct' trials
        pthhGthChcw(:, thetaStim > 0) = 0;     
        pthhGthChccw(:, thetaStim < 0) = 0;
        pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
        pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
        pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
        
        meanCW = th * pthhGthChcw_norm;
        meanCW(meanCW == 0) = [];
        meanCCW = th * pthhGthChccw_norm;
        meanCCW(meanCCW == 0) = [];
        weightCW = PChGtheta_lapse(1,thetaStim<=0);
        weightCW = weightCW / sum(weightCW);
        weightCCW = PChGtheta_lapse(2,thetaStim>=0);
        weightCCW = weightCCW / sum(weightCCW);
        
        estimateCCW_FlipDecision(nn, kk) = abs(meanCW * weightCW');
        estimateCW_FlipDecision(nn, kk) = meanCCW * weightCCW';

        %% Resample
        % Likelihood centers on mr, variance: sum of sensory and memory           
        pmrGth = exp(-((MR_th-THmr).^2)./(2*(stdSensory(kk)^2 + stdMemory^2)));
        pmrGth = pmrGth./(repmat(sum(pmrGth,1),nmr,1)); 

        pthGmrChcw = (pmrGth.*repmat(pthGC(2,:),nmr,1))';
        pthGmrChcw = pthGmrChcw./repmat(sum(pthGmrChcw,1),nth,1);
        pthGmrChcw(isnan(pthGmrChcw)) = 0;

        pthGmrChccw = (pmrGth.*repmat(pthGC(1,:),nmr,1))';
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
        pmrGmth = exp(-((MR_m-repmat(m, nmr, 1)).^2)./(2*(stdSensory(kk)^2 + stdMemory^2))); 

        pmrGmthChcw = pmrGmth;
        pmrGmthChcw(mr > 0, :) = 0;
        pmrGmthChcw = pmrGmthChcw./(repmat(sum(pmrGmthChcw,1),nmr,1));

        pmrGmthChccw = pmrGmth;
        pmrGmthChccw(mr < 0, :) = 0;
        pmrGmthChccw = pmrGmthChccw./(repmat(sum(pmrGmthChccw,1),nmr,1));

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
            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th'<= 0, :) = 0;
            pthhGthChcw(th'> 0, :) = 0;
        end

        % remove 'correct' trials
        pthhGthChcw(:, thetaStim > 0) = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
        pthhGthChccw(:, thetaStim < 0) = 0;   

        pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
        pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
        pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
        
        meanCW = th * pthhGthChcw_norm;
        meanCW(meanCW == 0) = [];
        meanCCW = th * pthhGthChccw_norm;
        meanCCW(meanCCW == 0) = [];
        weightCW = PChGtheta_lapse(1,thetaStim<=0);
        weightCW = weightCW / sum(weightCW);
        weightCCW = PChGtheta_lapse(2,thetaStim>=0);
        weightCCW = weightCCW / sum(weightCCW);
        
        estimateCCW_Resample(nn, kk) = abs(meanCW * weightCW');
        estimateCW_Resample(nn, kk) = meanCCW * weightCCW';  
        
        %% Variance only
        % Compute the estimate
        std_combined = sqrt(stdSensory(kk)^2 + stdMemory^2);
        pthhGthChcw = repmat(normpdf(th', -std_combined, stdMotor), 1, length(thetaStim));
        pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);   
        pthhGthChcw = pthhGthChcw  .* repmat(PChGtheta_lapse(1,:),nth,1);

        pthhGthChccw = repmat(normpdf(th', std_combined, stdMotor), 1, length(thetaStim)) .* repmat(PChGtheta_lapse(2,:),nth,1); 
        pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1); 
        pthhGthChccw =  pthhGthChccw .* repmat(PChGtheta_lapse(2,:),nth,1); 
           

        if includeIncongruentTrials == 0
            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th'<= 0, :) = 0;
            pthhGthChcw(th'> 0, :) = 0;
        end

        % remove 'correct' trials
        pthhGthChcw(:, thetaStim > 0) = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
        pthhGthChccw(:, thetaStim < 0) = 0;    
        
        pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
        pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
        pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
        
        meanCW = th * pthhGthChcw_norm;
        meanCW(meanCW == 0) = [];
        meanCCW = th * pthhGthChccw_norm;
        meanCCW(meanCCW == 0) = [];
        weightCW = PChGtheta_lapse(1,thetaStim<=0);
        weightCW = weightCW / sum(weightCW);
        weightCCW = PChGtheta_lapse(2,thetaStim>=0);
        weightCCW = weightCCW / sum(weightCCW);
        
        estimateCCW_Variance(nn, kk) = abs(meanCW * weightCW');
        estimateCW_Variance(nn, kk) = meanCCW * weightCCW'; 
        
        %% Flip decision + more memory noise
        pmmGth = exp(-((MM_th-THmm).^2)./(2*(stdSensory(kk)^2 + stdMemoryIncorrect^2))); % p(mm|th) = N(th, sm^2 + smm^2)
        pmmGth = pmmGth./(repmat(sum(pmmGth,1),nmm,1)); 

        pthGmmChcw = (pmmGth.*repmat(pthGC(2,:),nmm,1))';
        pthGmmChcw = pthGmmChcw./repmat(sum(pthGmmChcw,1),nth,1);
        pthGmmChcw(isnan(pthGmmChcw)) = 0;

        pthGmmChccw = (pmmGth.*repmat(pthGC(1,:),nmm,1))';
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
        % memory noise
        pmmGm = exp(-((MM_m-repmat(m, nmm, 1)).^2)./(2*stdMemoryIncorrect^2)); 
        pmmGm = pmmGm./(repmat(sum(pmmGm,1),nmm,1));   

        % attention marginalization: compute distribution only over those ms that lead to cw decision!
        pmmGthChcw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(1,:)',1,length(thetaStim)));
        b = repmat(a',1,length(thetaStim)) .* pmmGthChcw(indKeepCw, :);        

        pthhGthChcw = interp1(EthChcw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChcw = conv2(pthhGthChcw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChcw(pthhGthChcw < 0) = 0; 

        a = 1./gradient(EthChccw,dstep);
        % attention marginalization: compute distribution only over those ms that lead to cw decision!
        pmmGthChccw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(2,:)',1,length(thetaStim)));        
        b = repmat(a',1,length(thetaStim)) .* pmmGthChccw(indKeepCcw, :);        
        pthhGthChccw = interp1(EthChccw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChccw = conv2(pthhGthChccw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChccw(pthhGthChccw < 0) = 0; 
        pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
        pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            

        if includeIncongruentTrials == 0
            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th'<= 0, :) = 0;
            pthhGthChcw(th'> 0, :) = 0;
        end

        % remove 'correct' trials
        pthhGthChcw(:, thetaStim > 0) = 0;                                                                                                                         
        pthhGthChccw(:, thetaStim < 0) = 0;    
        
        pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
        pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
        pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
        
        meanCW = th * pthhGthChcw_norm;
        meanCW(meanCW == 0) = [];
        meanCCW = th * pthhGthChccw_norm;
        meanCCW(meanCCW == 0) = [];
        weightCW = PChGtheta_lapse(1,thetaStim<=0);
        weightCW = weightCW / sum(weightCW);
        weightCCW = PChGtheta_lapse(2,thetaStim>=0);
        weightCCW = weightCCW / sum(weightCCW);
        
        estimateCCW_FlipDecisionAddMem(nn, kk) = abs(meanCW * weightCW');
        estimateCW_FlipDecisionAddMem(nn, kk) = meanCCW * weightCCW';  
        
        %% Weight LH width by surprise (KL divergence)
        % Scale the LH width by KL divergence
        log_base = 3;        
        scale_factor = PCGm(2,:).*(log2(PCGm(2,:)./PCGm(1,:)) / log2(log_base)) + PCGm(1,:).*(log2(PCGm(1,:)./PCGm(2,:)) / log2(log_base));
        stdSensory_scale = sqrt(1+ scale_factor) * stdSensory(kk);
        pmGth = exp(-((M-THm).^2)./(2*stdSensory_scale.^2)');
    
        pmmGm = exp(-((MM_m-repmat(m, nmm, 1)).^2)./(2*stdMemory^2)); 
        pmmGm = pmmGm./(repmat(sum(pmmGm,1),nmm,1));   
        pmmGth = pmmGm * pmGth;
        
        pthGmmChcw = (pmmGth.*repmat(pthGC(2,:),nmm,1))';
        pthGmmChcw = pthGmmChcw./repmat(sum(pthGmmChcw,1),nth,1);
        pthGmmChcw(isnan(pthGmmChcw)) = 0;

        pthGmmChccw = (pmmGth.*repmat(pthGC(1,:),nmm,1))';
        pthGmmChccw = pthGmmChccw./repmat(sum(pthGmmChccw,1),nth,1);
        pthGmmChccw(isnan(pthGmmChccw)) = 0;

        EthChcw = th * pthGmmChcw;
        EthChccw = th * pthGmmChccw;
        % discard the correct part
        indKeepCw = find(mm>=0);      
        EthChcw = EthChcw(indKeepCw);
        while (sum(diff(EthChcw)>=0) > 0) 
            indDiscardCw = [diff(EthChcw)>=0];
            EthChcw(indDiscardCw) = [];
            indKeepCw(indDiscardCw) = [];
        end
        
        indKeepCcw = find(mm<=0);      
        EthChccw = EthChccw(indKeepCcw);
        while (sum(diff(EthChccw)>=0) >0)
            indDiscardCcw = [false diff(EthChccw)>=0];
            EthChccw(indDiscardCcw) = [];
            indKeepCcw(indDiscardCcw) = [];
        end

        a = abs(1./gradient(EthChcw,dstep));
        % memory noise
        pmmGm = exp(-((MM_m-repmat(m, nmm, 1)).^2)./(2*stdMemory^2)); 
        pmmGm = pmmGm./(repmat(sum(pmmGm,1),nmm,1));   

        % attention marginalization: compute distribution only over those ms that lead to cw decision!
        pmGth = exp(-((M-THm).^2)./(2*stdSensory(kk)^2));
        pmGth = pmGth./(repmat(sum(pmGth,1),nm,1)); 
        pmmGthChcw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(1,:)',1,length(thetaStim)));
        pmmGthChcw = pmmGthChcw./(repmat(sum(pmmGthChcw,1),nmm,1));  
        b = repmat(a',1,length(thetaStim)) .* pmmGthChcw(indKeepCw, :);        

        pthhGthChcw = interp1(EthChcw,b,th,'linear');
        pthhGthChcw(isnan(pthhGthChcw)) = 0;
        % add motor noise
        pthhGthChcw = conv2(pthhGthChcw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChcw(pthhGthChcw < 0) = 0; 

        a = abs(1./gradient(EthChccw,dstep));
        % attention marginalization: compute distribution only over those ms that lead to cw decision!
        pmmGthChccw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(2,:)',1,length(thetaStim)));        
        b = repmat(a',1,length(thetaStim)) .* pmmGthChccw(indKeepCcw, :);        
        pthhGthChccw = interp1(EthChccw,b,th,'linear');
        pthhGthChccw(isnan(pthhGthChccw)) = 0;
        % add motor noise
        pthhGthChccw = conv2(pthhGthChccw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChccw(pthhGthChccw < 0) = 0; 
        pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
        pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            

        if includeIncongruentTrials == 0
            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th'<= 0, :) = 0;
            pthhGthChcw(th'> 0, :) = 0;
        end

        % remove 'correct' trials
        pthhGthChcw(:, thetaStim > 0) = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        pthhGthChccw(:, thetaStim < 0) = 0;    
        
        pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
        pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
        pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
        
        meanCW = th * pthhGthChcw_norm;
        meanCW(meanCW == 0) = [];
        meanCCW = th * pthhGthChccw_norm;
        meanCCW(meanCCW == 0) = [];
        weightCW = PChGtheta_lapse(1,thetaStim<=0);
        weightCW = weightCW / sum(weightCW);
        weightCCW = PChGtheta_lapse(2,thetaStim>=0);
        weightCCW = weightCCW / sum(weightCCW);
        
        estimateCCW_Surprise(nn, kk) = abs(meanCW * weightCW');
        estimateCW_Surprise(nn, kk) = meanCCW * weightCCW';         
    end
end


ciCCW_Data(:, 1, 1) = estimateCCW_Data(:, 1) - ciCCW_Data(:, 1, 1);
ciCW_Data(:, 1, 1) = estimateCW_Data(:, 1) - ciCW_Data(:, 1, 1);
ciCCW_Data(:, 2, 1) = estimateCCW_Data(:, 2) - ciCCW_Data(:, 2, 1);
ciCW_Data(:, 2, 1) = estimateCW_Data(:, 2) - ciCW_Data(:, 2, 1);

ciCCW_Data(:, 1, 2) = ciCCW_Data(:, 1, 2) - estimateCCW_Data(:, 1);
ciCW_Data(:, 1, 2) = ciCW_Data(:, 1, 2) - estimateCW_Data(:, 1);
ciCCW_Data(:, 2, 2) = ciCCW_Data(:, 2, 2) - estimateCCW_Data(:, 2);
ciCW_Data(:, 2, 2) = ciCW_Data(:, 2, 2) - estimateCW_Data(:, 2);

%% Plot the mean estimates
figure
colorName = {'SlateGray', 'DarkMagenta', 'DarkOrange', 'DarkGoldenRod', 'SpringGreen', 'Teal', 'DodgerBlue', 'Navy'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
maxPlot = max([estimateCW_Data(:); estimateCCW_Data(:)]) + 3;
minPlot = min([estimateCW_Data(:); estimateCCW_Data(:)]) - 3;

hold on
set(gca, 'FontSize', 30)
nullX = zeros(1, 7);
nullY = zeros(1, 7);
h = NaN(1, 7);
for ii = 1 : 7
    h(ii) = scatter(nullX(ii), nullY(ii), 1, colorIndex(ii, :), 'filled');
end
[~, iconHandle]  = legend(h, 'S1', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', 'Location', 'NorthWest');
legend('boxoff')
for ii = 1 : 7
    iconHandle(ii+7).Children.MarkerSize = 15;
end
axis off

% Compare model and data
estCCW_Data = estimateCCW_Data(:);
estCW_Data = estimateCW_Data(:);
estData = [estCCW_Data; estCW_Data];

estCCW_Prior = estimateCCW_Prior(:);
estCW_Prior  = estimateCW_Prior(:);
estPrior = [estCCW_Prior ; estCW_Prior ];

estCCW_FlipDecision = estimateCCW_FlipDecision(:);
estCW_FlipDecision = estimateCW_FlipDecision(:);
estFlipDecision = [estCCW_FlipDecision; estCW_FlipDecision];

estCCW_Resample = estimateCCW_Resample(:);
estCW_Resample = estimateCW_Resample(:);
estResample = [estCCW_Resample; estCW_Resample];

estCCW_Variance = estimateCCW_Variance(:);
estCW_Variance = estimateCW_Variance(:);
estVariance = [estCCW_Variance; estCW_Variance];

estCCW_FlipDecisionAddMem = estimateCCW_FlipDecisionAddMem(:);
estCW_FlipDecisionAddMem = estimateCW_FlipDecisionAddMem(:);
estFlipDecisionAddMem = [estCCW_FlipDecisionAddMem; estCW_FlipDecisionAddMem];

estCCW_Surprise= estimateCCW_Surprise(:);
estCW_Surprise = estimateCW_Surprise(:);
estSurprise = [estCCW_Surprise; estCW_Surprise];

figure
subplot(2, 3, 1)
hold on
set(gca, 'FontSize', 12)
for ii = 1 : length(subjectIDAll)
    plot([estimateCCW_Data(ii, 1) estimateCW_Data(ii, 1)], [estimateCCW_Prior(ii, 1) estimateCW_Prior(ii, 1)], 'o', 'Color', colorIndex(ii, :))
    plot([estimateCCW_Data(ii, 2) estimateCW_Data(ii, 2)], [estimateCCW_Prior(ii, 2) estimateCW_Prior(ii, 2)], 'x', 'Color', colorIndex(ii, :), 'MarkerSize', 10)
end
plot([minPlot maxPlot], [minPlot maxPlot], 'k--')
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
xlabel('Mean estimate - data (deg)')
ylabel('Mean estimate - model (deg)')
r = round(corr(estData, estPrior, 'type', 'Pearson'), 2);
MSE = round(sum((estData - estPrior).^2) / length(estData), 1);
title (['Prior only, r: ' num2str(r) ', MSE: ' num2str(MSE)])

subplot(2, 3, 2)
hold on
set(gca, 'FontSize', 12)
for ii = 1 : length(subjectIDAll)
    plot([estimateCCW_Data(ii, 1) estimateCW_Data(ii, 1)], [estimateCCW_FlipDecision(ii, 1) estimateCW_FlipDecision(ii, 1)], 'o', 'Color', colorIndex(ii, :))
    plot([estimateCCW_Data(ii, 2) estimateCW_Data(ii, 2)], [estimateCCW_FlipDecision(ii, 2) estimateCW_FlipDecision(ii, 2)], 'x', 'Color', colorIndex(ii, :), 'MarkerSize', 10)
end
plot([minPlot maxPlot], [minPlot maxPlot], 'k--')
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
xlabel('Mean estimate - data (deg)')
ylabel('Mean estimate - model (deg)')
r = round(corr(estData, estFlipDecision, 'type', 'Pearson'), 2);
MSE = round(sum((estData - estFlipDecision).^2) / length(estData), 1);
title (['Flip Decision, r: ' num2str(r) ', MSE: ' num2str(MSE)])

subplot(2, 3, 3)
hold on
set(gca, 'FontSize', 12)
for ii = 1 : length(subjectIDAll)
    plot([estimateCCW_Data(ii, 1) estimateCW_Data(ii, 1)], [estimateCCW_Resample(ii, 1) estimateCW_Resample(ii, 1)], 'o', 'Color', colorIndex(ii, :))
    plot([estimateCCW_Data(ii, 2) estimateCW_Data(ii, 2)], [estimateCCW_Resample(ii, 2) estimateCW_Resample(ii, 2)], 'x', 'Color', colorIndex(ii, :), 'MarkerSize', 10)
end
plot([minPlot maxPlot], [minPlot maxPlot], 'k--')
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
xlabel('Mean estimate - data (deg)')
ylabel('Mean estimate - model (deg)')
r = round(corr(estData, estResample, 'type', 'Pearson'), 2);
MSE = round(sum((estData - estResample).^2) / length(estData), 1);
title (['Resample, r: ' num2str(r) ', MSE: ' num2str(MSE)])

subplot(2, 3, 4)
hold on
set(gca, 'FontSize', 12)
for ii = 1 : length(subjectIDAll)
    plot([estimateCCW_Data(ii, 1) estimateCW_Data(ii, 1)], [estimateCCW_Variance(ii, 1) estimateCW_Variance(ii, 1)], 'o', 'Color', colorIndex(ii, :))
    plot([estimateCCW_Data(ii, 2) estimateCW_Data(ii, 2)], [estimateCCW_Variance(ii, 2) estimateCW_Variance(ii, 2)], 'x', 'Color', colorIndex(ii, :), 'MarkerSize', 10)
end
plot([minPlot maxPlot], [minPlot maxPlot], 'k--')
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
xlabel('Mean estimate - data (deg)')
ylabel('Mean estimate - model (deg)')
r = round(corr(estData, estVariance, 'type', 'Pearson'), 2);
MSE = round(sum((estData - estVariance).^2) / length(estData), 1);
title (['Variance only, r: ' num2str(r) ', MSE: ' num2str(MSE)])

subplot(2, 3, 5)
hold on
set(gca, 'FontSize', 12)
for ii = 1 : length(subjectIDAll)
    plot([estimateCCW_Data(ii, 1) estimateCW_Data(ii, 1)], [estimateCCW_FlipDecisionAddMem(ii, 1) estimateCW_FlipDecisionAddMem(ii, 1)], 'o', 'Color', colorIndex(ii, :))
    plot([estimateCCW_Data(ii, 2) estimateCW_Data(ii, 2)], [estimateCCW_FlipDecisionAddMem(ii, 2) estimateCW_FlipDecisionAddMem(ii, 2)], 'x', 'Color', colorIndex(ii, :), 'MarkerSize', 10)
end
plot([minPlot maxPlot], [minPlot maxPlot], 'k--')
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
xlabel('Mean estimate - data (deg)')
ylabel('Mean estimate - model (deg)')
r = round(corr(estData, estFlipDecisionAddMem, 'type', 'Pearson'), 2);
MSE = round(sum((estData - estFlipDecisionAddMem).^2) / length(estData), 1);
title (['Flip decision + mem noise, r: ' num2str(r) ', MSE: ' num2str(MSE)])

subplot(2, 3, 6)
hold on
set(gca, 'FontSize', 12)
for ii = 1 : length(subjectIDAll)
    plot([estimateCCW_Data(ii, 1) estimateCW_Data(ii, 1)], [estimateCCW_Surprise(ii, 1) estimateCW_Surprise(ii, 1)], 'o', 'Color', colorIndex(ii, :))
    plot([estimateCCW_Data(ii, 2) estimateCW_Data(ii, 2)], [estimateCCW_Surprise(ii, 2) estimateCW_Surprise(ii, 2)], 'x', 'Color', colorIndex(ii, :), 'MarkerSize', 10)
end
plot([minPlot maxPlot], [minPlot maxPlot], 'k--')
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
xlabel('Mean estimate - data (deg)')
ylabel('Mean estimate - model (deg)')
r = round(corr(estData, estSurprise, 'type', 'Pearson'), 2);
MSE = round(sum((estData - estSurprise).^2) / length(estData), 1);
title (['Surprise-weighted, r: ' num2str(r) ', MSE: ' num2str(MSE)])

% subplot(2, 4, 1)
% hold on
% set(gca, 'FontSize', 12)
% for ii = 1 : length(subjectIDAll)
%     plot(estimateCCW_Data(ii, 1), estimateCW_Data(ii, 1), 'o', 'Color', colorIndex(ii, :))
%     errorbarxy(estimateCCW_Data(ii, 1), estimateCW_Data(ii, 1), ciCCW_Data(ii, 1, 1), ciCCW_Data(ii, 1, 2),...
%                                                                 ciCW_Data(ii, 1, 1), ciCW_Data(ii, 1, 2), {'.', 'k', 'k'})
%     plot(estimateCCW_Data(ii, 2), estimateCW_Data(ii, 2), 'x', 'Color', colorIndex(ii, :), 'MarkerSize', 10)
%     errorbarxy(estimateCCW_Data(ii, 2), estimateCW_Data(ii, 2), ciCCW_Data(ii, 2, 1), ciCCW_Data(ii, 2, 2), ...
%                                                                 ciCW_Data(ii, 2, 1), ciCW_Data(ii, 2, 2), {'.', 'k', 'k'})
% end
% plot([minPlot maxPlot], [minPlot maxPlot], 'k--')
% xlim([minPlot maxPlot])
% ylim([minPlot maxPlot])
% xlabel('Estimate CCW (deg)')
% ylabel('Estimate CW (deg)')
% title('Data')
% 
% subplot(2, 4, 2)
% hold on
% set(gca, 'FontSize', 12)
% for ii = 1 : length(subjectIDAll)
%     plot(estimateCCW_FlipEstimate(ii, 1), estimateCW_FlipEstimate(ii, 1), 'o', 'Color', colorIndex(ii, :))
%     plot(estimateCCW_FlipEstimate(ii, 2), estimateCW_FlipEstimate(ii, 2), 'x', 'Color', colorIndex(ii, :), 'MarkerSize', 10)
% end
% plot([minPlot maxPlot], [minPlot maxPlot], 'k--')
% xlim([minPlot maxPlot])
% ylim([minPlot maxPlot])
% xlabel('Estimate CCW (deg)')
% ylabel('Estimate CW (deg)')
% title('Flip Estimate')
% 
% subplot(2, 4, 3)
% hold on
% set(gca, 'FontSize', 12)
% for ii = 1 : length(subjectIDAll)
%     plot(estimateCCW_FlipDecision(ii, 1), estimateCW_FlipDecision(ii, 1), 'o', 'Color', colorIndex(ii, :))
%     plot(estimateCCW_FlipDecision(ii, 2), estimateCW_FlipDecision(ii, 2), 'x', 'Color', colorIndex(ii, :), 'MarkerSize', 10)
% end
% plot([minPlot maxPlot], [minPlot maxPlot], 'k--')
% xlim([minPlot maxPlot])
% ylim([minPlot maxPlot])
% xlabel('Estimate CCW (deg)')
% ylabel('Estimate CW (deg)')
% title('Flip Decision')
% 
% subplot(2, 4, 4)
% hold on
% set(gca, 'FontSize', 12)
% for ii = 1 : length(subjectIDAll)
%     plot(estimateCCW_Resample(ii, 1), estimateCW_Resample(ii, 1), 'o', 'Color', colorIndex(ii, :))
%     plot(estimateCCW_Resample(ii, 2), estimateCW_Resample(ii, 2), 'x', 'Color', colorIndex(ii, :), 'MarkerSize', 10)
% end
% plot([minPlot maxPlot], [minPlot maxPlot], 'k--')
% xlim([minPlot maxPlot])
% ylim([minPlot maxPlot])
% xlabel('Estimate CCW (deg)')
% ylabel('Estimate CW (deg)')
% title('Resample')
