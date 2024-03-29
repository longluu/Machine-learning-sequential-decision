%%%%%%%%%%%%%%%%%%%%%%% Compute the LLH of the models from bootstrap params %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Fix optimal binning model %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Fit params from correct - resample %%%%%%%%%%%%%%%%%%%%%%%

%% Plot the bootstrap params
subjID = {'ll'; 'an'; 'ep'; 'jp'; 'kc'; 'average'};
selectInd = [2 3 6 5];
lowerCI = NaN(length(subjID), length(selectInd));
upperCI = NaN(length(subjID), length(selectInd));
paramsBootstrap = cell(1, length(subjID));
maxIter = 0;
for kk = 1 : length(subjID);
    fileName = ['FitResult-Bootstrap-' subjID{kk} '.txt'];
%     fileName = ['FitResult-' subjID{kk} '.txt'];
    fileID = fopen(fileName);
    paramAll = textscan(fileID,'%f %f %f %f %f %f %f %f %f','CommentStyle','//');
    paramsBootstrap{kk} = cell2mat(paramAll(2:end-1));
    if maxIter < size(paramsBootstrap{kk}, 1)
        maxIter = size(paramsBootstrap{kk}, 1);
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
    for ii = 1:length(selectInd)
        % Read data
        tempParam = paramAll{selectInd(ii)};
        confInterval= prctile(tempParam,[2.5 97.5]);
        lowerCI(kk, ii) = confInterval(1);
        upperCI(kk, ii) = confInterval(2);
        
        % Plot histogram
        if plotFig
            subplot(1, 5, ii)
            hold on
            set(gca, 'FontSize', fontSize)
            hist(tempParam, 20)
            title(['Mean = ' num2str(round(mean(tempParam),2)) ', 95% CI = [' num2str(round(confInterval(1),2)) ' - ' num2str(round(confInterval(2),2)) ']'])
        end
    end
end

paramModel =       [2.6500    6.0895           0.0000     22.2852     1.6506   0.9414    2.0976;
                    3.0023    9.7384           0.0000     34.4053     0.0615   0.9480    3.1069;
                    4.6136   10.4165           0.0000     29.8375     0.1325   0.9940    3.8106;
                    7.7094   11.9114           0.0000     55.7419     0.0083   0.2850    3.8551;
                    5.1033   10.3703           0.0000     46.6421     4.7921   0.8187    3.3313;
                    5.2200   9.5650            0          48.1703     0.0548   0.9768    3.3234];

noiseSensoryExp1 = paramModel(:, 1:2);
noiseMemoryExp1 = paramModel(:, 5);
noiseAll = [noiseSensoryExp1 noiseMemoryExp1];
priorRange = paramModel(:, 4);

% Plot the parameters with bars
fontSize = 23;
colorName = {'Crimson', 'DarkOrange', 'Teal', 'DodgerBlue'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

figure;
hold on
subplot(1, 2, 1);
errorBarGraph(noiseAll, lowerCI(:, 1:3), upperCI(:, 1:3), colorIndex)
box off

subplot(1, 2, 2);
errorBarGraph(priorRange, lowerCI(:, 4), upperCI(:, 4), colorIndex)
box off

%% Compute the subject LLH
tic
logLH_AllModel = NaN(length(subjID), 6, maxIter+1);
logLH_Data = NaN(length(subjID), 1);
addMemNoise = 0;
paramsFit = paramModel;

dstep = 0.1;
rangeth = [-70 70];
th = rangeth(1):dstep:rangeth(2);
th = round(th*(1/dstep))/(1/dstep);
nth = length(th);
numBin = 31;                
for nn = 1 : length(subjID)
    subjectID = subjID{nn};
    
    paramSubj = [paramsBootstrap{nn}; paramsFit(nn, :)];
    for bb = 1 : maxIter+1 % run only the last iteration for no-bootstrap plot
        %% LLH of oracle model
        [~, ~, ~, estimateData, angleDiff, ~] = dataForBootstrap(subjectID, 0, 0);
        logLH = 0;
        binCenter = linspace(-30, 30, numBin);
        deltaBin = diff(binCenter(1:2));
        binEdge = [min(th) binCenter(1:end-1)+deltaBin max(th)];
        binCountAll = cell(size(estimateData));
        for ii = 1 : size(estimateData, 1)
            for jj = 1 : size(estimateData, 2)
                binCount = hist(estimateData{ii, jj}, binCenter);
                pEmpirical = binCount / sum(binCount);
                tempLogLH = binCount .* log(pEmpirical);
                tempLogLH = nansum(tempLogLH);
                if ~isnan(tempLogLH)
                    logLH = logLH +  tempLogLH;
                end
                binCountAll{ii, jj} = binCount;
            end
        end
        logLH_Data(nn, bb) = logLH;       
        
        %% LLH of other models
        flagSC = 1; % 1: self-conditioned model
                   % 0: standard Bayes
        includeIncongruentTrials = 0;

        paramsAll = paramSubj(bb, :);
        lapseRate = paramsAll(3);

        % stimulus orientation
        thetaStim = angleDiff; % 
        thetaStim = round(thetaStim*(1/dstep))/(1/dstep);

        % sensory noise
        stdSensory = paramsAll(1:2);

        % memory recall noise
        stdMemory = paramsAll(5);

        % motor noise;
        stdMotor = paramsAll(7);

        % priors
        smoothFactor = paramsAll(6);

        % LOOP - noise levels
        pC = [0.5, 0.5]'; % [cw ccw]
        pthcw = paramsAll(4);
        pthccw = -paramsAll(4); % paramsAll(4)


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

        logLH_Model = NaN(2, size(logLH_AllModel, 2));

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

            pthhGthChcw_Incorrect = pthhGthChcw;
            pthhGthChccw_Incorrect = pthhGthChccw;

            % remove correct trials
            pthhGthChcw_Incorrect(:, thetaStim > 0) = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          hcw(:, thetaStim < 0) = 0;
            pthhGthChccw_Incorrect(:, thetaStim < 0) = 0;

            % flip the estimate
            pthhGthChcw_Incorrect = flipud(pthhGthChcw_Incorrect);
            pthhGthChccw_Incorrect = flipud(pthhGthChccw_Incorrect);

            %% Model 1b (Flip estimate)
            pthhGthChcw = pthhGthChcw_Incorrect;
            pthhGthChccw = pthhGthChccw_Incorrect;
            pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
            pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
            pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
            pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
            pthhANDth_incorrect = pthhGthChcw_norm + pthhGthChccw_norm;

            % compute LLH
            rangeCollapse = round(length(thetaStim)/2);
            logLikelihoodEstimate = 0;
            tempEstimateModelX = th;

            for jj = rangeCollapse : length(thetaStim)
                tempEstimateModelY = pthhANDth_incorrect(:, jj); 
                tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
                tempestimateDataX = estimateData{kk,jj-rangeCollapse+1};
                pBin = NaN(1, length(binCenter));
                for ii = 1 : length(binCenter)
                    pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdge(ii) & tempEstimateModelX < binEdge(ii+1)));
                end
                binCount = binCountAll{kk,jj-rangeCollapse+1};
                pBin(pBin == 0) = NaN;
                logLikelihoodEstimate = logLikelihoodEstimate + nansum(binCount .* log(pBin));             
            end
            logLH_Model(kk, 1) = logLikelihoodEstimate;


            %% Model 2a1 (Flip decision, Resample: memory, No rejection)
            pmmGth = exp(-((MM_th-THmm).^2)./(2*(stdSensory(kk)^2 + stdMemory^2 + addMemNoise^2))); % p(mm|th) = N(th, sm^2 + smm^2)
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
            pmmGm = exp(-((MM_m-repmat(m, nmm, 1)).^2)./(2*stdMemory^2 + addMemNoise^2)); 
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
            pthhGthChcw(:, thetaStim > 0) = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          hcw(:, thetaStim < 0) = 0;
            pthhGthChccw(:, thetaStim < 0) = 0;
            pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
            pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
            pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
            pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
            pthhANDth_incorrect = pthhGthChcw_norm + pthhGthChccw_norm;

            % compute LLH
            rangeCollapse = round(length(thetaStim)/2);
            logLikelihoodEstimate = 0;
            tempEstimateModelX = th;

            for jj = rangeCollapse : length(thetaStim)
                tempEstimateModelY = pthhANDth_incorrect(:, jj); 
                tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
                tempestimateDataX = estimateData{kk,jj-rangeCollapse+1};
                pBin = NaN(1, length(binCenter));
                for ii = 1 : length(binCenter)
                    pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdge(ii) & tempEstimateModelX < binEdge(ii+1)));
                end
                binCount = binCountAll{kk,jj-rangeCollapse+1};
                pBin(pBin == 0) = NaN;
                logLikelihoodEstimate = logLikelihoodEstimate + nansum(binCount .* log(pBin));             
            end
            logLH_Model(kk, 2) = logLikelihoodEstimate;

            %% Model 2a2 (Resample: memory, With rejection)
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
            pmrGmth = exp(-((MR_m-repmat(m, nmr, 1)).^2)./(2*(stdMemory^2))); 

            pmrGmthChcw = pmrGmth;
            pmrGmthChcw(mr > 0, :) = 0;
            % put the tail with all 0 to 1 (deal with small memory noise)
            indZero = sum(pmrGmthChcw, 1) == 0;
            pmrGmthChcw(mr < 0, indZero) = 1;
            pmrGmthChcw = pmrGmthChcw./(repmat(sum(pmrGmthChcw,1),nmr,1));
            pmrGmthChcw(mr < 0, indZero) = 1e-50;

            pmrGmthChccw = pmrGmth;
            pmrGmthChccw(mr < 0, :) = 0;
            % put the tail with all 0 to 1 (deal with small memory noise)
            indZero = sum(pmrGmthChccw, 1) == 0;
            pmrGmthChccw(mr > 0, indZero) = 1;        
            pmrGmthChccw = pmrGmthChccw./(repmat(sum(pmrGmthChccw,1),nmr,1));
            pmrGmthChccw(mr > 0, indZero) = 1e-50;

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
            pthhANDth_incorrect = pthhGthChcw_norm + pthhGthChccw_norm;

            % compute LLH
            rangeCollapse = round(length(thetaStim)/2);
            logLikelihoodEstimate = 0;
            tempEstimateModelX = th;

            for jj = rangeCollapse : length(thetaStim)
                tempEstimateModelY = pthhANDth_incorrect(:, jj); 
                tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
                tempestimateDataX = estimateData{kk,jj-rangeCollapse+1};
                pBin = NaN(1, length(binCenter));
                for ii = 1 : length(binCenter)
                    pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdge(ii) & tempEstimateModelX < binEdge(ii+1)));
                end
                binCount = binCountAll{kk,jj-rangeCollapse+1};
                pBin(pBin == 0) = NaN;
                logLikelihoodEstimate = logLikelihoodEstimate + nansum(binCount .* log(pBin));             
            end
            logLH_Model(kk, 3) = logLikelihoodEstimate; 

            %% Model 2b1 (Resample: sensory + memory, No rejection)
            pmmGth = exp(-((MM_th-THmm).^2)./(2*(stdSensory(kk)^2 + stdMemory^2 + addMemNoise^2))); % p(mm|th) = N(th, sm^2 + smm^2)
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
            pmmGm = exp(-((MM_m-repmat(m, nmm, 1)).^2)./(2*(stdSensory(kk)^2 + stdMemory^2))); 
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
            pthhGthChcw(:, thetaStim > 0) = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          hcw(:, thetaStim < 0) = 0;
            pthhGthChccw(:, thetaStim < 0) = 0;
            pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
            pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
            pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
            pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
            pthhANDth_incorrect = pthhGthChcw_norm + pthhGthChccw_norm;

            % compute LLH
            rangeCollapse = round(length(thetaStim)/2);
            logLikelihoodEstimate = 0;
            tempEstimateModelX = th;

            for jj = rangeCollapse : length(thetaStim)
                tempEstimateModelY = pthhANDth_incorrect(:, jj); 
                tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
                tempestimateDataX = estimateData{kk,jj-rangeCollapse+1};
                pBin = NaN(1, length(binCenter));
                for ii = 1 : length(binCenter)
                    pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdge(ii) & tempEstimateModelX < binEdge(ii+1)));
                end
                binCount = binCountAll{kk,jj-rangeCollapse+1};
                pBin(pBin == 0) = NaN;
                logLikelihoodEstimate = logLikelihoodEstimate + nansum(binCount .* log(pBin));             
            end
            logLH_Model(kk, 4) = logLikelihoodEstimate;
            
            %% Model 2b2 (Resample: sensory + memory, With rejection)
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
            pthhANDth_incorrect = pthhGthChcw_norm + pthhGthChccw_norm;

            % compute LLH
            rangeCollapse = round(length(thetaStim)/2);
            logLikelihoodEstimate = 0;
            tempEstimateModelX = th;

            for jj = rangeCollapse : length(thetaStim)
                tempEstimateModelY = pthhANDth_incorrect(:, jj); 
                tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
                tempestimateDataX = estimateData{kk,jj-rangeCollapse+1};
                pBin = NaN(1, length(binCenter));
                for ii = 1 : length(binCenter)
                    pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdge(ii) & tempEstimateModelX < binEdge(ii+1)));
                end
                binCount = binCountAll{kk,jj-rangeCollapse+1};
                pBin(pBin == 0) = NaN;
                logLikelihoodEstimate = logLikelihoodEstimate + nansum(binCount .* log(pBin));             
            end
            logLH_Model(kk, 5) = logLikelihoodEstimate;    

            %% Model 1a (Prior only)
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
            pthhGthChcw(:, thetaStim > 0) = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          hcw(:, thetaStim < 0) = 0;
            pthhGthChccw(:, thetaStim < 0) = 0;  
            pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
            pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
            pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
            pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
            pthhANDth_incorrect = pthhGthChcw_norm + pthhGthChccw_norm;

            % compute LLH
            rangeCollapse = round(length(thetaStim)/2);
            logLikelihoodEstimate = 0;
            tempEstimateModelX = th;

            for jj = rangeCollapse : length(thetaStim)
                tempEstimateModelY = pthhANDth_incorrect(:, jj); 
                tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
                tempestimateDataX = estimateData{kk,jj-rangeCollapse+1};
                pBin = NaN(1, length(binCenter));
                for ii = 1 : length(binCenter)
                    pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdge(ii) & tempEstimateModelX < binEdge(ii+1)));
                end
                binCount = binCountAll{kk,jj-rangeCollapse+1};
                pBin(pBin == 0) = NaN;
                logLikelihoodEstimate = logLikelihoodEstimate + nansum(binCount .* log(pBin));             
            end
            logLH_Model(kk, 6) = logLikelihoodEstimate;         

        end
        logLH_AllModel(nn, :, bb) = nansum(logLH_Model, 1);
    end
end
toc

%% Plot the -LLH
% Normalize the logLH by data (upper bound) and Flip estimate (lower bound)
logLH_DataRep = repmat(logLH_Data, 1, 1, 6);
logLH_DataRep = permute(logLH_DataRep, [1 3 2]);
lowerBound = repmat(logLH_AllModel(:, 6, :), 1, 6, 1);
logLH_ModelNorm = (logLH_AllModel - lowerBound) ./ (logLH_DataRep - lowerBound);

meanLogLH = squeeze(logLH_ModelNorm(:, 1:5, end));
lowerConfLogLH = NaN(length(subjID), 5);
upperConfLogLH = NaN(length(subjID), 5);
for ii = 1 : length(subjID)
    for jj = 1:5
        tempLLH = squeeze(logLH_ModelNorm(ii, jj, :));
        tempLLH(isnan(tempLLH)) = [];
        confInterval= prctile(tempLLH(1:end-1),[2.5 97.5]);
        lowerConfLogLH(ii, jj) = confInterval(1);
        upperConfLogLH(ii, jj) = confInterval(2);
    end
end

maxY = 1;
minY = min(lowerConfLogLH(:)) - 0.01;
colorName = {'DarkOrange', 'Teal', 'Purple', 'Red', 'DodgerBlue', };
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

figure
errorBarGraph(meanLogLH, lowerConfLogLH, upperConfLogLH, colorIndex)
