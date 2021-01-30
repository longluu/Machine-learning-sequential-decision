%%%%%%%%%%%%%%%%%%%%%%% Compute the LLH of the models from bootstrap params %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Fix optimal binning model %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% ( the optimal bin is computed in computeOptBin.m) %%%%%%%%%%%%%%%%%%%%%%%

%% Plot the bootstrap params
subjectIDAll = {'ll', 'pw', 'eh', 'bh', 'ln', 'at', 'dh'}; % 
selectInd = [2 3 7 5 6 10];
lowerCI = NaN(length(subjectIDAll), length(selectInd));
upperCI = NaN(length(subjectIDAll), length(selectInd));
paramsBootstrap = cell(1, length(subjectIDAll));
path_bootstrap = 'C:\Users\longluu\Documents\GitHub\Machine-learning-sequential-decision\Working_memory_feedback\Experiment 2\Model fit\Fit result\Version 2 (free pC)\Bootstrap\Decision and estimate\';
maxIter = 0;
subjIter = NaN(1, length(subjectIDAll));
for kk = 1 : length(subjectIDAll)
    fileName = [path_bootstrap 'FitResult-Bootstrap-' subjectIDAll{kk} '.txt'];
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

% Extract fit parameter
paramModel = NaN(length(subjectIDAll), 9);
negLLH_1 = NaN(length(subjectIDAll), 1);
path_fitResult = 'C:\Users\longluu\Documents\GitHub\Machine-learning-sequential-decision\Working_memory_feedback\Experiment 2\Model fit\Fit result\Version 2 (free pC)\NoResample\';
for kk = 1 : length(subjectIDAll)
    fileName = [path_fitResult 'FitResult-' subjectIDAll{kk} '-extracted.txt'];
    fileID = fopen(fileName);
    paramAll = textscan(fileID,'%f %f %f %f %f %f %f %f %f %f','CommentStyle','//');
    result_mat = cell2mat(paramAll);
    paramModel(kk, :) = result_mat(end, 2:end);
    negLLH_1(kk, :) = result_mat(end, 1);
    fclose(fileID);
end

noiseSensoryExp1 = paramModel(:, 1:2);
noiseMemoryExp1 = paramModel(:, 6);
noiseAll = [noiseSensoryExp1 noiseMemoryExp1];
priorRange = paramModel(:, 4:5);
pCwAll = paramModel(:, end);

% Plot the parameters with bars
fontSize = 23;
colorName = {'Crimson', 'DarkOrange', 'Teal', 'DodgerBlue'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

figure;
hold on
subplot(1, 3, 1);
errorBarGraph(noiseAll, lowerCI(:, 1:3), upperCI(:, 1:3), colorIndex)
box off

subplot(1, 3, 2);
errorBarGraph(abs(priorRange), abs(lowerCI(:, 4:5)), abs(upperCI(:, 4:5)), colorIndex)
box off

subplot(1, 3, 3);
errorBarGraph(abs(pCwAll), abs(lowerCI(:, 6)), abs(upperCI(:, 6)), colorIndex)
box off

%% Compute the subject LLH
tic
logLH_AllModel = NaN(length(subjectIDAll), 6, maxIter+1);
logLH_Data = NaN(length(subjectIDAll), maxIter+1);
paramsFit =  paramModel;

dstep = 0.1;
rangeth = [-70 70];
th = rangeth(1):dstep:rangeth(2);
th = round(th, -log10(dstep));
nth = length(th);
numBin = 31;
for nn = 1 : length(subjectIDAll)
    subjectID = subjectIDAll{nn};
    
    paramSubj = [paramsBootstrap{nn}; paramsFit(nn, :)];
    for bb = 1: subjIter(nn)+1
        %% LLH of oracle model
        [~, ~, ~, estimateData, ~, ~] = dataForFitting(subjectID, 0, 0);
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

            pthhGthChcw_correct = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
            pthhGthChcw_correct = pthhGthChcw_correct  .* repmat(PChGtheta_lapse(2,:),nth,1);
            pthhGthChccw_correct = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
            pthhGthChccw_correct = pthhGthChccw_correct  .* repmat(PChGtheta_lapse(1,:),nth,1);

            if includeIncongruentTrials == 0
                % modify psychometric curve p(Chat|theta, Congruent) ~ p(Congruent| Chat, theta) * p(Chat|Theta)
                pCongruentGcwTh = sum(pthhGthChcw(th' >= 0, :));
                pCongruentGccwTh = sum(pthhGthChccw(th' <= 0, :));
                PChGtheta_lapse_new = PChGtheta_lapse .* [pCongruentGcwTh; pCongruentGccwTh];
                PChGtheta_lapse_new = PChGtheta_lapse_new ./ repmat(sum(PChGtheta_lapse_new, 1), 2, 1);

                % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
                pthhGthChccw_correct(th'>= 0, :) = 0;
                pthhGthChcw_correct(th'< 0, :) = 0;
            else
                PChGtheta_lapse_new = PChGtheta_lapse;
            end
            pthhGthChccw_correct(:, thetaStim > 0) = 0;
            pthhGthChcw_correct(:, thetaStim < 0) = 0; 

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
            pthhGthChcw = (1-lapseRate) * pthhGthChcw + lapseRate * pthhGthChccw_correct;
            pthhGthChccw(:, thetaStim < 0) = 0;  
            pthhGthChccw = (1-lapseRate) * pthhGthChccw + lapseRate * pthhGthChcw_correct;            
            pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
            pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
            pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
            pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
            pthhANDth_incorrect = pthhGthChcw_norm + pthhGthChccw_norm;
            pthhANDth_incorrect(:, thetaStim == 0) = pthhANDth_incorrect(:, thetaStim == 0)/2;

            % compute LLH
            logLikelihoodEstimate = 0;
            tempEstimateModelX = th;

            for jj = 1 : length(thetaStim)
                tempEstimateModelY = pthhANDth_incorrect(:, jj); 
                tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
                pBin = NaN(1, length(binCenter));
                for ii = 1 : length(binCenter)
                    pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdge(ii) & tempEstimateModelX < binEdge(ii+1)));
                end
                binCount = binCountAll{kk,jj};
                pBin(pBin == 0) = NaN;
                logLikelihoodEstimate = logLikelihoodEstimate + nansum(binCount .* log(pBin));             
            end
            logLH_Model(kk, 1) = logLikelihoodEstimate;
            
            %% Variance only
            % Compute the estimate
            std_combined = sqrt(stdSensory(kk)^2 + 0^2);
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
            pthhGthChcw = (1-lapseRate) * pthhGthChcw + lapseRate * pthhGthChccw_correct;
            pthhGthChccw(:, thetaStim < 0) = 0;  
            pthhGthChccw = (1-lapseRate) * pthhGthChccw + lapseRate * pthhGthChcw_correct;            
            pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
            pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
            pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
            pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
            pthhANDth_incorrect = pthhGthChcw_norm + pthhGthChccw_norm;
            pthhANDth_incorrect(:, thetaStim == 0) = pthhANDth_incorrect(:, thetaStim == 0)/2;

            % compute LLH
            logLikelihoodEstimate = 0;
            tempEstimateModelX = th;

            for jj = 1 : length(thetaStim)
                tempEstimateModelY = pthhANDth_incorrect(:, jj); 
                tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
                pBin = NaN(1, length(binCenter));
                for ii = 1 : length(binCenter)
                    pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdge(ii) & tempEstimateModelX < binEdge(ii+1)));
                end
                binCount = binCountAll{kk,jj};
                pBin(pBin == 0) = NaN;
                logLikelihoodEstimate = logLikelihoodEstimate + nansum(binCount .* log(pBin));             
            end
            logLH_Model(kk, 2) = logLikelihoodEstimate;
            
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
            pthhGthChcw = (1-lapseRate) * pthhGthChcw + lapseRate * pthhGthChccw_correct;
            pthhGthChccw(:, thetaStim < 0) = 0;  
            pthhGthChccw = (1-lapseRate) * pthhGthChccw + lapseRate * pthhGthChcw_correct;            
            pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
            pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
            pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
            pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
            pthhANDth_incorrect = pthhGthChcw_norm + pthhGthChccw_norm;
            pthhANDth_incorrect(:, thetaStim == 0) = pthhANDth_incorrect(:, thetaStim == 0)/2;
            pthhANDth_incorrect = pthhANDth_incorrect / dstep;

            % compute LLH
            logLikelihoodEstimate = 0;
            tempEstimateModelX = th;

            for jj = 1 : length(thetaStim)
                tempEstimateModelY = pthhANDth_incorrect(:, jj); 
                tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
                pBin = NaN(1, length(binCenter));
                for ii = 1 : length(binCenter)
                    pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdge(ii) & tempEstimateModelX < binEdge(ii+1)));
                end
                binCount = binCountAll{kk,jj};
                pBin(pBin == 0) = NaN;
                logLikelihoodEstimate = logLikelihoodEstimate + nansum(binCount .* log(pBin));             
            end
            logLH_Model(kk, 3) = logLikelihoodEstimate;

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
            pthhGthChcw = (1-lapseRate) * pthhGthChcw + lapseRate * pthhGthChccw_correct;
            pthhGthChccw(:, thetaStim < 0) = 0;  
            pthhGthChccw = (1-lapseRate) * pthhGthChccw + lapseRate * pthhGthChcw_correct;            
            pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
            pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
            pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
            pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
            pthhANDth_incorrect = pthhGthChcw_norm + pthhGthChccw_norm;
            pthhANDth_incorrect(:, thetaStim == 0) = pthhANDth_incorrect(:, thetaStim == 0)/2;
            pthhANDth_incorrect = pthhANDth_incorrect / dstep;

            % compute LLH
            logLikelihoodEstimate = 0;
            tempEstimateModelX = th;

            for jj = 1 : length(thetaStim)
                tempEstimateModelY = pthhANDth_incorrect(:, jj); 
                tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
                pBin = NaN(1, length(binCenter));
                for ii = 1 : length(binCenter)
                    pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdge(ii) & tempEstimateModelX < binEdge(ii+1)));
                end
                binCount = binCountAll{kk,jj};
                pBin(pBin == 0) = NaN;
                logLikelihoodEstimate = logLikelihoodEstimate + nansum(binCount .* log(pBin));             
            end
            logLH_Model(kk, 4) = logLikelihoodEstimate;
            
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
            pthhGthChcw = (1-lapseRate) * pthhGthChcw + lapseRate * pthhGthChccw_correct;
            pthhGthChccw(:, thetaStim < 0) = 0;  
            pthhGthChccw = (1-lapseRate) * pthhGthChccw + lapseRate * pthhGthChcw_correct;            
            pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
            pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
            pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
            pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
            pthhANDth_incorrect = pthhGthChcw_norm + pthhGthChccw_norm;
            pthhANDth_incorrect(:, thetaStim == 0) = pthhANDth_incorrect(:, thetaStim == 0)/2;
            pthhANDth_incorrect = pthhANDth_incorrect / dstep;

            % compute LLH
            logLikelihoodEstimate = 0;
            tempEstimateModelX = th;

            for jj = 1 : length(thetaStim)
                tempEstimateModelY = pthhANDth_incorrect(:, jj); 
                tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
                pBin = NaN(1, length(binCenter));
                for ii = 1 : length(binCenter)
                    pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdge(ii) & tempEstimateModelX < binEdge(ii+1)));
                end
                binCount = binCountAll{kk,jj};
                pBin(pBin == 0) = NaN;
                logLikelihoodEstimate = logLikelihoodEstimate + nansum(binCount .* log(pBin));             
            end
            logLH_Model(kk, 5) = logLikelihoodEstimate;  
            
            %% Weight LH width by surprise (KL divergence)
            % Scale the LH width by KL divergence
            log_base = exp(1);        
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
            pthhGthChcw = (1-lapseRate) * pthhGthChcw + lapseRate * pthhGthChccw_correct;
            pthhGthChccw(:, thetaStim < 0) = 0;  
            pthhGthChccw = (1-lapseRate) * pthhGthChccw + lapseRate * pthhGthChcw_correct;            
            pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
            pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
            pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
            pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
            pthhANDth_incorrect = pthhGthChcw_norm + pthhGthChccw_norm;

            % compute LLH
            logLikelihoodEstimate = 0;
            tempEstimateModelX = th;

            for jj = 1 : length(thetaStim)
                tempEstimateModelY = pthhANDth_incorrect(:, jj); 
                tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
                pBin = NaN(1, length(binCenter));
                for ii = 1 : length(binCenter)
                    pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdge(ii) & tempEstimateModelX < binEdge(ii+1)));
                end
                binCount = binCountAll{kk,jj};
                pBin(pBin == 0) = NaN;
                logLikelihoodEstimate = logLikelihoodEstimate + nansum(binCount .* log(pBin));             
            end
            logLH_Model(kk, 6) = logLikelihoodEstimate;              
        end
        logLH_AllModel(nn, :, bb) = nansum(logLH_Model, 1);
    end
end
toc

% Save LLH (for comparison with Flip Decision of the new noise model)
llh_resample = logLH_AllModel(:, 5, end);
save('llh_resample.mat', 'llh_resample')

%% Plot the -LLH
% Normalize the logLH by data (upper bound) and Flip estimate (lower bound)
logLH_DataRep = repmat(logLH_Data, 1, 1, size(logLH_AllModel, 2));
logLH_DataRep = permute(logLH_DataRep, [1 3 2]);
lowerBound = repmat(logLH_AllModel(:, 1, :), 1, size(logLH_AllModel, 2), 1);
logLH_ModelNorm = (logLH_AllModel - lowerBound) ./ (logLH_DataRep - lowerBound);

meanLogLH = squeeze(logLH_ModelNorm(:, 2:end, end));
lowerConfLogLH = NaN(length(subjectIDAll), size(logLH_AllModel, 2)-1);
upperConfLogLH = NaN(length(subjectIDAll), size(logLH_AllModel, 2)-1);
for ii = 1 : length(subjectIDAll)
    for jj = 1:size(logLH_AllModel, 2)-1
        tempLLH = squeeze(logLH_ModelNorm(ii, jj+1, :));
        tempLLH(isnan(tempLLH)) = [];
        confInterval= prctile(tempLLH(1:end-1),[2.5 97.5]);
        lowerConfLogLH(ii, jj) = confInterval(1);
        upperConfLogLH(ii, jj) = confInterval(2);
    end
end

maxY = 1;
minY = min(lowerConfLogLH(:)) - 0.01;
colorName = {'DarkOrange', 'Teal', 'Green', 'DodgerBlue', 'Red'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

figure
errorBarGraph(meanLogLH, lowerConfLogLH, upperConfLogLH, colorIndex)
