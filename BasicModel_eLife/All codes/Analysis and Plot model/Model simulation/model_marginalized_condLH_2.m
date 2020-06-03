% Density plot for combined incongruent error
% ********** New version - condition the prior AND the likelihood ********** 
% astocker - lluu
% 06.2017
flagSC = 1; % 1: self-conditioned model
           % 0: standard Bayes
flagDecisionGiven = 0;
incongruentType = 3; % 1: Wrong decision, right estimate
                     % 2: Wrong decision, forgot stimulus
                     % 3: Motor noise
                     % 4: Misrember given decision, right estimate
dstep = 0.1;
paramsAllSubject = [2.5611    4.8570    7.1153      0.0000    28.6016      0.8805    0.8461    2.6327;
                    4.0300    5.3505    8.1590      0.0000    41.7649   2.0957    0.0257    1.8431;
                    3.8692    6.4042   12.3882      0.0000    39.5010  14.0727    0.2322    2.9146;
                    3.7041    4.6566   10.9916      0.0000    36.4739   0.2433    0.6481    5.4045;
                    3.9226    5.7843   10.9330      0.0000    39.0534  14.2463    0.1492    5.8917];

% paramsAllSubject = [2.5611    4.8570    7.1153      0.0000    23.1044    0.8805    0.8461    2.6327;
%                     6.3630    8.4075   14.5464      0.0000    22.5672   5.9826    0.7652    2.9350;
%                     4.6736    6.1762    7.9466      0.0000    15.8459   1.2144    0.3681    2.9916;
%                     4.6585    5.5272    6.8023      0.0000    17.4750   4.1846    0.0934    2.5348;
%                     4.3172    4.5344    8.4083      0.0000    33.9063  14.5556    0.5305    2.0456];

                
incongruentErrorWeight = [0.1111 0.7778 5.7778 4.6111 3.8333];
incongruentErrorWeight = incongruentErrorWeight / sum(incongruentErrorWeight);
lapseExp2 = [0.13 5.06 13.6 0.58 2.85];
lapseExp3 = NaN(size(lapseExp2));

% stimulus orientation
thetaStim = -21:0.1:21; % 
thetaStim = round(thetaStim, -log10(dstep));

rangeth = [-42 42];
th = rangeth(1):dstep:rangeth(2);
th = round(th, -log10(dstep));
nth = length(th);

pthhGthAverage = zeros(nth, length(thetaStim));

for ii = 1 : size(paramsAllSubject, 1)
    paramsAll = paramsAllSubject(ii, :);
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
    pthGC_Discrimination(1,:) = TukeyWindow([0 pthcw], 0, smoothFactor, th);
    pthGC_Discrimination(2,:) = TukeyWindow([pthccw 0], 1, smoothFactor, th);

    %% LOOP - noise levels
    pthhGth_Collapse = zeros(nth, length(thetaStim));
    pCongruentCW = zeros(length(stdSensory) , round(length(thetaStim)/2));
%     figure
%     hold on
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

        if incongruentType == 1
            % Wrong decision, right estimate
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

            pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
            pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            

            % modify psychometric curve p(Chat|theta, Congruent) ~ p(Congruent| Chat, theta) * p(Chat|Theta)
            pCongruentGcwTh = sum(pthhGthChcw(th' >= 0, :));
            pCongruentGccwTh = sum(pthhGthChccw(th' <= 0, :));
            PChGtheta_lapse_new = PChGtheta_lapse .* [pCongruentGcwTh; pCongruentGccwTh];
            PChGtheta_lapse_new = PChGtheta_lapse_new ./ repmat(sum(PChGtheta_lapse_new, 1), 2, 1);

            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th'>= 0, :) = 0;
            pthhGthChcw(th'< 0, :) = 0;
        elseif incongruentType == 2
            % Wrong decision, forgot stimulus
            pthhGthChcw = repmat(normpdf(th', pthcw/2, stdMotor), 1, length(thetaStim));
            pthhGthChccw = repmat(normpdf(th', pthccw/2, stdMotor), 1, length(thetaStim));
            pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
            pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            
            
            PChGtheta_lapse_new = PChGtheta_lapse;
        elseif incongruentType == 3
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

            pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
            pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            

            PChGtheta_lapse_new = PChGtheta_lapse;

            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th' <= 0, :) = 0;
            pthhGthChcw(th' > 0, :) = 0;            
        elseif incongruentType == 4
            % Misremeber given decision, right estimate
            pthGC_new = flipud(pthGC);
            pmmGth = exp(-((MM_th-THmm).^2)./(2*(stdSensory(kk)^2 + stdMemory^2))); % p(mm|th) = N(th, sm^2 + smm^2)
            pmmGth = pmmGth./(repmat(sum(pmmGth,1),nmm,1)); 
            pthGmmChcw = (pmmGth.*repmat(pthGC_new(1,:),nmm,1))';
            pthGmmChcw = pthGmmChcw./repmat(sum(pthGmmChcw,1),nth,1);
            pthGmmChcw(isnan(pthGmmChcw)) = 0;

            pthGmmChccw = (pmmGth.*repmat(pthGC_new(2,:),nmm,1))';
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

            pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
            pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);     
            pthhGthChcw(:, thetaStim < 0) = 0;
            pthhGthChccw(:, thetaStim > 0) = 0;
            
            % modify psychometric curve p(Chat|theta, Congruent) ~ p(Congruent| Chat, theta) * p(Chat|Theta)
            pCongruentGcwTh = sum(pthhGthChcw(th' <= 0, :));
            pCongruentGccwTh = sum(pthhGthChccw(th' >= 0, :));
            PChGtheta_lapse_new = PChGtheta_lapse .* [pCongruentGcwTh; pCongruentGccwTh];
            PChGtheta_lapse_new = PChGtheta_lapse_new ./ repmat(sum(PChGtheta_lapse_new, 1), 2, 1);

            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th'<= 0, :) = 0;
            pthhGthChcw(th'> 0, :) = 0;
            pthhGthChcw(isnan(pthhGthChcw)) = 0;
            pthhGthChccw(isnan(pthhGthChccw)) = 0;
            
        end

        pthhGth = pthhGthChcw.*repmat(PChGtheta_lapse_new(1,:),nth,1) + pthhGthChccw.*repmat(PChGtheta_lapse_new(2,:),nth,1);
        pthhAndth = pthhGth;
        pthhGth_Collapse = pthhGth_Collapse + pthhAndth;
%         pCongruentCW(kk, :) = pCongruentGcwTh(round(length(thetaStim)/2):end);
        
        % Plot
%         xRange = [-22 22];
%         indX = find(thetaStim >= xRange(1) & thetaStim <= xRange(2));
%         xMax = length(indX);
%         yRange = [-35 35];
%         indY = find(th >= yRange(1) & th <= yRange(2));
%         yMax = length(indY);
%         thNew = th(indY);
%         indYStart = find(thNew == xRange(1));
%         indYEnd = find(thNew == xRange(2));
%         subplot(1, 3, kk)
%         imagesc(pthhAndth(indY, indX));
%         hold on;
%         axis xy;
%         colormap('gray');
%         plot([1 xMax],[round(yMax/2) round(yMax/2)],'w:', 'LineWidth', 1.5);
%         plot([round(xMax/2) round(xMax/2)],[1 yMax],'w:', 'LineWidth', 1.5);
%         plot([1 xMax],[indYStart indYEnd],'b:', 'LineWidth', 2);
%         set(gca, 'ylim', [1 yMax], 'xlim', [1 xMax], ...
%             'XTick', round(linspace(1,xMax,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
%             'YTick', round(linspace(1,yMax,5)), 'YTickLabel', num2cell(round(linspace(yRange(1),yRange(2),5))), ...
%             'FontSize', 20)
    end

%     lapseExp3(ii) =  lapseExp2(ii) * mean(pCongruentCW(:)); 
    
    % Plot
    figure
    xRange = [-22 22];
    indX = find(thetaStim >= xRange(1) & thetaStim <= xRange(2));
    xMax = length(indX);
    yRange = [-42 42];
    indY = find(th >= yRange(1) & th <= yRange(2));
    yMax = length(indY);
    thNew = th(indY);
    indYStart = find(thNew == xRange(1));
    indYEnd = find(thNew == xRange(2));

    matrixPlot = pthhGth_Collapse(indY, indX);
    imagesc(max(matrixPlot(:)) - matrixPlot);
    hold on;
    axis xy;
    colormap('gray');
    plot([1 xMax],[round(yMax/2) round(yMax/2)],'w:', 'LineWidth', 1.5);
    plot([round(xMax/2) round(xMax/2)],[1 yMax],'w:', 'LineWidth', 1.5);
    plot([1 xMax],[indYStart indYEnd],'b:', 'LineWidth', 2);
    set(gca, 'ylim', [1 yMax], 'xlim', [1 xMax], ...
        'XTick', round(linspace(1,xMax,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
        'YTick', round(linspace(1,yMax,5)), 'YTickLabel', num2cell(round(linspace(yRange(1),yRange(2),5))), ...
        'FontSize', 20)
    xlabel('True orientation (degree)')
    ylabel('Estimated orientation (degree)')        
    
    pthhGthAverage = pthhGthAverage + incongruentErrorWeight(ii) * pthhGth_Collapse;    
end

% Plot
figure
xRange = [-22 22];
indX = find(thetaStim >= xRange(1) & thetaStim <= xRange(2));
xMax = length(indX);
yRange = [-40 40];
indY = find(th >= yRange(1) & th <= yRange(2));
yMax = length(indY);
thNew = th(indY);
indYStart = find(thNew == xRange(1));
indYEnd = find(thNew == xRange(2));

matrixPlot = pthhGthAverage(indY, indX);
imagesc(max(matrixPlot(:)) - matrixPlot);

hold on;
axis xy;
colormap('gray');
plot([1 xMax],[round(yMax/2) round(yMax/2)],'k:', 'LineWidth', 1.5);
plot([round(xMax/2) round(xMax/2)],[1 yMax],'k:', 'LineWidth', 1.5);
plot([1 xMax],[indYStart indYEnd],'b:', 'LineWidth', 2);
set(gca, 'ylim', [1 yMax], 'xlim', [1 xMax], ...
    'XTick', round(linspace(1,xMax,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
    'YTick', round(linspace(1,yMax,5)), 'YTickLabel', num2cell(round(linspace(yRange(1),yRange(2),5))), ...
    'FontSize', 20)
xlabel('True orientation (degree)')
ylabel('Estimated orientation (degree)')        

