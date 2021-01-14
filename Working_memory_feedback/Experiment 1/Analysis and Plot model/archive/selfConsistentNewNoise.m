% self-consistent observer model with memory recall and motor noise
% for orientation stimulus (infinite space approximation)
% new memory model: m_m is sample from distribution centered on theta with
% variance being combined sensory and memory variance.
% correct and incorrect trials
%
% astocker
% 01.07.2021 

clf;

%% MODEL
model = {'conditioned','conditionedNewNoise'};
%model = {'conditionedNewNoise'};
%model = 'conditionedNewNoiseIncorr';
% model = 'conditionedANDresampled';


%% PARMS
% sensory noise
sm = [8 5 2];
% memory recall noise
smm = 8;
% motor noise;
smo = 2;
% lapse rate;
lapse = 0.002;

% fixed grids
rangeth = [-42 42];
stimrange = [-21 21];
dstep = 0.1; % fixed - adjustable depending on experiment
th = rangeth(1):dstep:rangeth(2);
nth = length(th);

% grid size for m and mm
nm = 500;
nmm = nm;

%% loop models
for kk=1:length(model)
    
    % generating priors
    PC = [0.5, 0.5]'; % [cw ccw]
    pthcw = 40;
    pthccw = -40;

    pthGC = zeros(2,nth);
    PCGth = zeros(2,nth);
    ind = find(th>=0 & th<=pthcw);
    pthGC(1,ind) = 0.5/length(ind);
    PCGth(1,ind) = 1;
    ind = find(th<=0 & th>=pthccw);
    pthGC(2,ind) = 0.5/length(ind);
    PCGth(2,ind) = 1;
    
%% loop noise levels
for k=1:length(sm)  
    
    % adaptive integration raster
    p = 0.999;
    a = norminv(1-p,rangeth(1),sm(k));
    b = norminv(p, rangeth(2),sm(k));
    m = a:(b-a)/(nm-1):b; 

    a = norminv(1-p,m(1),smm);
    b = norminv(p, m(end),smm);
    mm = a:(b-a)/(nmm-1):b;

    M = repmat(m',1,nth);
    MM = repmat(mm',1,nm);

    % generating measurment distributions
    pmGth = exp(-((M-th).^2)./(2*sm(k)^2));
    pmGth = pmGth./(repmat(sum(pmGth,1),nm,1));
    
  
    %% inference    
    % 1: categorical judgment
    PCGm = (pthGC * pmGth') .* repmat(PC,1,nm);
    PCGm = PCGm./(repmat(sum(PCGm,1),2,1));
    % max posterior decision
    PChGm = round(PCGm);
    % marginalization
    PChGth = PChGm * pmGth;

    % 2: estimation
    
    switch model{kk}
        case{'conditioned'}
            % memory model (old)
            pmmGm = exp(-((MM-m).^2)./(2*smm^2));
            pmmGm = pmmGm./(repmat(sum(pmmGm,1),nmm,1)); 
            pmmGth = pmmGm * pmGth;
            
            pmmGthChcw = pmmGm * (pmGth.*repmat(PChGm(1,:)',1,nth));
            pmmGthChcw = pmmGthChcw./repmat(sum(pmmGthChcw,1),nmm,1);% sum = 1

            pmmGthChccw = pmmGm * (pmGth.*repmat(PChGm(2,:)',1,nth));
            pmmGthChccw = pmmGthChccw./repmat(sum(pmmGthChccw,1),nmm,1);% sum = 1
            
            pthGCh = pthGC; % identical; defintion just for formal reasons

            % posterior cw
            pthGmmChcw = (pmmGth.*repmat(pthGCh(1,:),nmm,1))';
            pthGmmChcw = pthGmmChcw./repmat(sum(pthGmmChcw,1),nth,1);
            pthGmmChcw(isnan(pthGmmChcw)) = 0; % sum = 1

            [EthChcw,ind] = unique(th * pthGmmChcw); % function of mm and m! (via decision)
            % marginalization
            a = 1./gradient(EthChcw,dstep);
            b = repmat(a',1,nth) .* pmmGthChcw(ind,:);
            pthhGthChcw = interp1(EthChcw,b,th,'linear','extrap');
           
            % posterior ccw
            pthGmmChccw = (pmmGth.*repmat(pthGCh(2,:),nmm,1))';
            pthGmmChccw = pthGmmChccw./repmat(sum(pthGmmChccw,1),nth,1);
            pthGmmChccw(isnan(pthGmmChccw)) = 0; % sum = 1

            [EthChccw,ind] = unique(th * pthGmmChccw);
            % marginalization
            a = 1./gradient(EthChccw,dstep);
            b = repmat(a',1,nth) .* pmmGthChccw(ind,:);
            pthhGthChccw = interp1(EthChccw,b,th,'linear','extrap');
           
            
        case{'conditionedNewNoise'}  
             % memory model (new)
            pmmGth = exp(-((M-th).^2)./(2*(sm(k)^2+smm^2)));
            pmmGth = pmmGth./(repmat(sum(pmmGth,1),nmm,1));
                        
            pthGCh = pthGC; % identical; defintion just for formal reasons

            % posterior
            pthGmmChcw = (pmmGth.*repmat(pthGCh(1,:),nmm,1))';
            pthGmmChcw = pthGmmChcw./repmat(sum(pthGmmChcw,1),nth,1);
            pthGmmChcw(isnan(pthGmmChcw)) = 0; % sum = 1

            [EthChcw,ind] = unique(th * pthGmmChcw); %
            % marginalization
            a = 1./gradient(EthChcw,dstep);
            b = repmat(a',1,nth) .* pmmGth(ind,:); % marginalized over mm but not dependent on m!
            pthhGthChcw = interp1(EthChcw,b,th,'linear','extrap');

            % posterior
            pthGmmChccw = (pmmGth.*repmat(pthGCh(2,:),nmm,1))';
            pthGmmChccw = pthGmmChccw./repmat(sum(pthGmmChccw,1),nth,1);
            pthGmmChccw(isnan(pthGmmChccw)) = 0; % sum = 1

            [EthChccw,ind] = unique(th * pthGmmChccw);
            % marginalization
            a = 1./gradient(EthChccw,dstep);
            b = repmat(a',1,nth) .* pmmGth(ind,:);
            pthhGthChccw = interp1(EthChccw,b,th,'linear','extrap');

    end

    
    % add motor noise
    pdf('norm',th,0,smo);
    pthhGthChcw = conv2(pthhGthChcw,pdf('norm',th,0,smo)','same');
    pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not
    
    pthhGthChccw = conv2(pthhGthChccw,pdf('norm',th,0,smo)','same');
    pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);
   

    % total distributions
    % correct trials (i.e. decision correct)
    pthhGthCorr = pthhGthChcw.*repmat(PChGth(1,:),nth,1) + pthhGthChccw.*repmat(PChGth(2,:),nth,1); % sum = 1
    
    % incorrect trials (i.e. decision incorrect, then corrected according
    % to feedback)
   
    pthhGthIncorr = pthhGthChcw.*repmat(PChGth(2,:),nth,1) + pthhGthChccw.*repmat(PChGth(1,:),nth,1); % incorrect trials; sum 1
    pthhGthIncorrw = pthhGthChcw.*repmat(PChGth(2,:),nth,1).*repmat(PCGth(1,:),nth,1) + pthhGthChccw.*repmat(PChGth(1,:),nth,1).*repmat(PCGth(2,:),nth,1); % weighted with probablity of occuring; sum not 1


    
    %% plot
    figure(1);
    px = 3+length(sm);
    py = length(model)*2;
    
    ind = find(th >= stimrange(1) & th <= stimrange(2));
    nthstim = length(ind);

    % correct decision
    subplot(py,px,1+(kk-1)*px*2);
    title('p decision cw');
    plot(th,PChGth(1,:),'r-');
    hold on;
    if k==1
        plot([rangeth(1) rangeth(2)],[0.5 0.5],'k:');
        plot([0 0],[0 1],'k:');
    end
    axis([stimrange(1) stimrange(2) 0 1]);
    
    % incorrect decision
    subplot(py,px,(1+px)+(kk-1)*px*2);
    title('p decision ccw');
    plot(th,PChGth(2,:),'r-');
    hold on;
    if k==1
        plot([rangeth(1) rangeth(2)],[0.5 0.5],'k:');
        plot([0 0],[0 1],'k:');
    end
    axis([stimrange(1) stimrange(2) 0 1]);

    % densities - correct
    ax=subplot(py,px,2+(k-1)+(kk-1)*2*px);
    contour(pthhGthCorr(:,ind),5);
    colormap(ax,'parula');
    hold on;
    axis xy;
    xticklabels({});
    yticklabels({});
    plot([1 nthstim],[ind(1) ind(end)],'r--');
    plot([1 nthstim],[round(nth/2) round(nth/2)],'k:');
    plot([round(nthstim/2) round(nthstim/2)],[1 nth],'k:');
    
    % densities (weighted) - incorrect
    ax=subplot(py,px,2+px+(k-1)+(kk-1)*2*px);
    contour(pthhGthIncorrw(:,ind),5);
    colormap(ax,'parula');
    hold on;
    axis xy;
    xticklabels({});
    yticklabels({});
    plot([1 nthstim],[ind(1) ind(end)],'r--');
    plot([1 nthstim],[round(nth/2) round(nth/2)],'k:');
    plot([round(nthstim/2) round(nthstim/2)],[1 nth],'k:');
    
    % estimates - correct
    ax=subplot(py,px,2+length(sm)+(kk-1)*2*px);
    title('mean estimates')
    pthres = 0.1;
    ind = find(PChGth(1,:)>pthres);
    
    mthhGthChcw = th * pthhGthChcw;
    mthhGthChccw = th * pthhGthChccw;
    plot(th(ind),mthhGthChcw(ind),'c-','linewidth',2);
    hold on;
    ind = find(PChGth(2,:)>pthres);
    plot(th(ind),mthhGthChccw(ind),'g-','linewidth',2);
    axis([stimrange(1) stimrange(2) rangeth(1) rangeth(2)]);
    if k==1
        plot(th,zeros(1,nth),'k:');
        plot([rangeth(1) rangeth(2)],[rangeth(1) rangeth(2)],'k--');
        plot([0 0],[rangeth(1) rangeth(2)],'k:');
    end

    % estimates - incorrect 
    ax=subplot(py,px,2+length(sm)+px+(kk-1)*2*px);
    title('mean estimates')
    pthres = 0.1;
 
    
    a = pthhGthChcw.*repmat(PChGth(2,:),nth,1);
    a = a./repmat(sum(a,1),nth,1);
    
    %mthhGthChcw = th * (pthhGthChcw.*repmat(PChGth(2,:),nth,1));
     mthhGthChcw = th * a;
        
     
    b = pthhGthChccw.*repmat(PChGth(1,:),nth,1);
    b = b./repmat(sum(b,1),nth,1);
    mthhGthChccw = th * (pthhGthChccw.*repmat(PChGth(1,:),nth,1));
    mthhGthChccw = th * b;
    % normalize
%     pthhGthIncorrw = pthhGthIncorrw./repmat(sum(pthhGthIncorrw,1),nth,1); % normalize
%     q= find(isnan(pthhGthIncorrw));
%     pthhGhtIncorrw(q) = 0;
    
    %mthhGth = th * pthhGthIncorr;
       ind = find(th>0 & PChGth(2,:)>pthres);
    plot(th(ind),mthhGthChcw(ind),'c-','linewidth',2);
    hold on;
    ind = find(th<0 & PChGth(1,:)>pthres);
    plot(th(ind),mthhGthChccw(ind),'g-','linewidth',2);
    axis([stimrange(1) stimrange(2) rangeth(1) rangeth(2)]);
    if k==1
        plot(th,zeros(1,nth),'k:');
        plot([rangeth(1) rangeth(2)],[rangeth(1) rangeth(2)],'k--');
        plot([0 0],[rangeth(1) rangeth(2)],'k:');
    end
    
    % variance - correcct

    % variance - incorrect
end
end
