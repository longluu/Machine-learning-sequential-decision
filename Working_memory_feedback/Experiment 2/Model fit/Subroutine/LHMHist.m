% function to compute metrics  ...

function [zz, f0] = LHMHist(data_xx)

target_metric = 'L1';

data_xx = data_xx(:);
min_data_xx = min(data_xx);
max_data_xx = max(data_xx);
n_xx = numel(data_xx);

sq_n_xx = n_xx ^ 0.5;
v1 = (1:sq_n_xx).';
v2 = n_xx./v1;
Nbin_range = unique(round(cat(1,v1,v2)));

done = 0;
while(done == 0)
    
    done = 1;
    n_br = length(Nbin_range);
    metric_vv = zeros(n_br,2);
    
    for i_br = 1:n_br,
        
        nBins = Nbin_range(i_br);
        %ss = sprintf('\t\t# of Bins : %d',nBins); disp(ss);
        binE = linspace(min_data_xx,max_data_xx,(nBins+1));
        binCounts = myhistc(data_xx,binE);
        
        adata_xx_0 = myApproxDataFromHist(binE,binCounts,'nearest');
        adata_xx_1 = myApproxDataFromHist(binE,binCounts,'linear');
        
        if (strcmpi(target_metric,'L1'))
            m1 = sum(abs(data_xx - adata_xx_0));
            m2 = sum(abs(data_xx - adata_xx_1));
        end
        
        if (strcmpi(target_metric,'L2'))
            m1 = sum(abs(data_xx - adata_xx_0).^2);
            m2 = sum(abs(data_xx - adata_xx_1).^2);
        end
        
        if (strcmpi(target_metric,'LInf'))
            m1 = max(abs(data_xx - adata_xx_0));
            m2 = max(abs(data_xx - adata_xx_1));
        end
        
        metric_vv(i_br,1) = m1;
        metric_vv(i_br,2) = m2;
        
    end
    
end

bin_vv = Nbin_range;
nn_metric_vv = (metric_vv - min(metric_vv(:)))./(max(metric_vv(:)) - min(metric_vv(:)));
nn_bin_vv = (bin_vv - min(bin_vv))./(max(bin_vv) - min(bin_vv));

ndist_1 = (nn_bin_vv.^2 + nn_metric_vv(:,1).^2).^0.5;
ndist_2 = (nn_bin_vv.^2 + nn_metric_vv(:,2).^2).^0.5;

[~,ip_ii_NN] = min(ndist_1);
[~,ip_ii_LL] = min(ndist_2);

bn1 = min([bin_vv(ip_ii_LL) bin_vv(ip_ii_NN)]);
bn2 = max([bin_vv(ip_ii_LL) bin_vv(ip_ii_NN)]);
bin_vv1 = bn1:bn2;
metric_rr = zeros((bn2-bn1+1),1);

for bin_i=bn1:bn2,
    
    metric_rr((bin_i-bn1+1),1) = myHistRoughness(data_xx,bin_i);
    
end

[~,min_rr_ii] = min(metric_rr);
opt_bins = bin_vv1(min_rr_ii);
zz = opt_bins;

binE = linspace(min_data_xx,max_data_xx,(opt_bins+1));
binCounts = myhistc(data_xx,binE);

f0 = figure;

hold on;
bar(binE,binCounts,'histc');
hold off;

r_xx = max_data_xx - min_data_xx;
xlim([(min_data_xx - 0.1*r_xx)  (max_data_xx + 0.1*r_xx)]);
ylim([0 max(binCounts)*1.2]);
grid on;


% --------------------------------------------------------------

function [rd_xx] = myApproxDataFromHist(binE,binC,method)

nBins = length(binE)-1;
nData = sum(binC);
rd_xx = zeros(nData,1);

cur_ind = 0;
if (strcmpi(method,'nearest'))
    for i=1:nBins,
        if (binC(i)>0)
            temp_ind_vv = cur_ind + (1:binC(i)).';
            rd_xx(temp_ind_vv) = (binE(i)+binE(i+1))/2;
            cur_ind = cur_ind + binC(i);
        end
    end
end
if (strcmpi(method,'linear'))
    for i=1:nBins,
        if (binC(i)>0)
            temp_ind_vv = cur_ind + (1:binC(i)).';
            temp_data = linspace(binE(i),binE(i+1),(binC(i)+2));
            rd_xx(temp_ind_vv) = temp_data(2:(end-1));
            cur_ind = cur_ind + binC(i);
        end
    end
end

% --------------------------------------------------------------

% function to compute histogram roughness ...

function [hr] = myHistRoughness(data_xx, nBins)

min_data_xx = min(data_xx);
max_data_xx = max(data_xx);

bin_ww = (max_data_xx - min_data_xx)/nBins;
binE = linspace(min_data_xx,max_data_xx,(nBins+1)).';
binM = (binE(1:(end-1)) + binE(2:end))/2;
bin_cc = myhistc(data_xx,binE);
bin_cc = bin_cc(1:(end-1));
yy = bin_cc / (sum(bin_cc)*bin_ww);

yy_d2 = yy((1):(end-2)) + yy((3):(end)) - 2*yy((2):(end-1));
hr = sum(yy_d2.^2)*bin_ww;

% --------------------------------------------------------------

% slightly customized version of HISTC
% values matching last edge are mapped into the previous bin ...

function [BIN_COUNTS, DATA_LABELS] = myhistc(data_xx,binE)

data_xx = data_xx(:);
binE = binE(:);

n_xx = numel(data_xx);
nBins = numel(binE) - 1;

lbl_data = zeros(n_xx,1);
nn_data = zeros(nBins+1,1);

for i=1:nBins,
    if (i<nBins)
        temp_ii = ((data_xx >= binE(i)) & (data_xx < binE(i+1)));
    else
        temp_ii = ((data_xx >= binE(i)) & (data_xx <= binE(i+1)));
    end
    lbl_data(temp_ii) = i;
    nn_data(i) = sum(double(temp_ii));
end

BIN_COUNTS = nn_data;
DATA_LABELS = lbl_data;


