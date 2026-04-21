%% This script 

p.data = 'Z:\OMEGA\OMEGA_data';
p.results = 'Z:\Enrique\Waveform\Results';
p.figures = 'Z:\Enrique\Waveform\Figures';

% Toolbox paths
addpath('Z:\Toolbox\fieldtrip-20230118') % Path to Fieldtrip
addpath('Z:\Toolbox\NECOG')
ft_defaults

addpath(genpath('Z:\Toolbox\sourceBOSC'))

%% OMEGA 
% Participant and session list
sub = [1 2 3 4 5 6 7 8 9 11 12 14 15 16 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 37 39 40 41 42 44 45 46 47 48 49 50 51 52 55 56 57 58 59 60 61 62 63 64 65 67 68 69 70 71 72 73 74 75 76 77 78 79 80 84 85 87 88 89 90 91 92 94 95 96 97 98 99 101 102 103 104 105 106 134 145 146 148 149 150 151 152 154 155 156 157 158 159 160 161 165 166 167 168 169 170 171 175 176 177 179 181 184 185 195 197 200 207 208 210 212]';
ses = [1 1 1 1 1 1 1 2 1  2  1  2  1  1  1  2  3  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  3  1  2  1  2  1  2  1  1  2  1  1  1  1  1  1  1  1  2  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 ]';

subs = {} ; sess = {};
for s = 1:length(sub)
    if sub(s) < 10
        subs{s} = ['000' num2str(sub(s))];
    elseif sub(s) >= 10 & sub(s) < 100
        subs{s} = ['00' num2str(sub(s))];
    else
        subs{s} = ['0' num2str(sub(s))];
    end
    sess{s} = ['000' num2str(ses(s))];
end

Nsub  = length(sub);

load source_template_10mm_1925.mat

% Frex parameters 
frex = exp(0.6:0.1:3.7); % 1.8 Hz to 40 Hz
Nfrex = length(frex);

freqbands = [frex(1) frex(9) frex(16) frex(21) frex(25) frex(30)];
freqnames = {'delta', 'theta', 'alpha', 'lowbeta', 'highbeta'};

Nvoxin = length(find(source.inverse.inside));

%% Avgs
load([p.results '\delta\voxshape_s1.mat']) % Load a sample to obtain filenames
features = fieldnames(voxshape{1});

trls = zeros(length(freqnames), Nsub, Nvoxin); 
avgshape = zeros(Nsub, Nvoxin,length(features));
for fb = 1:length(freqnames)
    frqband = freqnames{fb};
    display(frqband)
    for s = 1:Nsub
        load([p.results '\' frqband '\voxshape_s' num2str(s) '.mat'])
        for vx = 1:length(voxshape)
            trls(fb,s,vx) = length(voxshape{vx}.band_amp);
                voxshape{vx}.time_rdsym = abs(voxshape{vx}.time_rdsym-0.5)*2; % 0.5 [0 1]. Rescale to average [0 1]
                voxshape{vx}.time_ptsym = abs(voxshape{vx}.time_ptsym-0.5)*2; % 0.5 [0 1]. Rescale to average [0 1]
                avgshape(s,vx,:) = structfun(@(x) median(double(x)), voxshape{vx});             
        end
    end
    filename = p.results + "\avgshape_" + freqnames{fb} + ".mat";
    save(filename, 'avgshape');
end

%% Statistics: ttest
addpath('Z:\Toolbox\slanCM\slanCM')
addpath('Z:\Toolbox\fieldtrip-20230118') 
freqnames = {'delta', 'theta', 'alpha', 'lowbeta', 'highbeta'};
paramlist = {'period', 'time_rdsym', 'time_ptsym', 'band_amp'};
[~, paramidx] = ismember(paramlist, features);

for fb = 1:length(freqnames)
    mkdir([p.results '/' freqnames{fb}])
    frqband = freqnames{fb};
    load([p.results '\avgshape_' freqnames{fb}])
    for param = 1:length(paramlist)
        % Load avg data
        tmpdata = squeeze(avgshape(:,:,paramidx(param)));
        zData = normalize(tmpdata, 2, 'zscore');

        % Ttest
        [~, pvals,~, stats] = arrayfun(@(v) ttest(zData(:,v), 0), 1:size(tmpdata, 2),'UniformOutput',false);
        pvals = cell2mat(pvals);
        tVals = cellfun(@(s) s.tstat, stats);
        pfdr = mafdr(pvals(:), 'BHFDR', true);
        tValsfdr = tVals;
        tValsfdr(pfdr>.05) = 0;
        
        % Permutation Based Analysis
        nperm = 1000;
        maxTval = zeros(nperm, 1);
        for pp = 1:nperm
            permsign = (rand(Nsub, 1) > 0.5)*2 - 1; 
            permdata = zData .* permsign;
            [~,~,~, permstats] = arrayfun(@(v) ttest(permdata(:,v), 0), 1:size(permdata, 2),'UniformOutput',false);
            permtVals = cellfun(@(s) s.tstat, permstats);
            maxTval(pp) = max(abs(permtVals));
        end

        permTcorr = arrayfun(@(t) (sum(maxTval >= t) + 1) / (nperm + 1), abs(tVals));
        tValspermcor = tVals;
        tValspermcor(find(permTcorr > 0.05)) = NaN;

        sBOSC_nii(tValspermcor,[p.figures '/' freqnames{fb} '/Perm12_' paramlist{param} '_' frqband])

        cfg = [];
        cfg.colmap = slanCM('vik');
        cfg.colim = [-max(abs(tVals)) max(abs(tVals))];
        cfg.interp = 'linear';
        cfg.savefig = 'yes';
        sBOSC_sourcefig(tValspermcor,[p.figures '/' freqnames{fb} '/PermTvals122_' paramlist{param} '_' frqband], cfg)
        close
    end
end
