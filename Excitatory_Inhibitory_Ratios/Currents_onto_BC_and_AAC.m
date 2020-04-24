%% Calculate All Excitatory/Inhibitory Ratios onto BCs and AACs
%  Melisa Gumus
%  May 2018

%% Load Data From Netclamp Results
clear all
close all
clc

g = fullfile('/Users','melisagumus','Documents',...
    'labs', 'skinnerlab','netclamp_results','pvbasket',{...
    'pvbasket_332810_1000';...
    'pvbasket_333500_1000';...
    'pvbasket_333776_1000';...
    'pvbasket_334466_1000';...
    'pvbasket_335018_1000';...
    'pvbasket_335432_1000';...
    'pvbasket_335846_1000';...
    'pvbasket_336260_1000';...
    'pvbasket_332948_1000';...
    'pvbasket_333086_1000';...
    'pvbasket_333224_1000';...
    'pvbasket_333638_1000';...
    'pvbasket_333914_1000';...
    'pvbasket_338192_1000';...
    'pvbasket_338054_1000'...
    },{...
    'mytrace_332810_syns.dat';...
    'mytrace_333500_syns.dat';...
    'mytrace_333776_syns.dat';...
    'mytrace_334466_syns.dat';...
    'mytrace_335018_syns.dat';...
    'mytrace_335432_syns.dat';...
    'mytrace_335846_syns.dat';...
    'mytrace_336260_syns.dat';...
    'mytrace_332948_syns.dat';...
    'mytrace_333086_syns.dat';...
    'mytrace_333224_syns.dat';...
    'mytrace_333638_syns.dat';...
    'mytrace_333914_syns.dat';...
    'mytrace_338192_syns.dat';...
    'mytrace_338054_syns.dat'...
    });

f = fullfile('/Users','melisagumus','Documents',...
    'labs', 'skinnerlab','netclamp_results','aac',{...
    'AAC_0_1000';...
    'AAC_36_1000';...
    'AAC_180_1000';...
    'AAC_288_1000';...
    'AAC_360_1000';...
    'AAC_468_1000';...
    'AAC_576_1000';...
    'AAC_720_1000';...
    'AAC_828_1000';...
    'AAC_900_1000';...
    'AAC_1008_1000';...
    'AAC_1152_1000';...
    'AAC_1224_1000';...
    'AAC_1332_1000';...
    'AAC_1404_1000'...
    },{...
    'mytrace_0_syns.dat';...
    'mytrace_36_syns.dat';...
    'mytrace_180_syns.dat';...
    'mytrace_288_syns.dat';...
    'mytrace_360_syns.dat';...
    'mytrace_468_syns.dat';...
    'mytrace_576_syns.dat';...
    'mytrace_720_syns.dat';...
    'mytrace_828_syns.dat';...
    'mytrace_900_syns.dat';...
    'mytrace_1008_syns.dat';...
    'mytrace_1152_syns.dat';...
    'mytrace_1224_syns.dat';...
    'mytrace_1332_syns.dat';...
    'mytrace_1404_syns.dat'...
    });

h = fullfile('/Users','melisagumus','Documents',...
    'labs', 'skinnerlab','netclamp_results','bic',{...
    'BiC_1470_1000';...
    'BiC_1580_1000';...
    'BiC_1635_1000';...
    'BiC_1855_1000';...
    'BiC_2020_1000';...
    'BiC_2130_1000';...
    'BiC_2350_1000';...
    'BiC_2460_1000';...
    'BiC_2570_1000';...
    'BiC_2900_1000';...
    'BiC_3010_1000';...
    'BiC_3175_1000';...
    'BiC_3340_1000';...
    'BiC_3505_1000';...
    'BiC_3615_1000'...
    },{...
    'mytrace_1470_syns.dat';...
    'mytrace_1580_syns.dat';...
    'mytrace_1635_syns.dat';...
    'mytrace_1855_syns.dat';...
    'mytrace_2020_syns.dat';...
    'mytrace_2130_syns.dat';...
    'mytrace_2350_syns.dat';...
    'mytrace_2460_syns.dat';...
    'mytrace_2570_syns.dat';...
    'mytrace_2900_syns.dat';...
    'mytrace_3010_syns.dat';...
    'mytrace_3175_syns.dat';...
    'mytrace_3340_syns.dat';...
    'mytrace_3505_syns.dat';...
    'mytrace_3615_syns.dat'...
    });


%% Write Data on Matrix
alldataBC = [];
for m = 1:1:15
    temp_data = readtable(g{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldataBC = [alldataBC temp_data];
end

dataBC = mat2cell(alldataBC, 40000, ...
    [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]);

alldataAAC = [];
for m = 1:1:15
    temp_data = readtable(f{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldataAAC = [alldataAAC temp_data];
end

dataAAC = mat2cell(alldataAAC, 40000, ...
    [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]);

alldataBiC = [];
for m = 1:1:15
    temp_data = readtable(h{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldataBiC = [alldataBiC temp_data];
end

dataBiC = mat2cell(alldataBiC, 40000, ...
    [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]);

%% Creates a big table consists of inputs from BiC, PYR, and BC... onto BC

M_on_BC = [];
current_BC_on_BC = [];
current_PYR_on_BC = [];
current_BiC_on_BC = [];
current_cck = [];
current_ivy = [];
current_ngf = [];
current_olm = [];
current_sca = [];

for m = 1:15  % number of cells
    for k = 2:12  % number of input
        if k == 3
            temp_current_BiC = dataBC{m}(:,k);
            [pks, locs] = findpeaks(dataBC{m}(:,k),'MinPeakDistance',3000); % peak detection
            temp_BiC = dataBC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BiC(element,:) = 0;
            end
            BiC_on_BC = temp_BiC;
        elseif k == 8
            temp_current_PYR = dataBC{m}(:,k);
            [pks, locs] = findpeaks(-dataBC{m}(:,k),'MinPeakDistance',3000); % peak detection
            temp_PYR = dataBC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_PYR(element,:) = 0;
            end
            PYR_on_BC = temp_PYR;
        elseif k == 9
            temp_current_BC = dataBC{m}(:,k);
            [pks, locs] = findpeaks(dataBC{m}(:,k),'MinPeakDistance',3000); % peak detection
            %findpeaks(data{m}(:,k),M
            temp_BC = dataBC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BC(element,:) = 0;
            end
            BC_on_BC = temp_BC;
        elseif k == 4
            temp_current_cck = dataBC{m}(:,k);
        elseif k == 5
            temp_current_ivy = dataBC{m}(:,k);
        elseif k == 6
            temp_current_ngf = dataBC{m}(:,k);
        elseif k == 7
            temp_current_olm = dataBC{m}(:,k);
        elseif k == 10
            temp_current_sca = dataBC{m}(:,k);

        end

    end
    current_PYR_on_BC = [current_PYR_on_BC temp_current_PYR];
    current_BiC_on_BC = [current_BiC_on_BC temp_current_BiC];
    current_BC_on_BC = [current_BC_on_BC temp_current_BC];
    current_cck_on_BC = [current_cck temp_current_cck];
    current_ivy_on_BC = [current_ivy temp_current_ivy];
    current_ngf_on_BC = [current_ngf temp_current_ngf];
    current_olm_on_BC = [current_olm temp_current_olm];
    current_sca_on_BC = [current_sca temp_current_sca];

    M_on_BC = [M_on_BC BiC_on_BC PYR_on_BC BC_on_BC];

end

allcellsBC = mat2cell(M_on_BC, 40000, ...
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]);

%% Creates a big table consists of inputs from BiC, PYR, and BC... onto AAC

M = [];
current_BC_on_AAC = [];
current_PYR_on_AAC = [];
current_BiC_on_AAC = [];
for m = 1:15  % number of cells
    for k = 2:12  % number of input
        if k == 3
            temp_current_BiC = dataAAC{m}(:,k);
            [pks, locs] = findpeaks(dataAAC{m}(:,k),'MinPeakDistance',3000); % peak detection
            temp_BiC = dataAAC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BiC(element,:) = 0;
            end
            BiC_on_AAC = temp_BiC;
        elseif k == 8
            temp_current_PYR = dataAAC{m}(:,k);
            [pks, locs] = findpeaks(-dataAAC{m}(:,k),'MinPeakDistance',3000); % peak detection
            temp_PYR = dataAAC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_PYR(element,:) = 0;
            end
            PYR_on_AAC = temp_PYR;
        elseif k == 9
            temp_current_BC = dataAAC{m}(:,k);
            [pks, locs] = findpeaks(dataAAC{m}(:,k),'MinPeakDistance',3000); % peak detection
            %findpeaks(data{m}(:,k),M
            temp_BC = dataAAC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BC(element,:) = 0;
            end
            BC_on_AAC = temp_BC;
            elseif k == 4
                temp_current_cck = dataAAC{m}(:,k);
            elseif k == 5
                temp_current_ivy = dataAAC{m}(:,k);
            elseif k == 6
                temp_current_ngf = dataAAC{m}(:,k);
            elseif k == 7
                temp_current_olm = dataAAC{m}(:,k);
            elseif k == 10
                temp_current_sca = dataAAC{m}(:,k);
        end
    end
    current_PYR_on_AAC = [current_PYR_on_AAC temp_current_PYR];
    current_BiC_on_AAC = [current_BiC_on_AAC temp_current_BiC];
    current_BC_on_AAC = [current_BC_on_AAC temp_current_BC];
    M = [M BiC_on_AAC PYR_on_AAC BC_on_AAC];
    current_cck_on_AAC = [current_cck temp_current_cck];
    current_ivy_on_AAC = [current_ivy temp_current_ivy];
    current_ngf_on_AAC = [current_ngf temp_current_ngf];
    current_olm_on_AAC = [current_olm temp_current_olm];
    current_sca_on_AAC = [current_sca temp_current_sca];
end

allcellsAAC = mat2cell(M, 40000, ...
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]);

%% Creates a big table consists of inputs from BiC, PYR, and BC... onto BiC

M_BiC = [];
current_BC_on_BiC = [];
current_PYR_on_BiC = [];
current_BiC_on_BiC = [];
for m = 1:15  % number of cells
    for k = 2:12  % number of input
        if k == 3
            temp_current_BiC = dataBiC{m}(:,k);
            [pks, locs] = findpeaks(dataBiC{m}(:,k),'MinPeakDistance',3000); % peak detection
            temp_BiC = dataBiC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BiC(element,:) = 0;
            end
            BiC_on_BiC = temp_BiC;
        elseif k == 8
            temp_current_PYR = dataBiC{m}(:,k);
            [pks, locs] = findpeaks(-dataBiC{m}(:,k),'MinPeakDistance',3000); % peak detection
            temp_PYR = dataBiC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_PYR(element,:) = 0;
            end
            PYR_on_BiC = temp_PYR;
        elseif k == 9
            temp_current_BC = dataBiC{m}(:,k);
            [pks, locs] = findpeaks(dataBiC{m}(:,k),'MinPeakDistance',3000); % peak detection
            %findpeaks(data{m}(:,k),M
            temp_BC = dataBiC{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BC(element,:) = 0;
            end
            BC_on_BiC = temp_BC;
            elseif k == 4
                temp_current_cck = dataAAC{m}(:,k);
            elseif k == 5
                temp_current_ivy = dataAAC{m}(:,k);
            elseif k == 6
                temp_current_ngf = dataAAC{m}(:,k);
            elseif k == 7
                temp_current_olm = dataAAC{m}(:,k);
            elseif k == 10
                temp_current_sca = dataAAC{m}(:,k);
        end
    end
    current_PYR_on_BiC = [current_PYR_on_BiC temp_current_PYR];
    current_BiC_on_BiC = [current_BiC_on_BiC temp_current_BiC];
    current_BC_on_BiC = [current_BC_on_BiC temp_current_BC];
    M_BiC = [M_BiC BiC_on_BiC PYR_on_BiC BC_on_BiC];
    current_cck_on_BiC = [current_cck temp_current_cck];
    current_ivy_on_BiC = [current_ivy temp_current_ivy];
    current_ngf_on_BiC = [current_ngf temp_current_ngf];
    current_olm_on_BiC = [current_olm temp_current_olm];
    current_sca_on_BiC = [current_sca temp_current_sca];
end

allcellsBiC = mat2cell(M_BiC, 40000, ...
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]);

%% TRY DIFFERENT COMBINATIONS - IPSCs Only From BC and AAC onto BC and AAC
% Sum all ipsc currents
all_ipsc_on_BC_AAC_combo1 = [];
all_epsc_on_BC_AAC_combo1 = [];

for i = 1:1:15
    for t = 1:1:15

        tot_cur_ipsc_on_BC_combo1 = current_BC_on_BC(:,i);
        tot_cur_ipsc_on_AAC_combo1 = current_BC_on_AAC(:,t);

        tot_cur_ipsc_on_BC_AAC_combo1 = tot_cur_ipsc_on_BC_combo1 + tot_cur_ipsc_on_AAC_combo1;
        all_ipsc_on_BC_AAC_combo1 = [all_ipsc_on_BC_AAC_combo1 tot_cur_ipsc_on_BC_AAC_combo1];

        tot_cur_epsc_on_BC_AAC_combo1 =  current_PYR_on_BC(:,i) + current_PYR_on_AAC(:,t);
        all_epsc_on_BC_AAC_combo1 = [all_epsc_on_BC_AAC_combo1 tot_cur_epsc_on_BC_AAC_combo1];
    end
end
%% Find the peaks of the summed IPSCs from BC and AAC onto BC, AAC - and EPSCs too

peaks_all_PV_on_BC_AAC_combo1 = [];
f1 = figure;
for k = 1:1:225
    [pks, locs] = findpeaks(all_ipsc_on_BC_AAC_combo1(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(all_ipsc_on_BC_AAC_combo1(:,k),'MinPeakDistance',3000);

    temp_cur = all_ipsc_on_BC_AAC_combo1(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_all_PV_on_BC_AAC_combo1 = [peaks_all_PV_on_BC_AAC_combo1 peaks_all];
end

%%
peaks_PYR_on_BC_AAC_combo1 = [];
for k = 1:1:225
    [pks, locs] = findpeaks(-all_epsc_on_BC_AAC_combo1(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(-all_epsc_on_BC_AAC_combo1(:,k),'MinPeakDistance',3000);

    temp_cur = all_epsc_on_BC_AAC_combo1(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_PYR_on_BC_AAC_combo1 = [peaks_PYR_on_BC_AAC_combo1 peaks_all];
end

%% IPSCs from BC and AAC onto BC and AAC

IPSC_all_on_BC_AAC_combo1 = [];
ipsc_all_on_BC_AAC_combo1 = [];
for i = 1:1:225
    pks_ipsc_all = peaks_all_PV_on_BC_AAC_combo1(:,i);
    pks_ipsc_all(pks_ipsc_all == 0) = [];
    ipsc_all_mean = mean(pks_ipsc_all);
    ipsc_all_std = std(pks_ipsc_all);
    ipsc_all_on_BC_AAC_combo1 = [ipsc_all_mean;ipsc_all_std];
    IPSC_all_on_BC_AAC_combo1 = [IPSC_all_on_BC_AAC_combo1 ipsc_all_on_BC_AAC_combo1];
end

%% Find Mean EPSC and SD onto BC, AAC

EPSC_on_BC_AAC_combo1 = [];
epsc_on_BC_AAC_combo1 = [];
for i = 1:1:225
    pks_epsc_all = peaks_PYR_on_BC_AAC_combo1(:,i);
    pks_epsc_all(pks_epsc_all == 0) = [];
    epsc_all_mean = mean(pks_epsc_all);
    epsc_all_std = std(pks_epsc_all);
    epsc_on_BC_AAC_combo1 = [epsc_all_mean;epsc_all_std];
    EPSC_on_BC_AAC_combo1 = [EPSC_on_BC_AAC_combo1 epsc_on_BC_AAC_combo1];
end

%% Excitatory/Inhibitory Ratios on BC, AAC population and BC, AAC population

IPSC_BC_AAC_on_BC_AAC_combo1 = IPSC_all_on_BC_AAC_combo1;

Ratios_BC_AAC =[];
temp_ratio_BC_AAC = [];

E_I_BC_AAC = abs(EPSC_on_BC_AAC_combo1(1,:)./IPSC_BC_AAC_on_BC_AAC_combo1(1,:))';

temp_ratio_BC_AAC = [temp_ratio_BC_AAC E_I_BC_AAC];
Ratios_BC_AAC = array2table(reshape(temp_ratio_BC_AAC,[15,15]));

%% E/I Ratio - BC, AAC to BC, AAC
fig = uitable('Data',Ratios_BC_AAC{:,:},...
    'RowName',[],...
    'ColumnName',[],...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);
%%

less_than_1 = sum(temp_ratio_BC_AAC >1)/225;
