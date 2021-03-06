%% Calculate Excitatory/Inhibitory Ratios onto PYR Cells
%  Melisa Gumus
%  May 2018


%% Load Data From Netclamp Results
clear all
close all
clc

f = fullfile('~',... # Change this path accordingly 
    'CA1_SimpleDetailed','Excitatory_Inhibitory_Ratios',...
    'Network_Clamp_Results','pyr',{...
    'pyr_29097_1000';...
    'pyr_36884_1000';...
    'pyr_52458_1000';...
    'pyr_68032_1000';...
    'pyr_83606_1000';...
    'pyr_99180_1000';...
    'pyr_106967_1000';...
    'pyr_114754_1000';...
    'pyr_75819_1000';...
    'pyr_44671_1000';...
    'pyr_60245_1000';...
    'pyr_21310_1000';...
    'pyr_91393_1000';...
    'pyr_301642_1000';...
    'pyr_325003_1000'...
    },{...
    'mytrace_29097_syns.dat';...
    'mytrace_36884_syns.dat';...
    'mytrace_52458_syns.dat';...
    'mytrace_68032_syns.dat';...
    'mytrace_83606_syns.dat';...
    'mytrace_99180_syns.dat';...
    'mytrace_106967_syns.dat';...
    'mytrace_114754_syns.dat';...
    'mytrace_75819_syns.dat';...
    'mytrace_44671_syns.dat';...
    'mytrace_60245_syns.dat';...
    'mytrace_21310_syns.dat';...
    'mytrace_91393_syns.dat';...
    'mytrace_301642_syns.dat';...
    'mytrace_325003_syns.dat'...
    });

%% Write Data on Matrix
alldata = [];
for m = 1:1:15
    temp_data = readtable(f{m},'Delimiter','\t');
    temp_data = table2array(temp_data);
    alldata = [alldata temp_data];
end

data = mat2cell(alldata, 40000, ...
    [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]);

%% Creates a big table consists of inputs from AAC, BiC, PYR, and BC... onto PYR

M = [];
current_AAC = [];
current_BC = [];
current_PYR = [];
current_BiC = [];
current_cck = [];
current_ivy = [];
current_ngf = [];
current_olm = [];
current_sca = [];
current_ca3 = [];
current_ec = [];
figure1 = figure;
figure2 = figure;
figure3 = figure;
figure4 = figure;
figure5 = figure;
figure6 = figure;
for m = 1:15  % number of cells
    for k = 2:12  % number of input - excluding first column (time)
        if k ==2
            temp_current_AAC = data{m}(:,k);
            figure(figure1);
            t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
            t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
            t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
            t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
            t.String = ['Peak Detection on IPSCs from AAC onto PYR'];
            subplot(5,3,m);
            [pks, locs] = findpeaks(data{m}(:,k),'MinPeakDistance',3000); % peak detection with interval based threshold
            findpeaks(data{m}(:,k),'MinPeakDistance',3000);
            hold on; 
            title (['PYR Number #' num2str(m)])
            xlabel('Time (ms)')
            ylabel('IPSCs from AACs')
            temp_AAC = data{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_AAC(element,:) = 0;
            end
            AAC = temp_AAC;
        elseif k == 3
            temp_current_BiC = data{m}(:,k);
            figure(figure2);
            t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
            t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
            t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
            t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
            t.String = ['Peak Detection on IPSCs from BiC onto PYR'];
            subplot(5,3,m);
            [pks, locs] = findpeaks(data{m}(:,k),'MinPeakDistance',3000); % peak detection
            findpeaks(data{m}(:,k),'MinPeakDistance',3000);
            hold on; 
            title (['PYR Number #' num2str(m)])
            xlabel('Time (msec)')
            ylabel('IPSC')
            temp_BiC = data{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BiC(element,:) = 0;
            end
            BiC = temp_BiC;
        elseif k == 8
            temp_current_PYR = data{m}(:,k);
            figure(figure3);
            t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
            t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
            t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
            t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
            t.String = ['Peak Detection on EPSCs from PYR onto PYR'];
            subplot(5,3,m);
            [pks, locs] = findpeaks(data{m}(:,k),'MinPeakDistance',3000); % peak detection
            findpeaks(data{m}(:,k),'MinPeakDistance',3000);
            hold on; 
            title (['PYR Number #' num2str(m)])
            xlabel('Time (msec)')
            ylabel('EPSC')
            temp_PYR = data{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_PYR(element,:) = 0;
            end
            PYR = temp_PYR;
        elseif k == 9
            temp_current_BC = data{m}(:,k);
            figure(figure4);
            t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
            t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
            t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
            t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
            t.String = ['Peak Detection on IPSCs from BC onto PYR'];
            subplot(5,3,m);
            [pks, locs] = findpeaks(data{m}(:,k),'MinPeakDistance',3000); % peak detection
            findpeaks(data{m}(:,k),'MinPeakDistance',3000);
            hold on; 
            title (['PYR Number #' num2str(m)])
            xlabel('Time (msec)')
            ylabel('IPSCs')
            temp_BC = data{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_BC(element,:) = 0;
            end
            BC = temp_BC;
        elseif k == 11
            temp_current_ca3 = data{m}(:,k);
            figure(figure5);
            t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
            t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
            t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
            t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
            t.String = ['Peak Detection on EPSCs from CA3 onto PYRs'];
            subplot(5,3,m);
            [pks, locs] = findpeaks(data{m}(:,k),'MinPeakDistance',3000); % peak detection
            findpeaks(data{m}(:,k),'MinPeakDistance',3000);
            hold on; 
            title (['BC Number #' num2str(m)])
            xlabel('Time')
            ylabel('EPSC from CA3')
            temp_ca3 = data{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_ca3(element,:) = 0;
            end 
            CA3 = temp_ca3;
        elseif k == 12
            temp_current_ec = data{m}(:,k);
            figure(figure6);
            t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
            t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
            t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
            t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
            t.String = ['Peak Detection on EPSCs from EC onto PYRs'];
            subplot(5,3,m);
            [pks, locs] = findpeaks(data{m}(:,k),'MinPeakDistance',3000); % peak detection
            findpeaks(data{m}(:,k),'MinPeakDistance',3000);
            hold on; 
            title (['BC Number #' num2str(m)])
            xlabel('Time')
            ylabel('EPSC from EC')
            temp_ec = data{m}(:,k);
            allrows = (1:40000)';
            notpeak = setdiff(allrows,locs);
            for t = 1:1:numel(notpeak)
                element = notpeak(t,:);
                temp_ec(element,:) = 0;
            end 
            EC = temp_ec;
        elseif k == 4
            temp_current_cck = data{m}(:,k);
        elseif k == 5
            temp_current_ivy = data{m}(:,k);
        elseif k == 6
            temp_current_ngf = data{m}(:,k);
        elseif k == 7
            temp_current_olm = data{m}(:,k);
        elseif k == 10
            temp_current_sca = data{m}(:,k);
        end
    end
    current_AAC = [current_AAC temp_current_AAC];
    current_PYR = [current_PYR temp_current_PYR];
    current_BiC = [current_BiC temp_current_BiC];
    current_BC = [current_BC temp_current_BC];
    current_cck = [current_cck temp_current_cck];
    current_ivy = [current_ivy temp_current_ivy];
    current_ngf = [current_ngf temp_current_ngf];
    current_olm = [current_olm temp_current_olm];
    current_sca = [current_sca temp_current_sca];
    current_ca3 = [current_ca3 temp_current_ca3];
    current_ec = [current_ec temp_current_ec];
    
    M = [M AAC BiC PYR BC CA3 EC];
end

allcells = mat2cell(M, 40000, ...
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]);

%% Graph of EPSCs onto PYR
f1 = figure;
for i = 1:1:15
    figure(f1);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Distribution of EPSCs on PYRs'];
    temp = allcells{i}(:,3); 
    temp(temp == 0) = [];   %get rid of zeros
    s1 = subplot(5,3,i);
    p = histfit(temp);
    set(p(1),'facecolor',[0.1 0.6 0.6]);
    set(p(2),'color','k');
    hold on;
    title (['Pyramidal Cell #' num2str(i)]); % name each graph with the corresponding pyr #
    xlabel('EPSC');
    ylabel('Number of EPSCs');
end


%% Graph of IPSCs from AAC onto PYR
f2 = figure;
for i = 1:1:15
    figure(f2);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Distribution of IPSCs from AACs onto PYRs'];
    temp = allcells{i}(:,1);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,3,i)
    p = histfit(temp);
    set(p(1),'facecolor',[0.1 0.6 0.6])
    set(p(2),'color','k')
    hold on
    title (['Pyramidal Cell #' num2str(i)]) % name each graph with the corresponding pyr #
    xlabel('IPSCs')
    ylabel('Number of IPSCs')
end

%suptitle('Distribution of IPSCs from AAC onto PYR Cells')

%% Graph of IPSCs from BiC onto PYR
f3 = figure;
for i = 1:1:15
    figure(f3);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Distribution of IPSCs from BiCs onto PYRs'];
    temp = allcells{i}(:,2);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,3,i)
    p = histfit(temp);
    set(p(1),'facecolor',[0.1 0.6 0.6])
    set(p(2),'color','k')
    hold on
    title (['Pyramidal Cell #' num2str(i)]) % name each graph with the corresponding pyr #
    xlabel('IPSCs')
    ylabel('Number of IPSCs')
end
%suptitle('Distribution of IPSCs from BiC onto PYR Cells')

%% Graph of IPSCs from BC onto PYR
f4 = figure;
for i = 1:1:15
    figure(f4);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Distribution of IPSCs from BCs onto PYRs'];
    temp = allcells{i}(:,2);
    temp(temp == 0) = [];   %get rid of zeros
    subplot(5,3,i)
    p = histfit(temp);
    set(p(1),'facecolor',[0.1 0.6 0.6])
    set(p(2),'color','k')
    hold on
    title (['Pyramidal Cell #' num2str(i)]) % name each graph with the corresponding pyr #
    xlabel('IPSCs')
    ylabel('Number of IPSCs')
end
%suptitle('Distribution of IPSCs from BC onto PYR Cells')

%% Find Mean EPSCs and SD onto PYR
EPSC = [];
epsc = [];
for i = 1:1:15 % number of PYR cells
    pks_epsc = allcells{i}(:,3);
    pks_epsc(pks_epsc==0)=[];
    epsc_mean = mean(pks_epsc);
    epsc_std = std(pks_epsc);
    epsc = [epsc_mean;epsc_std];
    EPSC = [EPSC epsc];
end

EPSC_table = EPSC;
num = (1:15)';
EPSC_table = array2table(EPSC_table');
EPSC_table.num = num;
EPSC_table = [EPSC_table(:,end) EPSC_table(:,1) EPSC_table(:,2)];

EPSC_table.Properties.VariableNames = {'Pyramidal_Cell_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(2,2,1)
EPSC_mean = EPSC(1,:);
EPSC_std = EPSC(2,:);
x = linspace(1,15,length(EPSC_mean));
scatter(x,EPSC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dEPSC_mean=0.1;
text(x+dx, EPSC_mean+dEPSC_mean, c);
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSCs','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_mean,EPSC_std,'b','LineStyle','none')
title('Mean Peak EPSCs onto PYR Cells','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of EPSC on PYR Cells
fig = uitable('Data',EPSC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'Pyramidal Cell Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);


%% IPSCs only from AAC onto PYR
IPSC_AAC = [];
ipsc_AAC = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_AAC = allcells{i}(:,1);
    pks_ipsc_AAC(pks_ipsc_AAC==0)=[];
    ipsc_AAC_mean = mean(pks_ipsc_AAC);
    ipsc_AAC_std = std(pks_ipsc_AAC);
    ipsc_AAC = [ipsc_AAC_mean;ipsc_AAC_std];
    IPSC_AAC = [IPSC_AAC ipsc_AAC];
end

IPSC_AAC_table = IPSC_AAC;
num = (1:15)';
IPSC_AAC_table = array2table(IPSC_AAC_table');
IPSC_AAC_table.num = num;
IPSC_AAC_table = [IPSC_AAC_table(:,end) IPSC_AAC_table(:,1) IPSC_AAC_table(:,2)];

IPSC_AAC_table.Properties.VariableNames = {'Pyramidal_Cell_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(2,2,2)
IPSC_AAC_mean = IPSC_AAC(1,:);
IPSC_AAC_std = IPSC_AAC(2,:);
x = linspace(0,15,length(IPSC_AAC_mean));
scatter(x,IPSC_AAC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_AAC_mean=0.1;
text(x+dx, IPSC_AAC_mean+dIPSC_AAC_mean, c);
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from AAC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_AAC_mean,IPSC_AAC_std,'b','LineStyle','none')
title('Mean Peak IPSC from AAC onto PYR Cells','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSC from AAC onto PYR Cells
fig = uitable('Data',IPSC_AAC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'Pyramidal Cell Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% IPSCs only from BiC onto PYR
IPSC_BiC = [];
ipsc_BiC = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_BiC = allcells{i}(:,2);
    pks_ipsc_BiC(pks_ipsc_BiC==0)=[];
    ipsc_BiC_mean = mean(pks_ipsc_BiC);
    ipsc_BiC_std = std(pks_ipsc_BiC);
    ipsc_BiC = [ipsc_BiC_mean;ipsc_BiC_std];
    IPSC_BiC = [IPSC_BiC ipsc_BiC];
end

IPSC_BiC_table = IPSC_BiC;
num = (1:15)';
IPSC_BiC_table = array2table(IPSC_BiC_table');
IPSC_BiC_table.num = num;
IPSC_BiC_table = [IPSC_BiC_table(:,end) IPSC_BiC_table(:,1) IPSC_BiC_table(:,2)];

IPSC_BiC_table.Properties.VariableNames = {'Pyramidal_Cell_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(2,2,3)
IPSC_BiC_mean = IPSC_BiC(1,:);
IPSC_BiC_std = IPSC_BiC(2,:);
x = linspace(0,15,length(IPSC_BiC_mean));
scatter(x,IPSC_BiC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_BiC_mean=0.1;
text(x+dx, IPSC_BiC_mean+dIPSC_BiC_mean, c);
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BiC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BiC_mean,IPSC_BiC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BiC onto PYR Cells','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSC from BiC onto PYR Cells
fig = uitable('Data',IPSC_BiC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'Pyramidal Cell Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% IPSCs only from BC onto PYR
IPSC_BC = [];
ipsc_BC = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_BC = allcells{i}(:,4);
    pks_ipsc_BC(pks_ipsc_BC==0)=[];
    ipsc_BC_mean = mean(pks_ipsc_BC);
    ipsc_BC_std = std(pks_ipsc_BC);
    ipsc_BC = [ipsc_BC_mean;ipsc_BC_std];
    IPSC_BC = [IPSC_BC ipsc_BC];
end

IPSC_BC_table = IPSC_BC;
num = (1:15)';
IPSC_BC_table = array2table(IPSC_BC_table');
IPSC_BC_table.num = num;
IPSC_BC_table = [IPSC_BC_table(:,end) IPSC_BC_table(:,1) IPSC_BC_table(:,2)];

IPSC_BC_table.Properties.VariableNames = {'Pyramidal_Cell_Number', 'Mean_Peak', 'Standard_Deviation'};

subplot(2,2,4)
IPSC_BC_mean = IPSC_BC(1,:);
IPSC_BC_std = IPSC_BC(2,:);
x = linspace(0,15,length(IPSC_BC_mean));
scatter(x,IPSC_BC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_BC_mean=0.1;
text(x+dx, IPSC_BC_mean+dIPSC_BC_mean, c);
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC from BC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_BC_mean,IPSC_BC_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC onto PYR Cells','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSC from BC onto PYR Cells
fig = uitable('Data',IPSC_BC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'Pyramidal Cell Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% Calculate Mean IPSCs Gathered from AAC, BiC and BC onto PYR

% Sum all IPSC currents
all_ipsc = [];
all_ipsc_together = [];

for i = 1:1:15
    tot_cur_ipsc = current_AAC(:,i) + current_BiC(:,i) + current_BC(:,i);
    tot_ipsc_together = current_BiC(:,i)...
        +current_AAC(:,i)...
        +current_BC(:,i)...
        +current_cck(:,i)...
        +current_ivy(:,i)...
        +current_ngf(:,i)...
        +current_olm(:,i)...
        +current_sca(:,i);
    
    all_ipsc = [all_ipsc tot_cur_ipsc];
    all_ipsc_together = [all_ipsc_together tot_ipsc_together];
end

% Find the peaks of the summed ipsc currents - from BC, BiC, AAC
peaks_all_PV = [];
f1 = figure;
for k = 1:1:15
    figure(f1);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on IPSCs from BC, BiC, and AAC onto PYR'];
    subplot(5,3,k);
    [pks, locs] = findpeaks(all_ipsc(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(all_ipsc(:,k),'MinPeakDistance',3000);
    hold on; 
    title (['PYR Number #' num2str(k)])
    xlabel('Time (mses)')
    ylabel('IPSC')
    temp_cur = all_ipsc(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    peaks_all_PV = [peaks_all_PV peaks_all];
end

peaks_all_PV_together = [];
f2 = figure;
for k = 1:1:15
    figure(f2);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on IPSCs from All Inhibitory Cells onto PYR'];
    subplot(5,3,k);
    [pks, locs] = findpeaks(all_ipsc_together(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(all_ipsc_together(:,k),'MinPeakDistance',3000);
    hold on; 
    title (['PYR Number #' num2str(k)])
    xlabel( 'Time (msec)')
    ylabel('IPSC')
    temp_cur_together = all_ipsc_together(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur_together(element,:) = 0;
    end
    peaks_all_together = temp_cur_together;
    peaks_all_PV_together = [peaks_all_PV_together peaks_all_together];
end

%% Mean IPSCs Gathered from AAC, BiC and BC onto PYR - Table and Graph

IPSC_all = [];
ipsc_all = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_all = peaks_all_PV(:,i);
    pks_ipsc_all(pks_ipsc_all == 0) = [];
    ipsc_all_mean = mean(pks_ipsc_all);
    ipsc_all_std = std(pks_ipsc_all);
    ipsc_all = [ipsc_all_mean;ipsc_all_std];
    IPSC_all = [IPSC_all ipsc_all];
end

IPSC_all_table = IPSC_all;
num = (1:15)';
IPSC_all_table = array2table(IPSC_all_table');
IPSC_all_table.num = num;
IPSC_all_table = [IPSC_all_table(:,end) IPSC_all_table(:,1) IPSC_all_table(:,2)];

IPSC_all_table.Properties.VariableNames = {'Pyramidal_Cell_Number', 'Mean_Peak', 'Standard_Deviation'};

IPSC_all_mean = IPSC_all(1,:);
IPSC_all_std = IPSC_all(2,:);
x = linspace(0,15,length(IPSC_all_mean));
figure
scatter(x,IPSC_all_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_all_mean=0.1;
text(x+dx, IPSC_all_mean+dIPSC_all_mean, c);
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_mean,IPSC_all_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC, BiC and AAC onto PYR Cells','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSCs from AAC, BiC, and BC onto PYR Cells
fig = uitable('Data',IPSC_all_table{:,:},...
    'RowName',[],...
    'ColumnName',{'Pyramidal Cell Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% All IPSCs from All Inhibitory Neurons onto PYR - graph and table

IPSC_all_together = [];
ipsc_all_together = [];
for i = 1:1:15% number of PYR cells
    pks_ipsc_all_together = peaks_all_PV_together(:,i);
    pks_ipsc_all_together(pks_ipsc_all_together == 0) = [];
    ipsc_all_mean_together = mean(pks_ipsc_all_together);
    ipsc_all_std_together = std(pks_ipsc_all_together);
    ipsc_all_together = [ipsc_all_mean_together;ipsc_all_std_together];
    IPSC_all_together = [IPSC_all_together ipsc_all_together];
end 

IPSC_all_together_table = IPSC_all_together;
num = (1:15)';
IPSC_all_together_table = array2table(IPSC_all_together_table');
IPSC_all_together_table.num = num;
IPSC_all_together_table = [IPSC_all_together_table(:,end) IPSC_all_together_table(:,1) IPSC_all_together_table(:,2)];

IPSC_all_together_table.Properties.VariableNames = {'Pyramidal_Cell_Number', 'Mean_Peak', 'Standard_Deviation'};


IPSC_all_together_mean = IPSC_all_together(1,:);
IPSC_all_std_together = IPSC_all_together(2,:);
x = linspace(0,14,length(IPSC_all_together_mean));
figure
scatter(x,IPSC_all_together_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_all_together_mean=0.1;
text(x+dx, IPSC_all_together_mean+dIPSC_all_together_mean, c);
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_together_mean,IPSC_all_std_together,'b','LineStyle','none')
title('Mean Peak IPSC from All Inhibitory Cells onto PYRs','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSCs from All Inhibitory Neurons onto PYR Cells
fig = uitable('Data',IPSC_all_together_table{:,:},...
    'RowName',[],...
    'ColumnName',{'Pyramidal Cell Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% IPSCs Only from AAC and BC onto PYR gathered
% Sum all ipsc currents
AAC_BC_ipsc = [];
for i = 1:1:15
    AAC_BC_cur_ipsc = current_AAC(:,i) + current_BC(:,i);
    AAC_BC_ipsc = [AAC_BC_ipsc AAC_BC_cur_ipsc];
end

% Find the peaks of the summed ipsc currents
peaks_AAC_BC = [];
f1 = figure;
for k = 1:1:15
    figure(f1);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on IPSCs from AACs, BCs onto PYRs'];
    subplot(5,3,k);
    [pks, locs] = findpeaks(AAC_BC_ipsc(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(AAC_BC_ipsc(:,k),'MinPeakDistance',3000);
    hold on; 
    title (['PYR Number #' num2str(k)])
    xlabel('Time (1/40 ms)')
    ylabel('IPSC')
    temp_cur = AAC_BC_ipsc(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    pk_AAC_BC = temp_cur;
    peaks_AAC_BC = [peaks_AAC_BC pk_AAC_BC];
end

%% IPSC from AAC and BC onto PYR - graph and table
IPSC_AAC_BC = [];
ipsc_AAC_BC = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_AAC_BC = peaks_AAC_BC(:,i);
    pks_ipsc_AAC_BC(pks_ipsc_AAC_BC == 0) = [];
    ipsc_AAC_BC_mean = mean(pks_ipsc_AAC_BC);
    ipsc_AAC_BC_std = std(pks_ipsc_AAC_BC);
    ipsc_AAC_BC = [ipsc_AAC_BC_mean;ipsc_AAC_BC_std];
    IPSC_AAC_BC = [IPSC_AAC_BC ipsc_AAC_BC];
end

IPSC_AAC_BC_table = IPSC_AAC_BC;
num = (1:15)';
IPSC_AAC_BC_table = array2table(IPSC_AAC_BC_table');
IPSC_AAC_BC_table.num = num;
IPSC_AAC_BC_table = [IPSC_AAC_BC_table(:,end) IPSC_AAC_BC_table(:,1) IPSC_AAC_BC_table(:,2)];

IPSC_AAC_BC_table.Properties.VariableNames = {'Pyramidal_Cell_Number', 'Mean_Peak', 'Standard_Deviation'};

IPSC_AAC_BC_mean = IPSC_AAC_BC(1,:);
IPSC_AAC_BC_std = IPSC_AAC_BC(2,:);
x = linspace(1,15,length(IPSC_AAC_BC_mean));
figure
scatter(x(1,:),IPSC_AAC_BC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_AAC_BC_mean=0.1;
text(x+dx, IPSC_AAC_BC_mean+dIPSC_AAC_BC_mean, c);
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_mean,IPSC_all_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC and AAC onto PYR Cells','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSCs Only from AAC and BC onto PYR Cells
fig = uitable('Data',IPSC_AAC_BC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'Pyramidal Cell Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);

%% IPSCs Only from BiC and BC onto PYR gathered
% Sum all ipsc currents
BiC_BC_ipsc = [];
for i = 1:1:15
    BiC_BC_cur_ipsc = current_BiC(:,i) + current_BC(:,i);
    BiC_BC_ipsc = [BiC_BC_ipsc BiC_BC_cur_ipsc];
end

% Find the peaks of the summed ipsc currents
peaks_BiC_BC = [];
f1 = figure;
for k = 1:1:15
    figure(f1);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on IPSCs from BiC and BC onto PYR'];
    subplot(5,3,k);
    [pks, locs] = findpeaks(BiC_BC_ipsc(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(BiC_BC_ipsc(:,k),'MinPeakDistance',3000);
    hold on; 
    title (['PYR Number #' num2str(k)])
    xlabel('Time (msec)')
    ylabel('IPSC')
    temp_cur = BiC_BC_ipsc(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    pk_BiC_BC = temp_cur;
    peaks_BiC_BC = [peaks_BiC_BC pk_BiC_BC];
end

%% IPSC from BiC and BC onto PYR - graph and table
IPSC_BiC_BC = [];
ipsc_BiC_BC = [];
for i = 1:1:15 % number of PYR cells
    pks_ipsc_BiC_BC = peaks_BiC_BC(:,i);
    pks_ipsc_BiC_BC(pks_ipsc_BiC_BC == 0) = [];
    ipsc_BiC_BC_mean = mean(pks_ipsc_BiC_BC);
    ipsc_BiC_BC_std = std(pks_ipsc_BiC_BC);
    ipsc_BiC_BC = [ipsc_BiC_BC_mean;ipsc_BiC_BC_std];
    IPSC_BiC_BC = [IPSC_BiC_BC ipsc_BiC_BC];
end

IPSC_BiC_BC_table = IPSC_BiC_BC;
num = (1:15)';
IPSC_BiC_BC_table = array2table(IPSC_BiC_BC_table');
IPSC_BiC_BC_table.num = num;
IPSC_BiC_BC_table = [IPSC_BiC_BC_table(:,end) IPSC_BiC_BC_table(:,1) IPSC_BiC_BC_table(:,2)];

IPSC_BiC_BC_table.Properties.VariableNames = {'Pyramidal_Cell_Number', 'Mean_Peak', 'Standard_Deviation'};

IPSC_BiC_BC_mean = IPSC_BiC_BC(1,:);
IPSC_BiC_BC_std = IPSC_BiC_BC(2,:);
x = linspace(1,15,length(IPSC_BiC_BC_mean));
figure
scatter(x(1,:),IPSC_BiC_BC_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dIPSC_BiC_BC_mean=0.1;
text(x+dx, IPSC_BiC_BC_mean+dIPSC_BiC_BC_mean, c);
xlabel('Individual PYR Cells','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak IPSC','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,IPSC_all_mean,IPSC_all_std,'b','LineStyle','none')
title('Mean Peak IPSC from BC and BiC onto PYR Cells','FontSize',15,'FontWeight','bold')

%% Mean Peak and Standard Deviation of IPSCs Only from BiC and BC onto PYR Cells
fig = uitable('Data',IPSC_BiC_BC_table{:,:},...
    'RowName',[],...
    'ColumnName',{'Pyramidal Cell Number','Mean Peak','Standard Deviation'},...
    'Units','Normalized',...
    'Position',[0, 0, 1, 1]);
%% Calculate Mean EPSCs Gathered from PYR and CA3, EC onto PYR
% Sum all ipsc currents
all_epsc = [];
for i = 1:1:15
    tot_cur_epsc =  current_PYR(:,i) + current_ca3(:,i) + current_ec(:,i);  
    all_epsc = [all_epsc tot_cur_epsc];
end

% Find the peaks of the summed epsc currents - only from BC and BiC
epsc_peaks_all = [];
f1 = figure;
for k = 1:1:15
    figure(f1);
    t = annotation('textbox','FontSize',18,'FontWeight','bold'); % This declares the textbox and the options for the character. 
    t.Position = [0 0 1 1]; % (0,0) is the point of the bottom-left corner of the textbox,
    t.HorizontalAlignment = 'center'; % This places the title in the center of the textbox horizontally
    t.VerticalAlignment = 'top'; % This places the title in the top of the textbox vertically
    t.String = ['Peak Detection on EPSCs from PYR, EC, and CA3 on PYR'];
    subplot(5,3,k);
    [pks, locs] = findpeaks(-all_epsc(:,k),'MinPeakDistance',3000); % peak detection
    findpeaks(-all_epsc(:,k),'MinPeakDistance',3000);
    hold on; 
    title (['PYR Number #' num2str(k)])
    xlabel('Time (msec)')
    ylabel('EPSC')
    temp_cur = all_epsc(:,k);
    allrows = (1:40000)';
    notpeak = setdiff(allrows,locs);
    for t = 1:1:numel(notpeak)
        element = notpeak(t,:);
        temp_cur(element,:) = 0;
    end
    peaks_all = temp_cur;
    epsc_peaks_all = [epsc_peaks_all peaks_all];
end

%% All EPSCs from All Excitatory Neurons onto PYR - graph and table

EPSC_all_together = [];
epsc_all_together = [];
for i = 1:1:15
    pks_epsc_all_together = epsc_peaks_all(:,i);
    pks_epsc_all_together(pks_epsc_all_together == 0) = [];
    epsc_all_mean_together = mean(pks_epsc_all_together);
    epsc_all_std_together = std(pks_epsc_all_together);
    epsc_all_together = [epsc_all_mean_together;epsc_all_std_together];
    EPSC_all_together = [EPSC_all_together epsc_all_together];
end 

EPSC_all_together_table = EPSC_all_together;
num = (1:15)';
EPSC_all_together_table = array2table(EPSC_all_together_table');
EPSC_all_together_table.num = num;
EPSC_all_together_table = [EPSC_all_together_table(:,end) EPSC_all_together_table(:,1) EPSC_all_together_table(:,2)];

EPSC_all_together_table.Properties.VariableNames = {'PYR_Number', 'Mean_Peak', 'Standard_Deviation'};
 
EPSC_all_together_mean = EPSC_all_together(1,:);
EPSC_all_std_together = EPSC_all_together(2,:);
x = linspace(0,15,length(EPSC_all_together_mean));
figure
scatter(x,EPSC_all_together_mean,'black','filled');
set(gca, 'XTickLabel',[]);
a = [1:15]'; b =num2str(a); c=cellstr(b);
dx=0.1; dEPSC_all_together_mean=0.1;
text(x+dx, EPSC_all_together_mean+dEPSC_all_together_mean, c);
xlabel('Individual PYRs','FontSize',13,'FontWeight','bold');
ylabel('Mean Peak EPSCs','FontSize',13,'FontWeight','bold');
hold on;
errorbar(x,EPSC_all_together_mean,EPSC_all_std_together,'b','LineStyle','none')
title('Mean Peak EPSCs from All Excitatory Cells onto BiCs','FontSize',15,'FontWeight','bold')

%% Excitatory/Inhibitory Ratios on PYR Cells without CA3/EC
Ratios_PYR = [];
E_I_AAC = abs(EPSC(1,:)./IPSC_AAC(1,:))';
E_I_BC = abs(EPSC(1,:)./IPSC_BC(1,:))';
E_I_BiC = abs(EPSC(1,:)./IPSC_BiC(1,:))';
E_I_AAC_BC = abs(EPSC(1,:)./IPSC_AAC_BC(1,:))';
E_I_BiC_BC = abs(EPSC(1,:)./IPSC_BiC_BC(1,:))';

E_I_all = abs(EPSC(1,:)./IPSC_all(1,:))';
E_I_all_together = abs(EPSC(1,:)./IPSC_all_together(1,:))';

%% Excitatory/Inhibitory Ratios on PYR Cells with CA3/EC
E_I_AAC_with_ca3 = abs(EPSC_all_together(1,:)./IPSC_AAC(1,:))';
E_I_BC_with_ca3 = abs(EPSC_all_together(1,:)./IPSC_BC(1,:))';
E_I_BiC_with_ca3 = abs(EPSC_all_together(1,:)./IPSC_BiC(1,:))';
E_I_AAC_BC_with_ca3 = abs(EPSC_all_together(1,:)./IPSC_AAC_BC(1,:))';
E_I_BiC_BC_with_ca3 = abs(EPSC_all_together(1,:)./IPSC_BiC_BC(1,:))';

E_I_all_with_ca3 = abs(EPSC_all_together(1,:)./IPSC_all(1,:))';
E_I_all_together_with_ca3 = abs(EPSC_all_together(1,:)./IPSC_all_together(1,:))';


%% E/I Ratio - Table without CA3/EC
pyr = 1:15;
Ratios_PYR = [pyr' E_I_AAC E_I_BC E_I_BiC E_I_AAC_BC E_I_BiC_BC E_I_all E_I_all_together];
Ratios_PYR = array2table(Ratios_PYR);
Ratios_PYR.Properties.VariableNames = {'pyr_no' 'Ratio_AAC_on_PYR'...
    'Ratio_BC_on_PYR' 'Ratio_BiC_on_PYR' ...
    'Ratio_AAC_BC_on_PYR'...
    'Ratio_BiC_BC_on_PYR'...
    'Ratio_AAC_BiC_BC_on_PYR' 'All_ipsc_onto_PYR'};
%% E/I Ratio - Table with CA3/EC
pyr = 1:15;
Ratios_PYR_with_ca3 = [pyr' E_I_AAC_with_ca3 E_I_BC_with_ca3 ...
    E_I_BiC_with_ca3 E_I_AAC_BC_with_ca3 E_I_BiC_BC_with_ca3 E_I_all_with_ca3...
    E_I_all_together_with_ca3];
Ratios_PYR_with_ca3 = array2table(Ratios_PYR_with_ca3);

Ratios_PYR_with_ca3.Properties.VariableNames = {'pyr_no' 'Ratio_AAC_on_PYR'...
    'Ratio_BC_on_PYR' 'Ratio_BiC_on_PYR' ...
    'Ratio_AAC_BC_on_PYR'...
    'Ratio_BiC_BC_on_PYR'...
    'Ratio_AAC_BiC_BC_on_PYR' 'All_ipsc_onto_PYR'};
%% Display the E/I table as a figure
% without CA3
uitable('Data',Ratios_PYR{:,:},...
    'RowName', [],...
    'ColumnName',{'Pyramidal Cell Number',...
    'AAC to PYR',...
    'BC to PYR',...
    'BiC to PYR',...
    'AAC, BC onto PYR',...
    'BiC, BC onto PYR',...
    'AAC, BiC and BC to PYR',...
    'All Inhibitory Neurons to PYR'},...
    'Units', 'Normalized',...
    'Position',[0, 0, 1, 1]);

%% With CA3/EC
uitable('Data',Ratios_PYR_with_ca3{:,:},...
    'RowName', [],...
    'ColumnName',{'Pyramidal Cell Number',...
    'AAC to PYR',...
    'BC to PYR',...
    'BiC to PYR',...
    'AAC, BC onto PYR',...
    'BiC, BC onto PYR',...
    'AAC, BiC and BC to PYR',...
    'All Inhibitory Neurons to PYR'},...
    'Units', 'Normalized',...
    'Position',[0, 0, 1, 1]);

%% Read in the files for voltage recordings

g = fullfile('~',... # Change this path accordingly 
    'CA1_SimpleDetailed','Excitatory_Inhibitory_Ratios',...
    'Network_Clamp_Results','pyr',{...
    'pyr_29097_1000';...
    'pyr_36884_1000';...
    'pyr_52458_1000';...
    'pyr_68032_1000';...
    'pyr_83606_1000';...
    'pyr_99180_1000';...
    'pyr_106967_1000';...
    'pyr_114754_1000';...
    'pyr_75819_1000';...
    'pyr_44671_1000';...
    'pyr_60245_1000';...
    'pyr_21310_1000';...
    'pyr_91393_1000';...
    'pyr_301642_1000';...
    'pyr_325003_1000'...
    },{...
    'mytrace_29097_soma.dat';...
    'mytrace_36884_soma.dat';...
    'mytrace_52458_soma.dat';...
    'mytrace_68032_soma.dat';...
    'mytrace_83606_soma.dat';...
    'mytrace_99180_soma.dat';...
    'mytrace_106967_soma.dat';...
    'mytrace_114754_soma.dat';...
    'mytrace_75819_soma.dat';...
    'mytrace_44671_soma.dat';...
    'mytrace_60245_soma.dat';...
    'mytrace_21310_soma.dat';...
    'mytrace_91393_soma.dat';...
    'mytrace_301642_soma.dat';...
    'mytrace_325003_soma.dat'...
    });


%%

allvol= [];
for m = 1:15
    temp_vol = readtable(g{m},'Delimiter','\t');
    temp_vol = table2array(temp_vol);
    allvol = [allvol temp_vol(:,2)];
end

vol = mat2cell(allvol, 40000,...
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]);

%% Graph of Voltage of PYR
figure 
for i = 1:1:15
    temp = vol{i};
    %temp(temp == 0) = [];   %get rid of zeros
    subplot(5,3,i)
    plot(temp);
    hold on
    title (['PYR Number #' num2str(i)]) 
    xlabel('Time (msec)')
    ylabel('Voltage')
end 

%%








