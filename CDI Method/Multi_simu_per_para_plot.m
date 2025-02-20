clear,clc
cd('....\MATLAB Code\CDI Variablen') %SET PATH WHERE YOU SAVED THE SIMULATION DATASETS
addpath('....\MATLAB Code') %ADD PATH WHERE THE CDI CODE IS AT

output_path = '....\MATLAB Code\CDI Variablen\Multiple EPK per set CSV'; %DEFINE OUTPUT PATH
    


x_min = 0.92;  % Lower limit for x-axis
x_max = 1.06;  % Upper limit for x-axis
y_min = 0;
y_max = 4;

% FUTURE JUMP
% datasets = {
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2021_BVola_15_New_30_2015_2025_001.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2021_BVola_15_New_30_2015_2025_002.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2021_BVola_15_New_30_2015_2025_003.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2022_BVola_15_New_30_2015_2022_001.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2022_BVola_15_New_30_2015_2025_002.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2022_BVola_15_New_30_2015_2025_003.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2022_BVola_15_New_33_2015_2025_001.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2022_BVola_15_New_33_2015_2025_002.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2022_BVola_15_New_33_2015_2025_003.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2022_BVola_15_New_33_2015_2025_004.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2022_BVola_15_New_33_2015_2025_005.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2022_BVola_15_New_33_2015_2025_006.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2023_BVola_15_New_30_2015_2025_001.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2023_BVola_15_New_30_2015_2025_002.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2023_BVola_15_New_30_2015_2025_003.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2024_BVola_15_New_30_2015_2025_001.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2024_BVola_15_New_30_2015_2025_002.mat',
%    'filtered_cdi_2020-01-02_Future_Jump_01_01_2024_BVola_15_New_30_2015_2025_003.mat',
%    'filtered_cdi_2020-01-02_Past_Jump_01_01_2010_BVola_30_New_15_2000_2024_003.mat',
%    'filtered_cdi_2020-01-02_Past_Jump_01_01_2018_BVola_30_New_15_2015_2025_001.mat',
%    'filtered_cdi_2020-01-02_Past_Jump_01_01_2018_BVola_30_New_15_2015_2025_002.mat',
%    'filtered_cdi_2020-01-02_Past_Jump_01_01_2018_BVola_30_New_15_2015_2025_003.mat'
% };

%Future Jump Long Case
% datasets = {
%     'filtered_cdi_2000-01-02_Future_Jump_01_01_2013_BVola_15_New_30_2000_2025_Date_2000_001',
%     'filtered_cdi_2000-01-02_Future_Jump_01_01_2013_BVola_15_New_30_2000_2025_Date_2000_002',
%     'filtered_cdi_2000-01-02_Future_Jump_01_01_2013_BVola_15_New_30_2000_2025_Date_2000_003',
%     'filtered_cdi_2000-01-02_Future_Jump_01_01_2013_BVola_15_New_30_2000_2025_Date_2000_004',
%     'filtered_cdi_2000-01-02_Future_Jump_01_01_2013_BVola_15_New_30_2000_2025_Date_2000_005',
%     'filtered_cdi_2000-01-02_Future_Jump_01_01_2013_BVola_15_New_30_2000_2025_Date_2000_006'
% };

%Constant 
% datasets = {
%    'filtered_cdi_2020-01-02_Constant_BVola_15_2015_2025_001.mat',
%    'filtered_cdi_2020-01-02_Constant_BVola_15_2015_2025_002.mat',
%    'filtered_cdi_2020-01-02_Constant_BVola_15_2015_2025_003.mat',
%    'filtered_cdi_2020-01-02_Constant_BVola_15_2015_2025_004.mat',
%    'filtered_cdi_2020-01-02_Constant_BVola_15_2015_2025_005.mat',
%    'filtered_cdi_2020-01-02_Constant_BVola_15_2015_2025_006.mat',
%    'filtered_cdi_2020-01-02_Constant_Vola_2000_2040_002.mat',
%    'filtered_cdi_2020-01-02_Constant_Vola_2000_2040_003.mat'
% };
%past
datasets = {
   'filtered_cdi_2020-01-02_Past_Jump_01_01_2010_BVola_30_New_15_2000_2024_003.mat',
   'filtered_cdi_2020-01-02_Past_Jump_01_01_2018_BVola_30_New_15_2015_2025_001.mat',
   'filtered_cdi_2020-01-02_Past_Jump_01_01_2018_BVola_30_New_15_2015_2025_002.mat',
   'filtered_cdi_2020-01-02_Past_Jump_01_01_2018_BVola_30_New_15_2015_2025_003.mat',
   'filtered_cdi_2000-01-02_Past_Jump_01_01_2010_BVola_30_New_15_2000_2024_001.mat',
   'filtered_cdi_2000-01-02_Past_Jump_01_01_2010_BVola_30_New_15_2000_2024_002.mat',
   'filtered_cdi_2000-01-02_Past_Jump_01_01_2010_BVola_30_New_15_2000_2024_003.mat'
};

results_4 = cell(length(datasets), 2); % Für moneyness und sampleestimate
results_5 = cell(length(datasets), 2);
results_6 = cell(length(datasets), 2);

for i = 1:length(datasets)
    data = load(datasets{i});
    
    realizedKhRet = cell(size(data.realizedKhRet, 1), 1);
    realizedQdenRet = cell(size(data.realizedQdenRet, 1), 1);
    
    for j = 1:size(data.realizedKhRet, 1)
        realizedKhRet{j} = data.realizedKhRet(j,:);
        realizedQdenRet{j} = data.realizedQdenRet(j,:);
    end
    
    
    % CDI Estimation
    [sampleestimate_4_4, returns_4_4] = CDI_estimator(realizedKhRet, realizedQdenRet, @OptSDF, 4, 4, 0.001);
    [sampleestimate_5_5, returns_5_5] = CDI_estimator(realizedKhRet, realizedQdenRet, @OptSDF, 5, 5, 0.001);
    [sampleestimate_6_6, returns_6_6] = CDI_estimator(realizedKhRet, realizedQdenRet, @OptSDF, 6, 6, 0.001);
    [sampleestimate_7_7, returns_7_7] = CDI_estimator(realizedKhRet, realizedQdenRet, @OptSDF, 7, 7, 0.001);
    [sampleestimate_8_8, returns_8_8] = CDI_estimator(realizedKhRet, realizedQdenRet, @OptSDF, 8, 8, 0.001);
    [sampleestimate_9_9, returns_9_9] = CDI_estimator(realizedKhRet, realizedQdenRet, @OptSDF, 9, 9, 0.001);
    
    returns_4_4 = returns_4_4(:);
    sampleestimate_4_4 = sampleestimate_4_4(:);
    returns_5_5 = returns_5_5(:);
    sampleestimate_5_5 = sampleestimate_5_5(:);
    returns_6_6 = returns_6_6(:);
    sampleestimate_6_6 = sampleestimate_6_6(:);
    returns_7_7 = returns_7_7(:);
    sampleestimate_7_7 = sampleestimate_7_7(:);
    returns_8_8 = returns_8_8(:);
    sampleestimate_8_8 = sampleestimate_8_8(:);
    returns_9_9 = returns_9_9(:);
    sampleestimate_9_9 = sampleestimate_9_9(:);
    
    % safe as
    results_4{i,1} = exp(returns_4_4);
    results_4{i,2} = sampleestimate_4_4;
    results_5{i,1} = exp(returns_5_5);
    results_5{i,2} = sampleestimate_5_5;
    results_6{i,1} = exp(returns_6_6);
    results_6{i,2} = sampleestimate_6_6;
    results_7{i,1} = exp(returns_7_7);
    results_7{i,2} = sampleestimate_7_7;
    results_8{i,1} = exp(returns_8_8);
    results_8{i,2} = sampleestimate_8_8;
    results_9{i,1} = exp(returns_9_9);
    results_9{i,2} = sampleestimate_9_9;
end

%PLOTTING

figure('Position', [100 100 1200 1600]);

subplot(3,2,1)
colors = lines(length(datasets));
for i = 1:length(datasets)
    plot(results_4{i,1}, results_4{i,2}, 'LineWidth', 2, 'Color', colors(i,:))
    hold on
end
title('4 Parameter Sets')
xlabel('Moneyness')
ylabel('EPK')
%legend(['Dataset ' num2str(1:length(datasets))])
xlim([x_min x_max])
ylim([y_min y_max])
grid off

subplot(3,2,2)
for i = 1:length(datasets)
    plot(results_5{i,1}, results_5{i,2}, 'LineWidth', 2, 'Color', colors(i,:))
    hold on
end
title('5 Parameter Sets')
xlabel('Moneyness')
ylabel('EPK')
xlim([x_min x_max])
ylim([y_min y_max])
grid off

subplot(3,2,3)
for i = 1:length(datasets)
    plot(results_6{i,1}, results_6{i,2}, 'LineWidth', 2, 'Color', colors(i,:))
    hold on
end
title('6 Parameter Sets')
xlabel('Moneyness')
ylabel('EPK')
xlim([x_min x_max])
ylim([y_min y_max])
grid off

% 7er Parameter
subplot(3,2,4)
for i = 1:length(datasets)
    plot(results_7{i,1}, results_7{i,2}, 'LineWidth', 2, 'Color', colors(i,:))
    hold on
end
title('7 Parameter Sets')
xlabel('Moneyness')
ylabel('EPK')
xlim([x_min x_max])
ylim([y_min y_max])
grid off

% 8er Parameter
subplot(3,2,5)
for i = 1:length(datasets)
    plot(results_8{i,1}, results_8{i,2}, 'LineWidth', 2, 'Color', colors(i,:))
    hold on
end
title('8 Parameter Sets')
xlabel('Moneyness')
ylabel('EPK')
xlim([x_min x_max])
ylim([y_min y_max])
grid off

% 9er Parameter
subplot(3,2,6)
for i = 1:length(datasets)
    plot(results_9{i,1}, results_9{i,2}, 'LineWidth', 2, 'Color', colors(i,:))
    hold on
end
title('9 Parameter Sets')
xlabel('Moneyness')
ylabel('EPK')
xlim([x_min x_max])
ylim([y_min y_max])
grid off

% safw pots
saveas(gcf, 'EPK_all_parameter_sets.png')



for i = 1:length(datasets)
    % Combine all results for this dataset into one matrix
    data_4 = [results_4{i,1}, results_4{i,2}];
    data_5 = [results_5{i,1}, results_5{i,2}];
    data_6 = [results_6{i,1}, results_6{i,2}];
    data_7 = [results_7{i,1}, results_7{i,2}];
    data_8 = [results_8{i,1}, results_8{i,2}];
    data_9 = [results_9{i,1}, results_9{i,2}];
    
    % Save each dataset to a CSV file
    writematrix(data_4, fullfile(output_path, sprintf('epk_data_4_4_dataset_%d.csv', i)));
    writematrix(data_5, fullfile(output_path, sprintf('epk_data_5_5_dataset_%d.csv', i)));
    writematrix(data_6, fullfile(output_path, sprintf('epk_data_6_6_dataset_%d.csv', i)));
    writematrix(data_7, fullfile(output_path, sprintf('epk_data_7_7_dataset_%d.csv', i)));
    writematrix(data_8, fullfile(output_path, sprintf('epk_data_8_8_dataset_%d.csv', i)));
    writematrix(data_9, fullfile(output_path, sprintf('epk_data_9_9_dataset_%d.csv', i)));
end



