clear,clc
cd('......MATLAB Code\CDI Variablen') %SET PATH WHERE YOU SAVED THE SIMULATION DATASETS
addpath('......\MATLAB Code') %ADD PATH WHERE THE CDI CODE IS AT

output_path = '\MATLAB Code\CDI Variablen\EPK per param set CSV\'; %DEFINE OUTPUT PATH


%---

x_min = 0.92;  
x_max = 1.06;  

y_min = 0;
y_max = 4;

%---

data = load("filtered_cdi_2020-01-01_Constant_20_2015_2020_Paper.mat")
%---

realizedKhRet = cell(size(data.realizedKhRet, 1), 1);
realizedQdenRet = cell(size(data.realizedQdenRet, 1), 1);

for i = 1:size(data.realizedKhRet, 1)
   realizedKhRet{i} = data.realizedKhRet(i,:); %   realizedKhRet{i} = data.realizedKhRet(i,:)';
   realizedQdenRet{i} = data.realizedQdenRet(i,:);
end



%Estimate EPK
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

% comvert to moneyness
moneyness_4_4 = exp(returns_4_4);
moneyness_5_5 = exp(returns_5_5);
moneyness_6_6 = exp(returns_6_6);
moneyness_7_7 = exp(returns_7_7);
moneyness_8_8 = exp(returns_8_8);
moneyness_9_9 = exp(returns_9_9);


if ~exist(output_path, 'dir')
    mkdir(output_path);
end

% save as .csv
data_4_4 = [moneyness_4_4, sampleestimate_4_4];
data_5_5 = [moneyness_5_5, sampleestimate_5_5];
data_6_6 = [moneyness_6_6, sampleestimate_6_6];
data_7_7 = [moneyness_7_7, sampleestimate_7_7];
data_8_8 = [moneyness_8_8, sampleestimate_8_8];
data_9_9 = [moneyness_9_9, sampleestimate_9_9];

writematrix(data_4_4, fullfile(output_path, 'epk_data_4_4.csv'));
writematrix(data_5_5, fullfile(output_path, 'epk_data_5_5.csv'));
writematrix(data_6_6, fullfile(output_path, 'epk_data_6_6.csv'));
writematrix(data_7_7, fullfile(output_path, 'epk_data_7_7.csv'));
writematrix(data_8_8, fullfile(output_path, 'epk_data_8_8.csv'));
writematrix(data_9_9, fullfile(output_path, 'epk_data_9_9.csv'));

%Plot EPK
figure('Position', [100 100 1000 600]);

subplot(1,1,1)
plot(moneyness_4_4, sampleestimate_4_4, 'b-', 'LineWidth', 2)
hold on
plot(moneyness_5_5, sampleestimate_5_5, 'r-', 'LineWidth', 2)
plot(moneyness_6_6, sampleestimate_6_6, 'g-', 'LineWidth', 2)
plot(moneyness_7_7, sampleestimate_7_7, 'k-', 'LineWidth', 2)
plot(moneyness_8_8, sampleestimate_8_8, 'm-', 'LineWidth', 2)
plot(moneyness_9_9, sampleestimate_9_9, 'c-', 'LineWidth', 2)

xlabel('Moneyness')
ylabel('EPK')
title('Estimated Pricing Kernel - All Estimates')
legend('4 bases/moments', '5 bases/moments', '6 bases/moments', '7 bases/moments', '8 bases/moments', '9 bases/moments')
grid off

xlim([x_min, x_max])
ylim([y_min, y_max])

 
all_estimates = [sampleestimate_4_4; sampleestimate_5_5; sampleestimate_6_6; sampleestimate_7_7];
valid_estimates = all_estimates(~isnan(all_estimates) & ~isinf(all_estimates));
median_est = median(valid_estimates);
std_est = std(valid_estimates);

saveas(gcf, 'CDI EPK all parameter sets.png')

