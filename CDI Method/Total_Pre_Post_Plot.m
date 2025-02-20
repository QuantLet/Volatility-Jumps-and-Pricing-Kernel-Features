
clear,clc
cd('....\MATLAB Code\CDI Variablen') %SET PATH WHERE YOU SAVED THE SIMULATION DATASETS
addpath('....\MATLAB Code') %ADD PATH WHERE THE CDI CODE IS AT

function [moneyness_all, sampleestimate_all] = process_data_all_kernels(data)
    realizedKhRet = cell(size(data.realizedKhRet, 1), 1);
    realizedQdenRet = cell(size(data.realizedQdenRet, 1), 1);
    
    for i = 1:size(data.realizedKhRet, 1)
        realizedKhRet{i} = data.realizedKhRet(i,:);
        realizedQdenRet{i} = data.realizedQdenRet(i,:);
    end
    
    [sampleestimate_4_4, returns_4_4] = CDI_estimator(realizedKhRet, realizedQdenRet, @OptSDF, 4, 4, 0.001);
    [sampleestimate_5_5, returns_5_5] = CDI_estimator(realizedKhRet, realizedQdenRet, @OptSDF, 5, 5, 0.001);
    [sampleestimate_6_6, returns_6_6] = CDI_estimator(realizedKhRet, realizedQdenRet, @OptSDF, 6, 6, 0.001);
    [sampleestimate_7_7, returns_7_7] = CDI_estimator(realizedKhRet, realizedQdenRet, @OptSDF, 7, 7, 0.001);
    [sampleestimate_8_8, returns_8_8] = CDI_estimator(realizedKhRet, realizedQdenRet, @OptSDF, 8, 8, 0.001);
    [sampleestimate_9_9, returns_9_9] = CDI_estimator(realizedKhRet, realizedQdenRet, @OptSDF, 9, 9, 0.001);
    
    moneyness_all = {exp(returns_4_4(:)), exp(returns_5_5(:)), exp(returns_6_6(:)), ...
                    exp(returns_7_7(:)), exp(returns_8_8(:)), exp(returns_9_9(:))};
    sampleestimate_all = {sampleestimate_4_4(:), sampleestimate_5_5(:), sampleestimate_6_6(:), ...
                         sampleestimate_7_7(:), sampleestimate_8_8(:), sampleestimate_9_9(:)};
end

function save_data_to_csv(moneyness, sampleestimate, prefix)

    output_path = '....\MATLAB Code\CDI Variablen\TotalPostPre CSV'; % SET OUTPUT PATH
    
    if ~exist(output_path, 'dir') %create if not there
        mkdir(output_path);
    end

    for i = 1:length(moneyness)
        data = [moneyness{i}, sampleestimate{i}];
        filename = fullfile(output_path, sprintf('%s_epk_data_%d_%d.csv', prefix, i+3, i+3));
        writematrix(data, filename);
    end
end

data_total = load('filtered_cdi_2000-01-02_Future_Jump_01_01_2013_BVola_15_New_30_2000_2025_Date_2000_006');
data_pre = load('filtered_cdi_2000-01-01_to_2012-12-30_Future_Jump_01_01_2013_BVola_15_New_30_2000_2013_Date_2000_006_PRE_JUMP');
data_post = load('filtered_cdi_2013-01-01_Future_Jump_01_01_2013_BVola_15_New_30_2000_2013_Date_2000_006_POST_JUMP');

[moneyness_total, sampleestimate_total] = process_data_all_kernels(data_total);
[moneyness_pre, sampleestimate_pre] = process_data_all_kernels(data_pre);
[moneyness_post, sampleestimate_post] = process_data_all_kernels(data_post);

% Csv
save_data_to_csv(moneyness_total, sampleestimate_total, 'total');
save_data_to_csv(moneyness_pre, sampleestimate_pre, 'pre');
save_data_to_csv(moneyness_post, sampleestimate_post, 'post');

figure('Position', [100 100 1200 800]);
set(gcf, 'Color', 'none'); 


colors = {'b-', 'r-', 'g-', 'k-', 'm-', 'c-'};
legend_labels = {'4 bases/moments', '5 bases/moments', '6 bases/moments', ...
                '7 bases/moments', '8 bases/moments', '9 bases/moments'};

% TOTAL case
subplot(2,1,1)
for i = 1:6
    plot(moneyness_total{i}, sampleestimate_total{i}, colors{i}, 'LineWidth', 2)
    hold on
end
xlabel('Moneyness')
ylabel('EPK')
title('Estimated Pricing Kernel - Total Period')
grid off
set(gca, 'Color', 'none')
xlim([0.92 1.06])
ylim([0 4])

% PRE JUMP
subplot(2,2,3)
for i = 1:6
    plot(moneyness_pre{i}, sampleestimate_pre{i}, colors{i}, 'LineWidth', 2)
    hold on
end
xlabel('Moneyness')
ylabel('EPK')
title('EPK - Pre Jump Period')
grid off
set(gca, 'Color', 'none')
xlim([0.92 1.06])
ylim([0 4])

% POST JUMP
subplot(2,2,4)
for i = 1:6
    plot(moneyness_post{i}, sampleestimate_post{i}, colors{i}, 'LineWidth', 2)
    hold on
end
xlabel('Moneyness')
ylabel('EPK')
title('EPK - Post Jump Period')
grid off
set(gca, 'Color', 'none')
xlim([0.92 1.06])
ylim([0 4])

% Save
set(gcf, 'InvertHardcopy', 'off'); 
savePath = fullfile('....\MATLAB Code\CDI Variablen\TotalPostPre CSV', 'CDI_EPK_Total_Pre_Post.png'); % OUTPUT PLOT SAVE AS
print(gcf, savePath, '-dpng', '-r300');