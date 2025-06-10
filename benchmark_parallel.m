%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    PARALLEL PERFORMANCE BENCHMARK SCRIPT                %
%    FOR THERMAL TOPOLOGY OPTIMIZATION, 2024              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

fprintf('=======================================================\n');
fprintf('    热传递拓扑优化并行性能基准测试\n');
fprintf('=======================================================\n');

%% 测试配置
test_configs = [
    struct('nely', 30, 'name', '小网格 (30x30)'),
    struct('nely', 60, 'name', '中等网格 (60x60)'),
    struct('nely', 90, 'name', '大网格 (90x90)')
];

% 限制迭代次数用于测试
max_test_iterations = 5;

%% 检查系统信息
fprintf('\n系统信息:\n');
fprintf('MATLAB版本: %s\n', version);
fprintf('可用内存: %.2f GB\n', memory().MemAvailableAllArrays / 1e9);

% 检查Parallel Computing Toolbox
if license('test', 'Distrib_Computing_Toolbox')
    fprintf('Parallel Computing Toolbox: 可用\n');
    
    % 获取并行池信息
    p = gcp('nocreate');
    if isempty(p)
        parpool('local');
        p = gcp;
    end
    fprintf('并行工作器数量: %d\n', p.NumWorkers);
    fprintf('并行池类型: %s\n', class(p));
else
    error('Parallel Computing Toolbox 不可用，无法进行并行计算测试');
end

%% 运行基准测试
results = [];

for i = 1:length(test_configs)
    config = test_configs(i);
    fprintf('\n=======================================================\n');
    fprintf('测试配置: %s\n', config.name);
    fprintf('=======================================================\n');
    
    % 串行性能测试
    fprintf('\n--- 串行版本测试 ---\n');
    serial_time = run_serial_benchmark(config.nely, max_test_iterations);
    
    % 并行性能测试
    fprintf('\n--- 并行版本测试 ---\n');
    parallel_time = run_parallel_benchmark(config.nely, max_test_iterations);
    
    % 计算加速比
    speedup = serial_time / parallel_time;
    efficiency = speedup / p.NumWorkers * 100;
    
    % 保存结果
    results = [results; struct(...
        'config', config.name, ...
        'nely', config.nely, ...
        'serial_time', serial_time, ...
        'parallel_time', parallel_time, ...
        'speedup', speedup, ...
        'efficiency', efficiency)];
    
    % 显示结果
    fprintf('\n结果总结:\n');
    fprintf('串行时间: %.2f 秒\n', serial_time);
    fprintf('并行时间: %.2f 秒\n', parallel_time);
    fprintf('加速比: %.2fx\n', speedup);
    fprintf('并行效率: %.1f%%\n', efficiency);
end

%% 显示完整结果表
fprintf('\n=======================================================\n');
fprintf('                   完整性能结果表\n');
fprintf('=======================================================\n');
fprintf('%-15s %10s %10s %8s %8s\n', '配置', '串行(s)', '并行(s)', '加速比', '效率(%)');
fprintf('-------------------------------------------------------\n');
for i = 1:length(results)
    r = results(i);
    fprintf('%-15s %10.2f %10.2f %8.2fx %7.1f%%\n', ...
        r.config, r.serial_time, r.parallel_time, r.speedup, r.efficiency);
end

%% 可视化结果
figure('Position', [100, 100, 1200, 400]);

% 子图1: 计算时间比较
subplot(1, 3, 1);
grid_sizes = [results.nely];
serial_times = [results.serial_time];
parallel_times = [results.parallel_time];

bar([serial_times; parallel_times]');
set(gca, 'XTickLabel', {results.config});
xlabel('网格配置');
ylabel('计算时间 (秒)');
title('串行 vs 并行计算时间');
legend('串行', '并行', 'Location', 'northwest');
grid on;

% 子图2: 加速比
subplot(1, 3, 2);
speedups = [results.speedup];
bar(speedups);
set(gca, 'XTickLabel', {results.config});
xlabel('网格配置');
ylabel('加速比');
title('并行加速比');
grid on;
% 添加理想加速比线
hold on;
plot([0.5, length(results)+0.5], [p.NumWorkers, p.NumWorkers], 'r--', 'LineWidth', 2);
legend('实际加速比', sprintf('理想加速比 (%dx)', p.NumWorkers), 'Location', 'northwest');

% 子图3: 并行效率
subplot(1, 3, 3);
efficiencies = [results.efficiency];
bar(efficiencies);
set(gca, 'XTickLabel', {results.config});
xlabel('网格配置');
ylabel('并行效率 (%)');
title('并行效率');
grid on;
ylim([0, 100]);
% 添加100%效率线
hold on;
plot([0.5, length(results)+0.5], [100, 100], 'r--', 'LineWidth', 2);
legend('实际效率', '理想效率 (100%)', 'Location', 'northwest');

sgtitle('热传递拓扑优化并行性能基准测试结果');

%% 性能分析和建议
fprintf('\n=======================================================\n');
fprintf('                   性能分析和建议\n');
fprintf('=======================================================\n');

avg_speedup = mean(speedups);
avg_efficiency = mean(efficiencies);

fprintf('平均加速比: %.2fx\n', avg_speedup);
fprintf('平均并行效率: %.1f%%\n', avg_efficiency);

if avg_efficiency > 80
    fprintf('\n✓ 并行性能优秀: 效率超过80%%\n');
elseif avg_efficiency > 60
    fprintf('\n◐ 并行性能良好: 效率在60-80%%之间\n');
    fprintf('建议：考虑优化内存访问模式或增加问题规模\n');
else
    fprintf('\n✗ 并行性能需要改进: 效率低于60%%\n');
    fprintf('建议：检查负载均衡、减少通信开销或优化算法\n');
end

% 网格规模对性能的影响
if length(results) > 1
    fprintf('\n网格规模分析:\n');
    for i = 2:length(results)
        size_ratio = (results(i).nely / results(i-1).nely)^2;
        time_ratio_serial = results(i).serial_time / results(i-1).serial_time;
        time_ratio_parallel = results(i).parallel_time / results(i-1).parallel_time;
        
        fprintf('网格从%dx%d到%dx%d (%.1fx元素):\n', ...
            results(i-1).nely, results(i-1).nely, ...
            results(i).nely, results(i).nely, size_ratio);
        fprintf('  串行时间增长: %.2fx\n', time_ratio_serial);
        fprintf('  并行时间增长: %.2fx\n', time_ratio_parallel);
        fprintf('  复杂度比: %.2f (理想为%.1f)\n', ...
            time_ratio_serial/size_ratio, size_ratio);
    end
end

fprintf('\n基准测试完成！\n');

%% 保存结果
save('benchmark_results.mat', 'results', 'p');
fprintf('结果已保存到 benchmark_results.mat\n');

%% 串行基准测试函数
function elapsed_time = run_serial_benchmark(nely, max_iterations)
    % 设置测试参数
    probtype = 3;
    Lx = 1.0; Ly = 1.0;
    nelx = nely*Lx/Ly;
    neltot = nelx*nely;
    
    % 模拟主要计算循环
    fprintf('  网格大小: %dx%d (%d 元素)\n', nelx, nely, neltot);
    fprintf('  最大迭代: %d\n', max_iterations);
    
    % 初始化数据
    dxv = ones(1,neltot); dyv = ones(1,neltot);
    muv = ones(1,neltot)*1e2; rhov = ones(1,neltot)*1e1;
    kv = ones(1,neltot)*0.8; Cpv = ones(1,neltot)*4180;
    Qv = ones(1,neltot)*1000; alpha = ones(1,neltot)*1e4;
    
    % 模拟状态变量
    uVars = rand(neltot, 8);
    pVars = rand(neltot, 4);
    TVars = rand(neltot, 4);
    
    % 开始计时
    tic;
    
    for iter = 1:max_iterations
        fprintf('    串行迭代 %d/%d\r', iter, max_iterations);
        
        % 模拟残差计算 (串行)
        sR = zeros(16, neltot);
        for i = 1:neltot
            % 模拟计算密集型操作
            temp_calc = sin(i) * cos(dxv(i)) + sqrt(muv(i)) * log(rhov(i)+1);
            sR(:,i) = temp_calc * ones(16,1);
        end
        
        % 模拟雅可比矩阵计算 (串行)
        sJ = zeros(256, neltot);
        for i = 1:neltot
            % 模拟计算密集型操作
            temp_calc = exp(-i/1000) * muv(i) / (rhov(i) + kv(i));
            sJ(:,i) = temp_calc * ones(256,1);
        end
        
        % 模拟目标函数计算 (串行)
        phiVals = zeros(1, neltot);
        for i = 1:neltot
            % 模拟计算密集型操作
            phiVals(i) = sum(uVars(i,:).^2) + sum(TVars(i,:).^2) * kv(i);
        end
        
        % 模拟敏感度计算 (串行)
        drdgVals = zeros(16, neltot);
        for i = 1:neltot
            % 模拟计算密集型操作
            temp_calc = alpha(i) * sqrt(sum(uVars(i,:).^2));
            drdgVals(:,i) = temp_calc * ones(16,1);
        end
    end
    
    elapsed_time = toc;
    fprintf('\n  串行完成时间: %.2f 秒\n', elapsed_time);
end

%% 并行基准测试函数
function elapsed_time = run_parallel_benchmark(nely, max_iterations)
    % 设置测试参数
    probtype = 3;
    Lx = 1.0; Ly = 1.0;
    nelx = nely*Lx/Ly;
    neltot = nelx*nely;
    
    % 模拟主要计算循环
    fprintf('  网格大小: %dx%d (%d 元素)\n', nelx, nely, neltot);
    fprintf('  最大迭代: %d\n', max_iterations);
    
    % 初始化数据
    dxv = ones(1,neltot); dyv = ones(1,neltot);
    muv = ones(1,neltot)*1e2; rhov = ones(1,neltot)*1e1;
    kv = ones(1,neltot)*0.8; Cpv = ones(1,neltot)*4180;
    Qv = ones(1,neltot)*1000; alpha = ones(1,neltot)*1e4;
    
    % 模拟状态变量
    uVars = rand(neltot, 8);
    pVars = rand(neltot, 4);
    TVars = rand(neltot, 4);
    
    % 开始计时
    tic;
    
    for iter = 1:max_iterations
        fprintf('    并行迭代 %d/%d\r', iter, max_iterations);
        
        % 模拟残差计算 (并行)
        sR = zeros(16, neltot);
        parfor i = 1:neltot
            % 模拟计算密集型操作
            temp_calc = sin(i) * cos(dxv(i)) + sqrt(muv(i)) * log(rhov(i)+1);
            sR(:,i) = temp_calc * ones(16,1);
        end
        
        % 模拟雅可比矩阵计算 (并行)
        sJ = zeros(256, neltot);
        parfor i = 1:neltot
            % 模拟计算密集型操作
            temp_calc = exp(-i/1000) * muv(i) / (rhov(i) + kv(i));
            sJ(:,i) = temp_calc * ones(256,1);
        end
        
        % 模拟目标函数计算 (并行)
        phiVals = zeros(1, neltot);
        parfor i = 1:neltot
            % 模拟计算密集型操作
            phiVals(i) = sum(uVars(i,:).^2) + sum(TVars(i,:).^2) * kv(i);
        end
        
        % 模拟敏感度计算 (并行)
        drdgVals = zeros(16, neltot);
        parfor i = 1:neltot
            % 模拟计算密集型操作
            temp_calc = alpha(i) * sqrt(sum(uVars(i,:).^2));
            drdgVals(:,i) = temp_calc * ones(16,1);
        end
    end
    
    elapsed_time = toc;
    fprintf('\n  并行完成时间: %.2f 秒\n', elapsed_time);
end 