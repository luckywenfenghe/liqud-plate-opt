%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           STEP SIZE GROWTH RATE COMPARISON              %
%    Comparing 1.08 vs 1.04 growth rates over iterations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

fprintf('=========================================================\n');
fprintf('      STEP SIZE GROWTH RATE COMPARISON\n');
fprintf('=========================================================\n');

% Parameters
mv_init = 0.01;
mv_max = 0.05;
warmup_iterations = 15;
iterations = 1:warmup_iterations;

% Growth rates
old_rate = 1.08;  % Original rate
new_rate = 1.04;  % New reduced rate

% Calculate step sizes
old_steps = mv_init * (old_rate .^ iterations);
new_steps = mv_init * (new_rate .^ iterations);

% Apply cap at mv_max * 0.5
cap = mv_max * 0.5;
old_steps_capped = min(old_steps, cap);
new_steps_capped = min(new_steps, cap);

% Display comparison table
fprintf('Iteration | Old Rate (1.08) | New Rate (1.04) | Difference\n');
fprintf('----------|------------------|------------------|----------\n');
for i = 1:warmup_iterations
    diff_percent = (old_steps_capped(i) - new_steps_capped(i)) / old_steps_capped(i) * 100;
    fprintf('%9d | %15.4f | %15.4f | %7.1f%%\n', ...
        i, old_steps_capped(i), new_steps_capped(i), diff_percent);
end

% Summary statistics
fprintf('\n=========================================================\n');
fprintf('SUMMARY:\n');
fprintf('Final step size (iteration 15):\n');
fprintf('  Old rate (1.08): %.4f\n', old_steps_capped(end));
fprintf('  New rate (1.04): %.4f\n', new_steps_capped(end));
fprintf('  Reduction: %.1f%%\n', (old_steps_capped(end) - new_steps_capped(end)) / old_steps_capped(end) * 100);

fprintf('\nAverage step size over warmup:\n');
fprintf('  Old rate (1.08): %.4f\n', mean(old_steps_capped));
fprintf('  New rate (1.04): %.4f\n', mean(new_steps_capped));
fprintf('  Reduction: %.1f%%\n', (mean(old_steps_capped) - mean(new_steps_capped)) / mean(old_steps_capped) * 100);

fprintf('\nMaximum step size reached:\n');
fprintf('  Old rate (1.08): %.4f (iteration %d)\n', max(old_steps_capped), find(old_steps_capped == max(old_steps_capped), 1));
fprintf('  New rate (1.04): %.4f (iteration %d)\n', max(new_steps_capped), find(new_steps_capped == max(new_steps_capped), 1));

% Plot comparison
figure('Position', [100, 100, 800, 500]);
plot(iterations, old_steps_capped, 'r-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Old Rate (1.08)');
hold on;
plot(iterations, new_steps_capped, 'b-s', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'New Rate (1.04)');
yline(cap, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Cap (mv\_max * 0.5)');

xlabel('Iteration');
ylabel('Move Limit (Step Size)');
title('Step Size Growth Rate Comparison');
legend('Location', 'southeast');
grid on;

% Annotate key points
text(8, old_steps_capped(8) + 0.002, sprintf('%.4f', old_steps_capped(8)), 'Color', 'red', 'FontWeight', 'bold');
text(8, new_steps_capped(8) - 0.002, sprintf('%.4f', new_steps_capped(8)), 'Color', 'blue', 'FontWeight', 'bold');

fprintf('\n=========================================================\n');
fprintf('BENEFITS OF REDUCED GROWTH RATE:\n');
fprintf('✓ More gradual step size increase\n');
fprintf('✓ Better stability in early iterations\n');
fprintf('✓ Reduced risk of oscillations\n');
fprintf('✓ More predictable convergence behavior\n');
fprintf('=========================================================\n'); 