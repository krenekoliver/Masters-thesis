% ================================================
% Asymptotic approximation λ_k for Coffey-Evans equation (β = 30)
% Dirichlet: [-π/2, π/2] -> [0, π]
% ================================================

clear; clc; close all;

k = 1:50;

n = k + 1;                    

P1 = 450 * pi;
P2 = 305550 * pi;

% Asymptotic approximation
sqrt_lambda = n + (P1 ./ (2*pi*n)) + (1./(4*pi*n.^3)) .* (P2/2 - P1.^2 / pi);
lambda_asymp = sqrt_lambda .^ 2;

% Numeric values from Matslise
x = [5 10 20 30 40 50];
y = [340.888, 637.682, 951.879, 1438.295, 2146.405, 3060.923];  

lambda_asymp_selected = lambda_asymp(x);

% Plot
figure('Position',[100 100 950 620],'Color','white');
hold on; grid on; box on;

semilogy(k, lambda_asymp, 'b-', 'LineWidth', 2);
semilogy(x, lambda_asymp_selected, 'bo','MarkerSize', 12, 'LineWidth', 2,'MarkerFaceColor', 'w');
semilogy(x, y, 'ro', 'MarkerSize', 9, 'LineWidth', 1.5,'MarkerFaceColor', 'r');


ylim([0 4000]); 

xlabel('k');
ylabel('\lambda_k');

h1 = plot(NaN, NaN, 'bo-', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'w');
h2 = plot(NaN, NaN, 'ro', 'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'r');


legend([h1 h2], ...
       {'Asymptotic approximation', 'Matslise'}, ...
       'Location', 'north','Box', 'on','FontSize',17);

set(gca, 'FontSize', 16);
set(gca, 'YMinorGrid', 'on');

print('-depsc', 'fix_ce_1_50.eps')

