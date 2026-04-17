% ================================================
% Asymptotic approximation λ_k for Coffey-Evans equation (β = 30)
% Dirichlet: [-π/2, π/2] -> [0, π]
% ================================================

clear; clc; close all;

k = 490:500;

n = k + 1;                    

P1 = 450 * pi;
P2 = 305550 * pi;

% Asymptotic approximation
sqrt_lambda = n + (P1 ./ (2*pi*n)) + (1./(4*pi*n.^3)) .* (P2/2 - P1.^2 / pi);
lambda_asymp = sqrt_lambda .^ 2;

% Numeric values from Matslise
x = [490 491 492 493 494 495 496 497 498 499 500];
y = [241531.107, 242514.106, 243499.106, 244486.106, 245475.105, 246466.105, 247459.104, 248454.104, 249451.103, 250450.103 ,251451.103 ];

[~, idx] = ismember(x, k);
lambda_asymp_selected = lambda_asymp(idx);

% Plot
figure('Position',[100 100 950 620],'Color','white');
hold on; grid on; box on;

semilogy(k, lambda_asymp, 'b-', 'LineWidth', 2);
semilogy(x, lambda_asymp_selected, 'bo','MarkerSize', 12, 'LineWidth', 2,'MarkerFaceColor', 'w');
semilogy(x, y, 'ro', 'MarkerSize', 9, 'LineWidth', 1.5,'MarkerFaceColor', 'r');

ylim([200000 300000]); 

xlabel('k');
ylabel('\lambda_k');

h1 = plot(NaN, NaN, 'bo-', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'w');
h2 = plot(NaN, NaN, 'ro', 'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'r');


legend([h1 h2], ...
       {'Asymptotic approximation', 'Matslise'}, ...
       'Location', 'north','Box', 'on','FontSize',17);

set(gca, 'FontSize', 16);
set(gca, 'YMinorGrid', 'on');

print('-depsc','fix_ce_490_500.eps')
