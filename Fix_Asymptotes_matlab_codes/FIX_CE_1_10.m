% ================================================
% Asymptotic approximation λ_k for Coffey-Evans equation (β = 30)
% Dirichlet: [-π/2, π/2] -> [0, π]
% ================================================

clear; clc; close all;

k = 1:10;

n = k + 1;                    

P1 = 450 * pi;
P2 = 305550 * pi;

% Asymptotic approximation
sqrt_lambda = n + (P1 ./ (2*pi*n)) + (1./(4*pi*n.^3)) .* (P2/2 - P1.^2 / pi);
lambda_asymp = sqrt_lambda .^ 2;

% Numeric values from Matslise
x = [1 2 3 4 5 6 7 8 9 10];
y = [117.946, 231.665, 231.665, 231.665, 340.888, 445.283, 445.283, 445.283, 544.418, 648.418];  


% Plot
figure('Position',[100 100 950 620],'Color','white');
hold on; grid on; box on;

semilogy(k, lambda_asymp, 'bo-','MarkerSize', 12, 'LineWidth', 2,'MarkerFaceColor', 'w', 'DisplayName', 'Asymptotic approximation');
semilogy(x, y, 'ro', 'MarkerSize', 9, 'LineWidth', 1.5,'MarkerFaceColor', 'r', 'DisplayName', 'Matslise');
ylim([0 900]);
xlabel('k');
ylabel('\lambda_k');

legend('Location', 'northeast','FontSize',17);

set(gca, 'FontSize', 16);
set(gca, 'YMinorGrid', 'on');

print('-depsc', 'fix_ce_1_10.eps')

