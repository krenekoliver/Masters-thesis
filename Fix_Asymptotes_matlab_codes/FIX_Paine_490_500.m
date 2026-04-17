% ================================================
% Asymptotic approximation λ_k for p(x) = e^x
% Dirichlet: u(0)=u(π)=0
% ================================================

clear; clc; close all;

n = 490:500;

P1 = exp(pi) - 1;
P2 = 0.5*(exp(2*pi) - 1);

% Asymptotic approximation
sqrt_lambda = n + (2*pi*n).^(-1)*P1 + (4*pi*n.^3).^(-1)*(P2/2 - P1^2/pi + exp(pi)/2 - 1/2);
lambda_asymp = sqrt_lambda.^2;

% Numeric values from Matslise
x = [490 491 492 493 494 495 496 497 498 499 500];
y = [240107.047 241088.047 242071.047 243056.047 244043.047 245032.047 246023.047 247016.047 248011.047 249008.047 250007.047];

[~, idx] = ismember(x, n);
lambda_asymp_selected = lambda_asymp(idx);

% Plot
figure('Position',[100 100 950 620],'Color','white');
hold on; grid on; box on;

semilogy(n, lambda_asymp, 'b-', 'LineWidth', 2);
semilogy(x, lambda_asymp_selected, 'bo','MarkerSize', 12, 'LineWidth', 2,'MarkerFaceColor', 'w');
semilogy(x, y, 'ro', 'MarkerSize', 9, 'LineWidth', 1.5,'MarkerFaceColor', 'r');

xlabel('k');
ylabel('\lambda_k');

h1 = plot(NaN, NaN, 'bo-', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'w');
h2 = plot(NaN, NaN, 'ro', 'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'r');


legend([h1 h2], ...
       {'Asymptotic approximation', 'Matslise'}, ...
       'Location', 'northwest','Box', 'on','FontSize',17);

set(gca, 'FontSize', 16);
set(gca, 'YMinorGrid', 'on');

print('-depsc', 'fix_paine_490_500.eps')