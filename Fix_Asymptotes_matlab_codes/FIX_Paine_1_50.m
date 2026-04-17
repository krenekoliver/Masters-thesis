% ================================================
% Asymptotic approximation λ_k for p(x) = e^x
% Dirichlet: u(0)=u(π)=0
% ================================================

clear; clc; close all;

n = 1:50;

P1 = exp(pi) - 1;
P2 = 0.5*(exp(2*pi) - 1);

% Asymptotic approximation
sqrt_lambda = n + (2*pi*n).^(-1)*P1 + (4*pi*n.^3).^(-1)*(P2/2 - P1^2/pi + exp(pi)/2 - 1/2);
lambda_asymp = sqrt_lambda.^2;

% Numerical values from Matslise
x = [1 10 20 30 40 50];
y = [4.89666939017 107.1166761105273 407.06523522670701 907.0554605783833 1607.0520261907873 2507.0504344089];

lambda_asymp_selected = lambda_asymp(x);

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

print('-depsc', 'fix_paine_1_50.eps')