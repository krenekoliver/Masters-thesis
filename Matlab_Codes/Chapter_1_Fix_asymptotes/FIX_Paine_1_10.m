% ================================================
% Asymptotic approximation λ_k for p(x) = e^x
% Dirichlet: u(0)=u(π)=0
% ================================================

clear; clc; close all;

n = 1:10;

P1 = exp(pi) - 1;
P2 = 0.5*(exp(2*pi) - 1);

% Asymptotic approximation
sqrt_lambda = n + (2*pi*n).^(-1)*P1 + (4*pi*n.^3).^(-1)*(P2/2 - P1^2/pi + exp(pi)/2 - 1/2);
lambda_asymp = sqrt_lambda.^2;

% Numerical values from Matslise
x = [1 2 3 4 5 6 7 8 9 10];
y = [4.89666939017 10.04518988883256 16.0192672359282 23.2662709260240 32.2637070872960 43.2200196633159 56.181594044615 71.1529975070128 88.1321191545589 107.1166761105273];

% Plot
figure('Position',[100 100 950 620],'Color','white');
hold on; grid on; box on;

semilogy(n, lambda_asymp, 'bo-','MarkerSize', 12, 'LineWidth', 2,'MarkerFaceColor', 'w', 'DisplayName', 'Asymptotic approximation');
semilogy(x, y, 'ro', 'MarkerSize', 9, 'LineWidth', 1.5,'MarkerFaceColor', 'r', 'DisplayName', 'Matslise');

xlabel('k');
ylabel('\lambda_k');

legend('Location', 'northwest','FontSize',17);

set(gca, 'FontSize', 16);
set(gca, 'YMinorGrid', 'on');

print('-depsc', 'fix_paine_1_10.eps')