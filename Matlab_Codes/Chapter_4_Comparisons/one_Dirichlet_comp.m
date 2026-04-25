clear all; close all; clc;

function str = formatThousands(num)
    str = num2str(round(num), '%d');
    str = regexprep(str, '(\d)(?=(\d{3})+$)', '$1 ');
end

k_wave = 8;
m2 = 0.000000001;               
plot_time = 500000.007;
N = 256;
I = 2*pi;
k1 = 5;
k2 = 50;
m1 = 1;

h = I / (N + 1);
x_inner = h * (1:N)';
x_full  = [0; x_inner; I];

c_limit = I * sqrt(k1 / m1);

%% initial conditions
u0 = sin(k_wave * pi * x_inner / I);  
v0 = -c_limit * (k_wave * pi / I) * cos(k_wave * pi * x_inner / I);
a0 = zeros(N,1);

y0_exact = [u0; a0];
ydot0    = [v0; zeros(N,1)];

tspan = linspace(0, 500000, 1000);

%% Wave equation (DST)
k_vec = pi * (1:N)' / I;
omega_classic = c_limit * k_vec;

u0_hat = dst(u0);
v0_hat = dst(v0);

u_classic_inner = @(t) idst( u0_hat .* cos(omega_classic*t) + (v0_hat ./ max(omega_classic,1e-12)) .* sin(omega_classic*t) );

u_classic_fun = @(t) [0; u_classic_inner(t); 0];

%% Rate-type metamaterial (DST)
D1 = diag(1./m1);
D2 = diag(1./m2);

Lambda = diag(-k_vec.^2);
F = dst(eye(N));
A = idst(Lambda * F);

TopLeft     = (k1 * I^2) * (D1 * A);
TopRight    = D1;
BottomLeft  = -k2 * (k1 * I^2) * (D1 * A);
BottomRight = -k2 * (D1 + D2);

M = [TopLeft, eye(N)*TopRight; BottomLeft, eye(N)*BottomRight];

[V, D] = eig(full(M));
omega_meta = sqrt(-diag(D) + 1e-12i);

coeff_cos = V \ y0_exact;
coeff_sin = (V \ ydot0) ./ omega_meta;

u_meta_inner = @(t) real(V(1:N,:) * (coeff_cos .* cos(omega_meta*t) + coeff_sin .* sin(omega_meta*t)));

u_meta_fun = @(t) [0; u_meta_inner(t); 0];

%% PostScript
figure('Position',[80 80 1480 760],'Color','white');
hold on;
plot(x_full, u_classic_fun(plot_time), 'b-', 'LineWidth',2.5, 'DisplayName','Classical wave');
plot(x_full, u_meta_fun(plot_time), 'r-', 'LineWidth',2.5, 'DisplayName','Metamaterial');

xlabel('x'); ylabel('u(x,t)');
title(sprintf('m_2 = 10^{-9} | t = %s', formatThousands(plot_time)));
ax = gca;
ax.XTick = 0 : pi/2 : I;
ax.XTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};
ax.YTick = -1.5 : 0.5 : 1.5;
ax.YTickLabel = {'-3/2','-1','-1/2','0','1/2','1','3/2'};
ax.FontSize = 25;
ylim([-1.5 1.5]);
xlim([0 I]);
grid on;
% legend('Location','northwest');

print('-depsc2', '-r300', sprintf('comparison_one_dirichlet_10_9.eps'));

%% Video
video_filename = fullfile('videos', sprintf('comparison_dirichlet_initial.avi'));
vid = VideoWriter(video_filename, 'Uncompressed AVI');
vid.FrameRate = 30;
open(vid);

figure('Position',[80 80 1480 760],'Color','white');
for i = 1:length(tspan)
    t = tspan(i);
    clf; hold on;
    plot(x_full, u_classic_fun(t), 'b-', 'LineWidth',2.5);
    plot(x_full, u_meta_fun(t), 'r-', 'LineWidth',2.5);
    
    xlabel('x'); ylabel('u(x,t)');
    title(sprintf('Dirichlet BC | t = %s', formatThousands(t)));
    ax = gca;
    ax.XTick = 0 : pi/2 : I;
    ax.XTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};
    ax.YTick = -1.5 : 0.5 : 1.5;
    ax.YTickLabel = {'-3/2','-1','-1/2','0','1/2','1','3/2'};
    ax.FontSize = 25;
    ylim([-1.5 1.5]);
    xlim([0 I]);
    grid on;
    %legend('Classical wave', 'Metamaterial', 'Location','northwest');
    
    frame = getframe(gcf);
    writeVideo(vid, frame);
end
close(vid);