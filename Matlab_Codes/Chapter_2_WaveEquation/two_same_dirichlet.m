clear all; close all; clc;

function str = formatThousands(num)
    str = num2str(round(num), '%d');           
    str = regexprep(str, '(\d)(?=(\d{3})+$)', '$1 ');  
end

I       = pi;
c       = 1.0;
k1      = 16;                     
k2      = 24;                     
T_final = 50;                
N       = 50;                   

t_postscript = 12.57;
postscript_filename = sprintf('pic_two_dirichlet_same_FD.eps');

x_sin   = (1:N)' * I / (N + 1);

%% initial conditions
period = pi/2;

u0_sin  = zeros(N,1);
ut0_sin = zeros(N,1);

mask = x_sin <= period;

u0_sin(mask)  = sin(k1 * x_sin(mask)) + sin(k2 * x_sin(mask));
ut0_sin(mask) = -c*k1 * cos(k1 * x_sin(mask)) - c*k2 * cos(k2 * x_sin(mask));

%% DST
b0 = dst(u0_sin);
g0 = dst(ut0_sin);
m_vec = (1:N)';

%% FD
dx = I / (N + 1);
T = spdiags(ones(N,1)*[1 -2 1], -1:1, N, N);
A_fd = (c^2 / dx^2) * T;
[V_fd, D_fd] = eig(full(A_fd));
omega_fd = sqrt(-diag(D_fd));
coeff_u0_fd = V_fd \ u0_sin;
coeff_ut0_fd = V_fd \ ut0_sin;

%% Numerov
K = spdiags(ones(N,1)*[1 -2 1], -1:1, N, N);
M = spdiags(ones(N,1)*[1 10 1], -1:1, N, N);
const = 12 / dx^2;
A_num = M \ (const * K);
[V_num, D_num] = eig(full(A_num));
omega_num = c * sqrt(-diag(D_num));
coeff_u0_num  = V_num \ u0_sin;
coeff_ut0_num = V_num \ ut0_sin;

%% Chebyshev methods (CDM, MBCDM, LGCC)
dom = [0 I];
D2_cheb = diffmat(N+2, 2, dom);
D2_int  = D2_cheb(2:end-1, 2:end-1);

[x_cheb, ~] = chebpts(N+2, dom);
x_cheb_int = x_cheb(2:end-1);

u0_cdm  = zeros(size(x_cheb_int));
ut0_cdm = zeros(size(x_cheb_int));
mask_cheb = x_cheb_int <= period;

u0_cdm(mask_cheb)  = sin(k1 * x_cheb_int(mask_cheb)) + sin(k2 * x_cheb_int(mask_cheb));
ut0_cdm(mask_cheb) = -c*k1 * cos(k1 * x_cheb_int(mask_cheb)) - c*k2 * cos(k2 * x_cheb_int(mask_cheb));

% CDM
[V_cdm, D_cdm] = eig(full(D2_int));
omega_cdm = c * sqrt(-diag(D_cdm));
coeff_u0_cdm  = V_cdm \ u0_cdm;
coeff_ut0_cdm = V_cdm \ ut0_cdm;

% MBCDM
[V_mbcd, D_mbcd] = eig(full(D2_int));
omega_mbcd = c * sqrt(-diag(D_mbcd));
coeff_u0_mbcd  = V_mbcd \ u0_cdm;
coeff_ut0_mbcd = V_mbcd \ ut0_cdm;

% LGCC
[V_lgcc, D_lgcc] = eig(full(D2_int));
omega_lgcc = c * sqrt(-diag(D_lgcc));
coeff_u0_lgcc  = V_lgcc \ u0_cdm;
coeff_ut0_lgcc = V_lgcc \ ut0_cdm;

%% Exact solution
N_exact = 1024;
x_exact = (1:N_exact)' * I / (N_exact + 1);
u0_exact = zeros(N_exact,1);
ut0_exact = zeros(N_exact,1);

mask_exact = x_exact <= period;

u0_exact(mask_exact)  = sin(k1 * x_exact(mask_exact)) + sin(k2 * x_exact(mask_exact));
ut0_exact(mask_exact) = -c*k1 * cos(k1 * x_exact(mask_exact)) - c*k2 * cos(k2 * x_exact(mask_exact));

b0_exact = dst(u0_exact);
g0_exact = dst(ut0_exact);
m_exact_vec = (1:N_exact)';

%% Video
video_filename = fullfile('videos', 'two_same_dirichlet.avi');
vid = VideoWriter(video_filename, 'Uncompressed AVI');
vid.FrameRate = 40;
open(vid);

dt_plot = 10;
n_frames = round(T_final / dt_plot);

figure('Position', [80 80 1480 760],'Color','white');

for n = 1:n_frames
    t = (n-1) * dt_plot;
    
    u_dst = idst( b0 .* cos(c * m_vec * t) + (g0 ./ (c * m_vec)) .* sin(c * m_vec * t) );
    u_fd = V_fd * (coeff_u0_fd .* cos(omega_fd * t) + (coeff_ut0_fd ./ omega_fd) .* sin(omega_fd * t));
    u_num = V_num * (coeff_u0_num .* cos(omega_num * t) + (coeff_ut0_num ./ omega_num) .* sin(omega_num * t));
    
    u_cdm  = V_cdm  * (coeff_u0_cdm  .* cos(omega_cdm  * t) + (coeff_ut0_cdm  ./ omega_cdm)  .* sin(omega_cdm  * t));
    u_mbcd = V_mbcd * (coeff_u0_mbcd .* cos(omega_mbcd * t) + (coeff_ut0_mbcd ./ omega_mbcd) .* sin(omega_mbcd * t));
    u_lgcc = V_lgcc * (coeff_u0_lgcc .* cos(omega_lgcc * t) + (coeff_ut0_lgcc ./ omega_lgcc) .* sin(omega_lgcc * t));
    
    u_exact = idst( b0_exact .* cos(c * m_exact_vec * t) + (g0_exact ./ (c * m_exact_vec)) .* sin(c * m_exact_vec * t) );
    
    clf; hold on;
    
    plot(x_exact, u_exact,'k--', 'LineWidth', 2.5, 'DisplayName', 'Exact');
    
    x_dst_plot = [0; x_sin; I];
    u_dst_plot = [0; u_dst; 0];
    % plot(x_dst_plot, u_dst_plot, 'm-', 'LineWidth', 2.5, 'DisplayName', 'DST');
    
    x_fd_plot = [0; x_sin; I];
    u_fd_plot = [0; u_fd; 0];
    % plot(x_fd_plot, u_fd_plot, 'b-', 'LineWidth', 2.5, 'DisplayName', 'FD');
    
    x_num_plot = [0; x_sin; I];
    u_num_plot = [0; u_num; 0];
     plot(x_num_plot, u_num_plot, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Numerov');
    
    % plot(x_cheb_int, u_cdm,  'Color', [0.6 0.4 0.2], 'LineWidth', 2.5, 'DisplayName', 'CDM');
    % plot(x_cheb_int, u_mbcd, 'Color', [0.93 0.45 0.15], 'LineWidth', 2.5, 'DisplayName', 'MBCD');
    % plot(x_cheb_int, u_lgcc, 'Color', [0 0.6 0.3], 'LineWidth', 2.5, 'DisplayName', 'LGCC');
    
    xlabel('x');
    ylabel('u(x,t)');
    title(sprintf('DFT | t = %s', formatThousands(t)));
    
    ax = gca;
    ax.XTick = 0 : pi/4 : I;
    ax.XTickLabel = {'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'};
    ax.YTick = -2 : 1 : 2;
    ax.YTickLabel = {'-2', '-1','0', '1', '2'};
    
    ax.FontSize = 25;
    ax.LineWidth = 1.2;

    ylim([-2 2]);
    xlim([0 I]);
    grid on;
    hold off;
    
    frame = getframe(gcf);
    writeVideo(vid, frame);
end

close(vid);

% === PostScript (EPS) ===
t = t_postscript;

u_dst = idst( b0 .* cos(c * m_vec * t) + (g0 ./ (c * m_vec)) .* sin(c * m_vec * t) );
u_fd = V_fd * (coeff_u0_fd .* cos(omega_fd * t) + (coeff_ut0_fd ./ omega_fd) .* sin(omega_fd * t));
u_num = V_num * (coeff_u0_num .* cos(omega_num * t) + (coeff_ut0_num ./ omega_num) .* sin(omega_num * t));

u_cdm  = V_cdm  * (coeff_u0_cdm  .* cos(omega_cdm  * t) + (coeff_ut0_cdm  ./ omega_cdm)  .* sin(omega_cdm  * t));
u_mbcd = V_mbcd * (coeff_u0_mbcd .* cos(omega_mbcd * t) + (coeff_ut0_mbcd ./ omega_mbcd) .* sin(omega_mbcd * t));
u_lgcc = V_lgcc * (coeff_u0_lgcc .* cos(omega_lgcc * t) + (coeff_ut0_lgcc ./ omega_lgcc) .* sin(omega_lgcc * t));

u_exact = idst( b0_exact .* cos(c * m_exact_vec * t) + (g0_exact ./ (c * m_exact_vec)) .* sin(c * m_exact_vec * t) );

figure(1); clf; hold on;

plot(x_exact, u_exact,'k--', 'LineWidth', 2.5, 'DisplayName', 'Exact');

x_dst_plot = [0; x_sin; I];
u_dst_plot = [0; u_dst; 0];
% plot(x_dst_plot, u_dst_plot, 'm-', 'LineWidth', 2.5, 'DisplayName', 'DST');

x_fd_plot = [0; x_sin; I];
u_fd_plot = [0; u_fd; 0];
 plot(x_fd_plot, u_fd_plot, 'b-', 'LineWidth', 2.5, 'DisplayName', 'FD');

x_num_plot = [0; x_sin; I];
u_num_plot = [0; u_num; 0];
% plot(x_num_plot, u_num_plot, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Numerov');

% plot(x_cheb_int, u_cdm,  'Color', [0.6 0.4 0.2], 'LineWidth', 2.5, 'DisplayName', 'CDM');
% plot(x_cheb_int, u_mbcd, 'Color', [0.93 0.45 0.15], 'LineWidth', 2.5, 'DisplayName', 'MBCD');
% plot(x_cheb_int, u_lgcc, 'Color', [0 0.6 0.3], 'LineWidth', 2.5, 'DisplayName', 'LGCC');

xlabel('x');
ylabel('u(x,t)');
title(sprintf('FD | t = %s', formatThousands(t)));

ax = gca;
ax.XTick = 0 : pi/4 : I;
ax.XTickLabel = {'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'};
ax.YTick = -2 : 1 : 2;
ax.YTickLabel = {'-2', '-1','0', '1', '2'};

ax.FontSize = 25;
ax.LineWidth = 1.2;

ylim([-2 2]);
xlim([0 I]);
grid on;
hold off;

% === PostScript (EPS) ===
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Renderer', 'Painters');
print(gcf, postscript_filename, '-depsc');