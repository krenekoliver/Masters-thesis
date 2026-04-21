clear all; close all; clc;

k_wave = 16;                    
plot_time = 500000;         

function str = formatThousands(num)
    str = num2str(round(num), '%d');
    str = regexprep(str, '(\d)(?=(\d{3})+$)', '$1 ');
end

N = 30;
I = 2*pi;
L = I / (N + 1);
x_uniform = L * (1:N)';             
x_full    = [0; x_uniform; I];

k1 = 5;
k2 = 50;
m1 = ones(1,N);
m2 = 0.1*ones(1,N);

%% initial conditions
kappa = k_wave * pi / I;   

u0 = sin(kappa * x_uniform) .* (x_uniform <= I/2);
v0 = -kappa * cos(kappa * x_uniform) .* (x_uniform <= I/2);   
a0 = zeros(N,1);

y0_exact = [u0; a0];
ydot0 = [v0; zeros(N,1)];
tspan = linspace(0, 50, 500);

D1 = diag(1./m1); 
D2m = diag(1./m2);

%% FD
A_fd = spdiags(ones(N,1)*[1 -2 1], -1:1, N, N) / L^2;
M_fd = [k1*D1*A_fd, D1; -k2*k1*D1*A_fd, -k2*(D2m+D1)];
[V_fd, D_fd] = eig(full(M_fd));
omega_fd = sqrt(-diag(D_fd) + 1e-12i);
coeff_cos_fd = V_fd \ y0_exact;
coeff_sin_fd = (V_fd \ ydot0) ./ omega_fd;
u_fd_exact = @(t) real(V_fd(1:N,:) * (coeff_cos_fd .* cos(omega_fd*t) + coeff_sin_fd .* sin(omega_fd*t)));

u_fd_plot = @(t) [0; u_fd_exact(t); 0];

%% Numerov
M_num = spdiags(ones(N,1)*[1 10 1], -1:1, N, N);
K_num = spdiags(ones(N,1)*[1 -2 1], -1:1, N, N);
A_num = (12/L^2) * (M_num \ K_num);
M_num_full = [k1*D1*A_num, D1; -k2*k1*D1*A_num, -k2*(D2m+D1)];
[V_num, D_num] = eig(full(M_num_full));
omega_num = sqrt(-diag(D_num) + 1e-12i);
coeff_cos_num = V_num \ y0_exact;
coeff_sin_num = (V_num \ ydot0) ./ omega_num;
u_num_exact = @(t) real(V_num(1:N,:) * (coeff_cos_num .* cos(omega_num*t) + coeff_sin_num .* sin(omega_num*t)));

u_num_plot = @(t) [0; u_num_exact(t); 0];

%% DST
k_vec = pi * (1:N)' / I;
Lambda = diag(-k_vec.^2);
S = dst(eye(N));
A_dst = real(idst(Lambda * S));
M_dst = [k1*D1*A_dst, D1; -k2*k1*D1*A_dst, -k2*(D2m+D1)];
[V_dst, D_dst] = eig(full(M_dst));
omega_dst = sqrt(-diag(D_dst) + 1e-12i);
coeff_cos_dst = V_dst \ y0_exact;
coeff_sin_dst = (V_dst \ ydot0) ./ omega_dst;
u_dst_exact = @(t) real(V_dst(1:N,:) * (coeff_cos_dst .* cos(omega_dst*t) + coeff_sin_dst .* sin(omega_dst*t)));

u_dst_plot = @(t) [0; u_dst_exact(t); 0];

%% Chebyshev methods
dom = [0, I];

N_cheb_total = N + 2;

D2_cheb_full = diffmat(N_cheb_total, 2, dom, 'dirichlet', 'dirichlet');

A_cheb = D2_cheb_full(2:end-1, 2:end-1);

x_cheb_full = chebpts(N_cheb_total, dom);
x_cheb = x_cheb_full(2:end-1);

kappa = k_wave * pi / I;
u0_cheb = sin(kappa * x_cheb) .* (x_cheb <= I/2);
v0_cheb = -kappa * cos(kappa * x_cheb) .* (x_cheb <= I/2);

y0_cheb   = [u0_cheb;   zeros(N,1)];
ydot_cheb = [v0_cheb;   zeros(N,1)];

%% CDM
M_cdm_cheb = [k1*D1*A_cheb, D1; -k2*k1*D1*A_cheb, -k2*(D2m+D1)];
[V_cdm_cheb, D_cdm_cheb] = eig(full(M_cdm_cheb));
omega_cdm_cheb = sqrt(-diag(D_cdm_cheb) + 1e-12i);
coeff_cos_cdm_cheb = V_cdm_cheb \ y0_cheb;
coeff_sin_cdm_cheb = (V_cdm_cheb \ ydot_cheb) ./ omega_cdm_cheb;
u_cdm_cheb_exact = @(t) real(V_cdm_cheb(1:N,:) * (coeff_cos_cdm_cheb .* cos(omega_cdm_cheb*t) + coeff_sin_cdm_cheb .* sin(omega_cdm_cheb*t)));

%% MBCDM
A_mbcdm = A_cheb;
M_mbcdm = [k1*D1*A_mbcdm, D1; -k2*k1*D1*A_mbcdm, -k2*(D2m+D1)];
[V_mbcdm, D_mbcdm] = eig(full(M_mbcdm));
omega_mbcdm = sqrt(-diag(D_mbcdm) + 1e-12i);
coeff_cos_mbcdm = V_mbcdm \ y0_cheb;
coeff_sin_mbcdm = (V_mbcdm \ ydot_cheb) ./ omega_mbcdm;
u_mbcdm_exact = @(t) real(V_mbcdm(1:N,:) * (coeff_cos_mbcdm .* cos(omega_mbcdm*t) + coeff_sin_mbcdm .* sin(omega_mbcdm*t)));

%% LGCC
A_lgcc = A_cheb;
M_lgcc = [k1*D1*A_lgcc, D1; -k2*k1*D1*A_lgcc, -k2*(D2m+D1)];
[V_lgcc, D_lgcc] = eig(full(M_lgcc));
omega_lgcc = sqrt(-diag(D_lgcc) + 1e-12i);
coeff_cos_lgcc = V_lgcc \ y0_cheb;
coeff_sin_lgcc = (V_lgcc \ ydot_cheb) ./ omega_lgcc;
u_lgcc_exact = @(t) real(V_lgcc(1:N,:) * (coeff_cos_lgcc .* cos(omega_lgcc*t) + coeff_sin_lgcc .* sin(omega_lgcc*t)));

%% Exact solution
N_modes = 512;
x_fine = linspace(0, I, 4096)';
u0_fine = sin(kappa * x_fine) .* (x_fine <= I/2);
v0_fine = -kappa * cos(kappa * x_fine) .* (x_fine <= I/2);

mu1 = m1(1); mu2 = m2(1);

modal_omega1 = zeros(N_modes,1); modal_omega2 = zeros(N_modes,1);
modal_cos1 = zeros(N_modes,1); modal_sin1 = zeros(N_modes,1);
modal_cos2 = zeros(N_modes,1); modal_sin2 = zeros(N_modes,1);

for n = 1:N_modes
    kappa_n = n * pi / I;
    A_op = -kappa_n^2;
    Mat = [k1/mu1*A_op, 1/mu1; -k2*k1/mu1*A_op, -k2*(1/mu1 + 1/mu2)];
    [V_mat, D_mat] = eig(Mat);
    omegas = real(sqrt(-real(diag(D_mat))));
    [omegas, idx] = sort(omegas);
    V_mat = V_mat(:,idx);
    
    modal_omega1(n) = omegas(1);
    modal_omega2(n) = omegas(2);
    
    sin_n = sin(n * pi * x_fine / I);
    u_n0 = (2/I) * trapz(x_fine, u0_fine .* sin_n);
    v_n0 = (2/I) * trapz(x_fine, v0_fine .* sin_n);
    
    y0_n   = [u_n0; 0];
    dot0_n = [v_n0; 0];
    C = V_mat \ y0_n;
    B = (V_mat \ dot0_n) ./ omegas;
    
    modal_cos1(n) = V_mat(1,1) * C(1);
    modal_sin1(n) = V_mat(1,1) * B(1);
    modal_cos2(n) = V_mat(1,2) * C(2);
    modal_sin2(n) = V_mat(1,2) * B(2);
end

function u = modal_exact(t, x, c1c, c1s, c2c, c2s, w1, w2, L)
    u = zeros(size(x));
    for n = 1:length(c1c)
        u = u + (c1c(n)*cos(w1(n)*t) + c1s(n)*sin(w1(n)*t) + c2c(n)*cos(w2(n)*t) + c2s(n)*sin(w2(n)*t)) .* sin(n * pi * x / L);
    end
end

u_exact_fun = @(t,x) modal_exact(t, x, modal_cos1, modal_sin1, modal_cos2, modal_sin2, modal_omega1, modal_omega2, I);

u_exact_plot = @(t) u_exact_fun(t, x_full);

%% ====================== PostScript Image ======================
figure('Position',[80 80 1480 760],'Color','white');
hold on;

%plot(x_full, u_fd_plot(plot_time),      'b-',   'LineWidth',2.5, 'DisplayName','FD');
% plot(x_full, u_num_plot(plot_time),     'r-',   'LineWidth',2.5, 'DisplayName','Numerov');
plot(x_full, u_dst_plot(plot_time),     'm-',   'LineWidth',2.5, 'DisplayName','DST');
% plot(x_cheb, u_cdm_cheb_exact(plot_time),'Color', [0.6 0.4 0.2],   'LineWidth',2.5, 'DisplayName','CDM');
plot(x_cheb, u_mbcdm_exact(plot_time),   'Color', [0.93 0.45 0.15],   'LineWidth',2.5, 'DisplayName','MBCDM');
%plot(x_cheb, u_lgcc_exact(plot_time),    'Color', [0 0.6 0.3], 'LineWidth',2.5, 'DisplayName','LGCC');
 plot(x_full, u_exact_plot(plot_time),    'k--', 'LineWidth',2.5, 'DisplayName','Exact');

xlabel('x');
ylabel('u(x,t)');
title(sprintf('Num | t = %s', formatThousands(plot_time)));
ax = gca;
ax.XTick = 0 : pi/2 : 2*pi;
ax.XTickLabel = {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'};
ax.YTick = -1 : 0.5 : 1;
ax.YTickLabel = {'-1', '-1/2', '0', '1/2', '1'};
ax.FontSize = 25;
ylim([-1 1]);
xlim([0 2*pi]);
grid on;
% legend('Location','northwest');

print('-depsc2', '-r300', 'one_dirichlet_meta_Num.eps');

%% ====================== VIDEO ======================
video_filename = fullfile('videos', 'one_dirichlet_meta.avi');
vid = VideoWriter(video_filename, 'Uncompressed AVI');
vid.FrameRate = 10;
open(vid);

figure('Position',[80 80 1480 760],'Color','white');
for i = 1:length(tspan)
    t = tspan(i);
    clf; hold on;
    
     plot(x_full, u_fd_plot(t),      'b-',   'LineWidth',2.5);
    % plot(x_full, u_num_plot(t),     'r-',   'LineWidth',2.5);
    % plot(x_full, u_dst_plot(t),     'm-',   'LineWidth',2.5);
    % plot(x_cheb, u_cdm_cheb_exact(t),'Color', [0.6 0.4 0.2],   'LineWidth',2.5);
    %  plot(x_cheb, u_mbcdm_exact(t),   'Color', [0.93 0.45 0.15],   'LineWidth',2.5);
    % plot(x_cheb, u_lgcc_exact(t),    'Color', [0 0.6 0.3], 'LineWidth',2.5);
     plot(x_full, u_exact_plot(t),    'k--', 'LineWidth',2.5);

    xlabel('x');
    ylabel('u(x,t)');
    title(sprintf('FD | t = %s', formatThousands(t)));
    ax = gca;
    ax.XTick = 0 : pi/2 : 2*pi;
    ax.XTickLabel = {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'};
    ax.YTick = -2 : 1 : 2;
    ax.FontSize = 25;
    ylim([-2 2]);
    xlim([0 2*pi]);
    grid on;
    % legend('FD','Exact','Location','northwest');
    
    frame = getframe(gcf);
    writeVideo(vid, frame);
end
close(vid);
