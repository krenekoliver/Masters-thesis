clear all; close all; clc;

plot_time = 499999.831;

function str = formatThousands(num)
    str = num2str(round(num), '%d');
    str = regexprep(str, '(\d)(?=(\d{3})+$)', '$1 ');
end

N = 128;
I = 2*pi;
L = I/N;

k1 = 5;
k2 = 50;
m1 = ones(1,N);
m2 = 0.1*ones(1,N);

x = I/N * (-N/2 : N/2-1)';

k_wave1 = 4;
k_wave2 = 8;
amp1    = 1.0;
amp2    = 1.0;

u0 = amp1 * sin(k_wave1 * x) + amp2 * sin(k_wave2 * x);
v0 = -amp1*k_wave1*cos(k_wave1*x) - amp2*k_wave2*cos(k_wave2*x);
a0 = zeros(N,1);      

%% initial conditions
y0_exact = [u0; a0];                
ydot0     = [v0; zeros(N,1)];      

tspan = linspace(0, 30, 1000);

%% FD
A_fd = toeplitz([-2 1 zeros(1,N-3) 1]) / L^2;
D1 = diag(1./m1); D2 = diag(1./m2);
M_fd = [k1*D1*A_fd, D1; -k2*k1*D1*A_fd, -k2*(D2+D1)];
[V_fd, D_fd] = eig(full(M_fd));
lambda_fd = diag(D_fd);
omega_fd = sqrt(-lambda_fd + 1e-12i);
coeff_cos_fd = V_fd \ y0_exact;
coeff_sin_fd = (V_fd \ ydot0) ./ omega_fd;
u_fd_exact = @(t) real(V_fd(1:N,:) * (coeff_cos_fd .* cos(omega_fd*t) + coeff_sin_fd .* sin(omega_fd*t)));

%% Numerov
M_num = spdiags(ones(N,1)*[1 10 1], -1:1, N, N); 
M_num(1,N)=1; M_num(N,1)=1;
K_num = spdiags(ones(N,1)*[1 -2 1], -1:1, N, N); 
K_num(1,N)=1; K_num(N,1)=1;
A_num = (12/L^2) * (M_num \ K_num);
D1 = diag(1./m1); D2 = diag(1./m2);
M_num_full = [k1*D1*A_num, D1; -k2*k1*D1*A_num, -k2*(D2+D1)];
[V_num, D_num] = eig(full(M_num_full));
lambda_num = diag(D_num);
omega_num = sqrt(-lambda_num + 1e-12i);
coeff_cos_num = V_num \ y0_exact;
coeff_sin_num = (V_num \ ydot0) ./ omega_num;
u_num_exact = @(t) real(V_num(1:N,:) * (coeff_cos_num .* cos(omega_num*t) + coeff_sin_num .* sin(omega_num*t)));

%% DFT
k_vec = (2*pi/I) * [0:(N/2-1), -N/2:-1]';
Lambda = diag(-k_vec.^2);
F = fft(eye(N));
A_dft = ifft(Lambda * F);
A_dft = real(A_dft);
D1 = diag(1./m1); D2 = diag(1./m2);
M_dft = [k1*D1*A_dft, D1; -k2*k1*D1*A_dft, -k2*(D2+D1)];
[V_dft, D_dft] = eig(full(M_dft));
lambda_dft = diag(D_dft);
omega_dft = sqrt(-lambda_dft + 1e-12i);
coeff_cos_dft = V_dft \ y0_exact;
coeff_sin_dft = (V_dft \ ydot0) ./ omega_dft;
u_dft_exact = @(t) real(V_dft(1:N,:) * (coeff_cos_dft .* cos(omega_dft*t) + coeff_sin_dft .* sin(omega_dft*t)));

%% Exact
mu1 = m1(1);
mu2 = m2(1);

function u_mode = exact_mode(t, x, k, amp, mu1, mu2, k1, k2)
    kappa = k;
    A_op = -kappa^2;
    Mat = [k1/mu1 * A_op,          1/mu1;
           -k2*k1/mu1 * A_op,     -k2*(1/mu1 + 1/mu2)];

    [V_mat, D_mat] = eig(Mat);
    lambdas = diag(D_mat);
    omegas = real(sqrt(-real(lambdas)));
    [omegas, sort_idx] = sort(omegas);
    V_mat = V_mat(:, sort_idx);
    omega1 = omegas(1);
    omega2 = omegas(2);
    v1 = V_mat(:,1);
    v2 = V_mat(:,2);


    C_sin = V_mat \ [1; 0];

    B_cos = (V_mat \ [-kappa; 0]) ./ omegas;

    u_sin = C_sin(1)*v1(1)*cos(omega1*t) + C_sin(2)*v2(1)*cos(omega2*t);
    u_cos = B_cos(1)*v1(1)*sin(omega1*t) + B_cos(2)*v2(1)*sin(omega2*t);

    u_mode = amp * (u_sin .* sin(kappa * x) + u_cos .* cos(kappa * x));
end

u_exact_fun = @(t, x) exact_mode(t, x, k_wave1, amp1, mu1, mu2, k1, k2) + exact_mode(t, x, k_wave2, amp2, mu1, mu2, k1, k2);

%% === PostScript Image ===
figure('Position',[80 80 1480 760],'Color','white');
hold on;

%plot(x, u_fd_exact(plot_time),   'b-', 'LineWidth',2.5, 'DisplayName','FD');
%plot(x, u_num_exact(plot_time),  'r-', 'LineWidth',2.5, 'DisplayName','Numerov');
plot(x, u_dft_exact(plot_time),  'm-', 'LineWidth',2.5, 'DisplayName','DFT');
plot(x, u_exact_fun(plot_time, x), 'k--', 'LineWidth',2.5, 'DisplayName','Exact');

xlabel('x');
ylabel('u(x,t)');
title(sprintf('DFT | t = %s', formatThousands(plot_time)));

ax = gca;
ax.XTick = -pi : pi/2 : pi;
ax.XTickLabel = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
ax.YTick = -2 : 1 : 2;
ax.YTickLabel = {'-2', '-1','0', '1', '2'};

ax.FontSize = 25;
ax.LineWidth = 1.2;

ylim([-2 2]);
xlim([-pi pi]);
grid on;
%legend('Location','northwest','FontSize',20);

print('-depsc2', '-r300', 'two_periodic_same_meta_DFT.eps');

%% Video
video_filename = fullfile('videos', 'two_periodic_same_meta.avi');
vid = VideoWriter(video_filename, 'Uncompressed AVI');
vid.FrameRate = 10;
open(vid);

figure('Position',[80 80 1480 760],'Color','white');

for i = 1:length(tspan)
    t = tspan(i);
    
    clf; hold on;
    plot(x, u_fd_exact(t),   'b-', 'LineWidth',3.0, 'DisplayName','FD');
    %plot(x, u_num_exact(t),  'r-', 'LineWidth',3.0, 'DisplayName','Numerov ');
    %plot(x, u_dft_exact(t),  'm-', 'LineWidth',3.0, 'DisplayName','DFT ');
    plot(x, u_exact_fun(t, x), 'k--', 'LineWidth',3.0, 'DisplayName','Exact');
    
    xlabel('x');
    ylabel('u(x,t)');
    title(sprintf('Exact | t = %s', formatThousands(t)));

    ax = gca;
    ax.XTick = -pi : pi/2 : pi;
    ax.XTickLabel = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
    ax.YTick = -2 : 1 : 2;
    ax.YTickLabel = {'-2', '-1','0', '1', '2'};
    
    ax.FontSize = 25;
    ax.LineWidth = 1.2;
    
    ylim([-2 2]);
    xlim([-pi pi]);
    grid on;
    %legend('Location','northwest');
    hold off;
    
    frame = getframe(gcf);
    writeVideo(vid, frame);
end

close(vid);