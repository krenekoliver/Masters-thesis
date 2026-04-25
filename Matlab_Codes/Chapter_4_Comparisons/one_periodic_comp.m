clear all; close all; clc;

function str = formatThousands(num)
    str = num2str(round(num), '%d');
    str = regexprep(str, '(\d)(?=(\d{3})+$)', '$1 ');
end

k_wave   = 4;
m2       = 0.1;
plot_time = 0;

N = 256;
I = 2*pi;
k1 = 5;
k2 = 50;
m1 = 1;

x = I/N * (-N/2 : N/2-1)';

c_limit = I * sqrt(k1 / m1);

%% initial conditions
u0 = sin(k_wave * x);
v0 = - c_limit * k_wave * cos(k_wave * x);
a0 = zeros(N,1);

y0_exact = [u0; a0];
ydot0    = [v0; zeros(N,1)];

tspan = linspace(0, 300, 1000);

%% Wave equation (DFT)
kx = (2*pi/I) * [0:(N/2-1), -N/2:-1]';
omega_classic = c_limit * abs(kx);

u0_hat = fft(u0);
v0_hat = fft(v0);

u_classic_fun = @(t) real(ifft( u0_hat .* cos(omega_classic*t) + (v0_hat ./ max(omega_classic,1e-12)) .* sin(omega_classic*t) ));

%% Rate-type metamaterial (DFT)
D1 = diag(1./m1); 
D2 = diag(1./m2);

k_vec = (2*pi/I) * [0:(N/2-1), -N/2:-1]';
Lambda = diag(-k_vec.^2);
F = fft(eye(N));
A = real(ifft(Lambda * F));

TopLeft     = (k1 * I^2) * (D1 * A);
TopRight    = D1;
BottomLeft  = -k2 * (k1 * I^2) * (D1 * A);
BottomRight = -k2 * (D1 + D2);

M = [ TopLeft,     eye(N)*TopRight;
      BottomLeft,  eye(N)*BottomRight ];

[V, D] = eig(full(M));
omega_meta = sqrt(-diag(D) + 1e-12i);

coeff_cos = V \ y0_exact;
coeff_sin = (V \ ydot0) ./ omega_meta;

u_meta_fun = @(t) real(V(1:N,:) * (coeff_cos .* cos(omega_meta*t) + coeff_sin .* sin(omega_meta*t)));

%% PostScript
figure('Position',[80 80 1480 760],'Color','white');
hold on;
plot(x, u_classic_fun(plot_time), 'b-', 'LineWidth',2.5);
plot(x, u_meta_fun(plot_time), 'r-', 'LineWidth',2.5);

xlabel('x'); ylabel('u(x,t)');
title(sprintf('Initial | t = %s', formatThousands(plot_time)));
ax = gca;
ax.XTick = -pi : pi/2 : pi;
ax.XTickLabel = {'-\pi','-\pi/2','0','\pi/2','\pi'};
ax.YTick = -1 : 0.5 : 1;
ax.YTickLabel = {'-1','-1/2','0','1/2','1'};
ax.FontSize = 25;
ylim([-1 1]);
xlim([-pi pi]);
grid on;
%legend('Location','northwest');

print('-depsc2', '-r300', sprintf('comparison_periodic_one_initial.eps'));

%% Video
video_filename = fullfile('videos', sprintf('comparison_m.avi'));
vid = VideoWriter(video_filename, 'Uncompressed AVI');
vid.FrameRate = 30;
open(vid);

figure('Position',[80 80 1480 760],'Color','white');
for i = 1:length(tspan)
    t = tspan(i);
    clf; hold on;
    
    plot(x, u_classic_fun(t), 'b-', 'LineWidth',2.5);
    plot(x, u_meta_fun(t),    'r-', 'LineWidth',2.5);
    
    xlabel('x'); ylabel('u(x,t)');
    title(sprintf('Periodic | t = %s', formatThousands(t)));
    ax = gca;
    ax.XTick = -pi : pi/2 : pi;
    ax.XTickLabel = {'-\pi','-\pi/2','0','\pi/2','\pi'};
    ax.FontSize = 25;
    ylim([-1 1]);
    xlim([-pi pi]);
    grid on;
    %legend('Classical wave', 'Metamaterial', 'Location','northwest');
    
    frame = getframe(gcf);
    writeVideo(vid, frame);
end
close(vid);
