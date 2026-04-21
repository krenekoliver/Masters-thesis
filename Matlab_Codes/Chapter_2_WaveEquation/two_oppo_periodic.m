clear all; close all; clc;

function str = formatThousands(num)
    str = num2str(round(num), '%d');
    str = regexprep(str, '(\d)(?=(\d{3})+$)', '$1 ');
end

N = 128;            
I = 2*pi;                
dx = I / N;
x = linspace(0, I - dx, N)'; 
c = 1;
k = 5;
l = 6;

omega_k = c * k;
omega_l = c * l;

T_final = 5000001.7878;

t_postscript = 5000001.7878;
postscript_filename = sprintf('pic_two_periodic_oppo_DFT.eps');

%% initial conditions
u0 = sin(k * x) + sin(l * x);
v0 = -omega_k * cos(k * x) + omega_l * cos(l * x);

%% FD
K_fd = spdiags([ones(N,1), -2*ones(N,1), ones(N,1)], -1:1, N, N);
K_fd(1,N) = 1; K_fd(N,1) = 1;
A_fd = (c^2 / dx^2) * K_fd;

[V_fd, D_fd] = eig(full(A_fd));
lambda_fd = diag(D_fd);

u_fd_fun = @(t) real(V_fd * (cos(t*sqrt(-lambda_fd)).*(V_fd\u0) + (sin(t*sqrt(-lambda_fd))./sqrt(-lambda_fd)).*(V_fd\v0)));

%% Numerov
M = spdiags([ones(N,1), 10*ones(N,1), ones(N,1)], -1:1, N, N);
M(1,N) = 1; M(N,1) = 1;

K_num = spdiags([ones(N,1), -2*ones(N,1), ones(N,1)], -1:1, N, N);
K_num(1,N) = 1; K_num(N,1) = 1;

const = 12 / dx^2;
A_num = c^2 * (M \ (const * K_num));

[V_num, D_num] = eig(full(A_num));
lambda_num = diag(D_num);

u_num_fun = @(t) real(V_num * (cos(t*sqrt(-lambda_num)).*(V_num\u0) + (sin(t*sqrt(-lambda_num))./sqrt(-lambda_num)).*(V_num\v0)));

%% DFT
kx = [0:N/2-1, -N/2:-1]' * (2*pi / I);
omega_k_all = c * abs(kx);

u0_hat = fft(u0);
v0_hat = fft(v0);

u_fft_fun = @(t) real(ifft( u0_hat .* cos(omega_k_all * t) + (v0_hat ./ max(omega_k_all, 1e-12)) .* sin(omega_k_all * t) ));

%% Exact solution
u_exact_fun = @(t) sin(k * x - omega_k * t) + sin(l * x + omega_l * t);

%% Video
video_filename = fullfile('videos', 'two_opposite_periodic.avi');
vid = VideoWriter(video_filename, 'Uncompressed AVI');
vid.FrameRate = 40;
open(vid);

t_vec = linspace(0, T_final, 10);

figure('Position', [80 80 1480 760],'Color','white');

for i = 1:length(t_vec)
    t = t_vec(i);
    
    u_fd   = u_fd_fun(t);
    u_num  = u_num_fun(t);
    u_fft  = u_fft_fun(t);
    u_an   = u_exact_fun(t);
    
    figure(1); clf; hold on;
    
    % plot(x, u_fd,   'b-',  'LineWidth', 2.5, 'DisplayName', 'FD');
    % plot(x, u_num,  'r-',  'LineWidth', 2.5, 'DisplayName', 'Numerov');
     plot(x, u_fft,  'm-',  'LineWidth', 2.5, 'DisplayName', 'DFT');
    plot(x, u_an,   'k--', 'LineWidth', 2.5, 'DisplayName', 'Exact');
    
    xlabel('x');
    ylabel('u(x,t)');
    title(sprintf('DFT | t = %s', formatThousands(t)));
    
    ax = gca;
    ax.XTick = 0 : pi/2 : I;
    ax.XTickLabel = {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'};
    ax.YTick = -2 : 1 : 2;
    ax.YTickLabel = {'-2', '-1','0', '1', '2'};
    
    ax.FontSize = 25;
    ax.LineWidth = 1.2;

    % title(sprintf('Num | t = %s', formatThousands(t)), 'FontSize', 25);
    % % legend('Location', 'best');
    ylim([-2 2]);
    xlim([0 I]);
    grid on;
    hold off;
    
    frame = getframe(gcf);
    writeVideo(vid, frame);
    
    % === PostScript (EPS) ===
    if abs(t - t_postscript) < (t_vec(2)-t_vec(1))/2
        set(gcf, 'PaperPositionMode', 'auto');
        set(gcf, 'Renderer', 'Painters');
        print(gcf, postscript_filename, '-depsc');
    end
end

close(vid);