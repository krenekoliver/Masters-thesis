clear all; close all; clc;

function str = formatThousands(num)
    str = num2str(round(num), '%d');
    str = regexprep(str, '(\d)(?=(\d{3})+$)', '$1 ');
end

N = 128;
I = pi;
dx = I / N;
x = linspace(0, I - dx, N)';

c = 1;
k = 4;
omega = c * k;

T_final = 30;

t_postscript = 5000017.495;
postscript_filename = sprintf('pic_one_periodic_DFT.eps');

%% initial conditions
u0 = sin(k * x);
v0 = -omega * cos(k * x);

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
omega_k = c * abs(kx);

u0_hat = fft(u0);  
v0_hat = fft(v0); 

u_fft_fun = @(t) real(ifft( u0_hat .* cos(omega_k * t) + (v0_hat ./ max(omega_k, 1e-12)) .* sin(omega_k * t) ));

%% Exact solution
u_exact_fun = @(t) sin(k * x - omega * t);

%% Video
video_filename = fullfile('videos', 'one_periodic.avi');
vid = VideoWriter(video_filename, 'Uncompressed AVI');
vid.FrameRate = 40;
open(vid);

t_vec = linspace(0, T_final, 2000);

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
    plot(x, u_an,   'k--', 'LineWidth', 2.5, 'DisplayName', 'Exact solution');

    xlabel('x');
    ylabel('u(x,t)');
    title(sprintf('DFT | t = %s', formatThousands(t)));
    
    ax = gca;
    ax.XTick = 0 : pi/4 : I;
    ax.XTickLabel = {'0', '\pi/4','\pi/2', '3\pi/4', '\pi'};
    ax.YTick = -1 : 0.5 : 1;
    ax.YTickLabel = {'-1', '-1/2','0', '1/2', '1'};
    
    ax.FontSize = 25;
    ax.LineWidth = 1.2;

    % title(sprintf('Exact | t = %s', formatThousands(t)), 'FontSize', 27);
    % % legend('Location', 'best');
    ylim([-1 1]);
    xlim([0 I]);

    grid on;
    hold off;
    
    frame = getframe(gcf);
    writeVideo(vid, frame);
    
    % === PostScript (EPS) ===
    if t == t_postscript
        set(gcf, 'PaperPositionMode', 'auto');
        set(gcf, 'Renderer', 'Painters');
        print(gcf, postscript_filename, '-depsc');
    end
end

close(vid);