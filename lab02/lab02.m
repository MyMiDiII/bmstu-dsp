function lab02()
% Общее
N = 200;
interval = 5;

% Импульсы
x = linspace(-interval,interval,N);
% Прямоугольный
L = 2;
u_rect = zeros(size(x));
u_rect(abs(x / L) <= 1) = 1;

% Гаусса
sigma = 0.5;
u_gauss = exp(- x.^2 / sigma.^2);

n = 0:length(x)-1;

% Дискретное преобразование Фурье
disp("ДПФ: Прямоугольный");
tic();
twins_dft_v_rect = dft(u_rect);
toc();
dft_v_rect = fftshift(twins_dft_v_rect);

disp("ДПФ: Гаусс");
tic();
twins_dft_v_gauss = dft(u_gauss);
toc();
dft_v_gauss = fftshift(twins_dft_v_gauss);

abs_twins_dft_v_rect = abs(twins_dft_v_rect) / length(n);
abs_dft_v_rect = abs(dft_v_rect) / length(n);
abs_twins_dft_v_gauss = abs(twins_dft_v_gauss) / length(n);
abs_dft_v_gauss = abs(dft_v_gauss) / length(n);

% Быстрое преобразование Фурье
% Дискретное преобразование Фурье
disp("БПФ: Прямоугольный");
tic();
twins_fft_v_rect = fft(u_rect);
toc();
fft_v_rect = fftshift(twins_fft_v_rect);

disp("БПФ: Гаусс");
tic();
twins_fft_v_gauss = fft(u_gauss);
toc();
fft_v_gauss = fftshift(twins_fft_v_gauss);

abs_twins_fft_v_rect = abs(twins_fft_v_rect) / length(n);
abs_fft_v_rect = abs(fft_v_rect) / length(n);
abs_twins_fft_v_gauss = abs(twins_fft_v_gauss) / length(n);
abs_fft_v_gauss = abs(fft_v_gauss) / length(n);

figure;

subplot(2, 2, 1);
title("ДПФ: Прямоугольный");
hold on;
grid on;
plot(n, abs_twins_dft_v_rect);
plot(n, abs_dft_v_rect);
legend("С эффектом Близнецов", "Без эффекта Близнецов");

subplot(2, 2, 2);
title("ДПФ: Гаусс");
hold on;
grid on;
plot(n, abs_twins_fft_v_gauss);
plot(n, abs_fft_v_gauss);
legend("С эффектом Близнецов", "Без эффекта Близнецов");

subplot(2, 2, 3);
title("БПФ: Прямоугольный");
hold on;
grid on;
plot(n, abs_twins_fft_v_rect);
plot(n, abs_fft_v_rect);
legend("С эффектом Близнецов", "Без эффекта Близнецов");

subplot(2, 2, 4);
title("БПФ: Гаусс");
hold on;
grid on;
plot(n, abs_twins_fft_v_gauss);
plot(n, abs_fft_v_gauss);
legend("С эффектом Близнецов", "Без эффекта Близнецов");
end


function [y] = dft(x)
  y = zeros(1, length(x));
  coef = -2 * pi * sqrt(-1) / length(x);
  for k = 1:length(y)
    for n = 0:length(x)-1
      y(k) = y(k) + x(n+1) * exp(coef * k * n);
    endfor
  endfor
end

% 19
% u_rect_no_twins = zeros(size(x));
% for i = 1:length(x)
%   u_rect_no_twins(i) = u_rect(i) * (-1).^x(i);
% endfor
% 25
% dft_v_rect = dft(u_rect_no_twins);
