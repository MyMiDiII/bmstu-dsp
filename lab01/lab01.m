% Исходные сигналы
step = 0.005;
% Гаусса
sigma_num = 3;
sigma = 2;
x_gauss = -sigma_num * sigma:step:sigma_num * sigma;
u_gauss = exp(- x_gauss.^2 / sigma.^2);
% Прямоугольный
indent = 8;
L = 2;
x_rect = -L - indent:step:L + indent;
u_rect = zeros(size(x_rect))
u_rect(abs(x_rect / L) <= 1) = 1;

% Дикретизация
delta_x = 1;
% Гаусс
x_gauss_descr = -sigma_num * sigma:delta_x:sigma_num * sigma;
u_gauss_descr = exp(- x_gauss_descr.^2 / sigma.^2);
% Прямоугольный
x_rect_descr = -L - indent:delta_x:L + indent;
u_rect_descr = zeros(size(x_rect_descr))
u_rect_descr(abs(x_rect_descr / L) <= 1) = 1;

% Восстановление сигнала
function [restored]=restore_signal(x_restored, x_discr, discr, delta_x)
  restored = zeros(size(x_restored));
  disp("size restored")
  disp(size(restored))
  for i = 1:length(x_restored)
    for j = 1:length(discr)
      % sinc(x) = sin(pi * x) / (pi * x)
      restored(i) = restored(i) + discr(j) * sinc(delta_x * (x_restored(i) - x_discr(j)));
    endfor
  endfor
end

u_gauss_restored = restore_signal(x_gauss, x_gauss_descr, u_gauss_descr, delta_x);
u_gauss_restored;
u_rect_restored = restore_signal(x_rect, x_rect_descr, u_rect_descr, delta_x);

% Графики
f = figure();
set(0, "defaultlinelinewidth", 3, "defaultaxesfontsize", 16);

% Гаусс
subplot(2, 1, 2);
hold on;
grid on;
plot(x_gauss, u_gauss);
plot(x_gauss_descr, u_gauss_descr, ".r", "markersize", 15);
plot(x_gauss, u_gauss_restored, "g");
legend("Исходный", "Дискретный", "Восстановленный");

% Прямоугольный
subplot(2, 1, 1);
hold on;
grid on;
plot(x_rect, u_rect);
plot(x_rect_descr, u_rect_descr, ".r", "markersize", 15);
plot(x_rect, u_rect_restored, "g");
legend("Исходный", "Дискретный", "Восстановленный");
