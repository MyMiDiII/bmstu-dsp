% Исходные сигналы
step = 0.001;
% Гаусса
sigma = 3;
x_gauss = -3 * sigma:step:3 * sigma;
u_gauss = exp(- x_gauss.^2 / sigma.^2);
% Прямоугольный
L = 3;
x_rect = -L - 1:step:L + 1;
u_rect = zeros(size(x_rect))
u_rect(abs(x_rect / L) <= 1) = 1;

% Графики
f = figure();
set(0, "defaultlinelinewidth", 3, "defaultaxesfontsize", 16);

% Гаусс
subplot(1, 2, 1);
hold on;
plot(x_gauss, u_gauss);

% Прямоугольный
subplot(1, 2, 2);
hold on;
plot(x_rect, u_rect);
