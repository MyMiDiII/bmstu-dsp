function lab06()
  step = 0.01;
  half_range = 5;
  x = -half_range:step:half_range;

  sigma = 2;
  amplitude = 1;
  u0 = gauss(x, sigma, amplitude);

  xi_imp = impulse_noise(length(x), 10, amplitude);
  u_imp = u0 + xi_imp;

  xi_gauss = gauss_noise(length(x));
  u_gauss = u0 + xi_gauss;

  % Фильтр Винера
  wiener_imp = wiener_filter(u_imp, xi_imp);
  wiener_gauss = wiener_filter(u_gauss, xi_gauss);

  u_wiener_imp = filtering(u_imp, wiener_imp);
  u_wiener_gauss = filtering(u_gauss, wiener_gauss);

  figure("name", "Фильтр Винера");

  graph(2, 1, 1, x, [u_imp; u_wiener_imp],
        ["Сигнал с импульсной помехой";
         "Отфильтрованный сигнал";],
        "Импульсная помеха");

  graph(2, 1, 2, x, [u_gauss; u_wiener_gauss],
        ["Сигнал с помехой Гаусса";
         "Отфильтрованный сигнал";],
        "Помеха Гаусса");

  % Режекторный фильтр
  rejection_imp = rejection_filter(u_imp, xi_imp);
  rejection_gauss = rejection_filter(u_gauss, xi_gauss);

  u_rejection_imp = filtering(u_imp, rejection_imp);
  u_rejection_gauss = filtering(u_gauss, rejection_gauss);

  figure("name", "Режекторный фильтр");

  graph(2, 1, 1, x, [u_imp; u_rejection_imp],
        ["Сигнал с импульсной помехой";
         "Отфильтрованный сигнал";],
        "Импульсная помеха");

  graph(2, 1, 2, x, [u_gauss; u_rejection_gauss],
        ["Сигнал с помехой Гаусса";
         "Отфильтрованный сигнал";],
        "Помеха Гаусса");
end

function y = wiener_filter(u, xi)
  y = 1 - (fft(xi)./fft(u)).^2;
end

function y = rejection_filter(u, xi)
  u = fft(u);
  xi = fft(xi);
  y = zeros(size(u));
  for i = 1:length(u)
    if abs(u(i)) - abs(xi(i)) > 1
      y(i) = 1;
    else
      y(i) = 0;
    endif
  endfor
end

function y = filtering(signal, filter)
  y = ifft(fft(signal).*filter);
end

function noise = impulse_noise(points_num, impulses_num, sig_amplitude)
  begin = randi([fix(points_num / 6), fix(points_num / 3)], 1, 1);
  steps = randi([30, 100], impulses_num - 1, 1);
  positions = zeros(1, impulses_num);
  positions(1) = mod(begin, points_num) + 1;
  for i = 2:length(positions)
    positions(i) = mod(positions(i - 1) + steps(i - 1), points_num) + 1;
  endfor

  noise = zeros(1, points_num);
  amplitudes = sig_amplitude * randi([15, 20], impulses_num, 1) / 100;

  for i = 1:impulses_num;
    noise(positions(i)) = amplitudes(i);
  endfor
end

function noise = gauss_noise(points_num)
  mu = 0;
  sigma = 0.05;
  noise = normrnd(mu,sigma,[1 points_num]);
end

function u = gauss(x, sigma, A)
  u = A * exp(- (x / sigma).^2);
end

function graph(rows, columns, num, x, us, lbls, title_string)
  subplot(rows, columns, num);
  title(title_string);
  hold on;
  grid on;
  for signal = us'
    plot(x, signal);
    axis([-6,6.1,-0.2,1.41], "ticy");
    set(gca,'YTick',[-0.2:0.2:1.6]);
    set(gca,'XTick',[-6:2:6]);
  endfor
  legend(lbls);
end
