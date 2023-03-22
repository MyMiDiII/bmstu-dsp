function lab04()
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

  % Фильтры низких частот
  lbf = low_butterworth_filter(6, 20);
  lgf = low_gauss_filter(4, 20);

  u_imp_lbf = filtfilt(lbf, 1, u_imp);
  u_gauss_lbf = filtfilt(lbf, 1 ,u_gauss);

  u_imp_lgf = filtfilt(lgf, 1, u_imp);
  u_gauss_lgf = filtfilt(lgf, 1 ,u_gauss);

  figure("name", "Низкие частоты");
  graph(2, 4, 1:2, x, [u_gauss; u_imp; u0],
        [
         "Сигнал с помехой Гаусса";
         "Сигнал с импульсной помехой";
         "Идеальный сигнал"],
        "Исходные сигналы");

  graph(2, 4, 3:4, x, [u_gauss_lbf; u_imp_lbf; u0],
        [
         "Сигнал с помехой Гаусса";
         "Сигнал с импульсной помехой";
         "Идеальный сигнал"],
        "Фильтр Баттеруорта");

  graph(2, 4, 6:7, x, [u_gauss_lgf; u_imp_lgf; u0],
        [
         "Сигнал с помехой Гаусса";
         "Сигнал с импульсной помехой";
         "Идеальный сигнал"],
        "Фильтр Гаусса");

  % Фильтры высоких частот
  hbf = high_butterworth_filter(6, 20);
  hgf = high_gauss_filter(4, 20);

  u_imp_hbf = u_imp - filtfilt(hbf, 1, u_imp);
  u_gauss_hbf = u_gauss - filtfilt(hbf, 1 ,u_gauss);

  u_imp_hgf = u_imp - filtfilt(hgf, 1, u_imp);
  u_gauss_hgf = u_gauss - filtfilt(hgf, 1 ,u_gauss);

  figure("name", "Высокие частоты");
  graph(2, 4, 1:2, x, [u_gauss; u_imp; u0],
        [
         "Сигнал с помехой Гаусса";
         "Сигнал с импульсной помехой";
         "Идеальный сигнал"],
        "Исходные сигналы");

  graph(2, 4, 3:4, x, [u_gauss_hbf; u_imp_hbf; u0],
        [
         "Сигнал с помехой Гаусса";
         "Сигнал с импульсной помехой";
         "Идеальный сигнал"],
        "Фильтр Баттеруорта");

  graph(2, 4, 6:7, x, [u_gauss_hgf; u_imp_hgf; u0],
        [
         "Сигнал с помехой Гаусса";
         "Сигнал с импульсной помехой";
         "Идеальный сигнал"],
        "Фильтр Гаусса");
end

function y = low_butterworth_filter(F, size)
  f = linspace(-size/2,size/2,size);
  y = 1./(1 + (f./F).^4);
  y = y / sum(y);
end

function y = low_gauss_filter(sigma, size)
  f = linspace(-size/2,size/2,size);
  y = exp(-f.^2/(2*sigma^2));
  y = y / sum(y);
end

function y = high_butterworth_filter(F, size)
  f = linspace(-size/2,size/2,size);
  y = 1./(1 + (F./f).^4);
  y = y / sum(y);
end

function y = high_gauss_filter(sigma, size)
  f = linspace(-size/2,size/2,size);
  y = 1 - exp(-f.^2/(2*sigma^2));
  y = y / sum(y);
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
  endfor
  legend(lbls);
end
