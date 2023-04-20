function lab08()
  epsv = 0.05;

  step = 0.05;
  half_range = 5;
  x = -half_range:step:half_range;
  N = length(x);

  sigma = 2;
  amplitude = 1;
  u0 = gauss(x, sigma, amplitude);

  xi_imp = impulse_noise(length(x), 10, amplitude);
  u_base = u0 + xi_imp;

  u_med_sm = smooth_imp(@med, u_base, epsv);
  u_mean_sm = smooth_imp(@mean, u_base, epsv);

  graph(1, 2, 1, x, [u_base; u_med_sm],
        ["Исходный сигнал",
         "Сглаженный сигнал"],
        "MED-функция фильтрации");

  graph(1, 2, 2, x, [u_base; u_mean_sm],
        ["Исходный сигнал",
         "Сглаженный сигнал"],
        "MEAN-функция фильтрации");
end

function u = smooth_imp(smth, v, epsv)
  u = v;

  for i = 1:1:length(u)
    smth_val = smth(v, i);
    if (abs(u(i) - smth_val) > epsv)
      u(i) = smth_val;
    endif
  endfor
end

function y = med(u, i)
  if (i == 1)
    m = [u(i), u(i+1)];
  elseif (i == length(u))
    m = [u(i), u(i - 1)];
  else
    m = [u(i - 1), u(i), u(i + 1)];
  endif
  m = sort(m);
  y = m(ceil((length(m) + 1) / 2));
end

function y = mean(u, i)
  step = 3;
  i_from = i - step;
  i_to = i + step;

  s = 0;
  num = 0;
  for k = i_from:1:i_to
    if (k > 0 && k <= length(u))
      s += u(k);
      num += 1;
    endif
  endfor

  y = s / num;
end

function noise = impulse_noise(points_num, impulses_num, sig_amplitude)
  begin = randi([fix(points_num / 6), fix(points_num / 5)], 1, 1);
  steps = randi([40, 100], impulses_num - 1, 1);
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

function u = gauss(x, sigma, A)
  u = A * exp(- (x / sigma).^2);
end

function graph(rows, columns, num, x, us, lbls, title_string)
  set(0, "defaultlinelinewidth", 3, "defaultaxesfontsize", 16);
  subplot(rows, columns, num);
  title(title_string);
  hold on;
  grid on;
  for signal = us'
    plot(x, signal);
  endfor
  legend(lbls);
end
