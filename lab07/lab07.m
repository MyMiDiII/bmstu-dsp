function lab07()
  step = 0.01;
  half_range = 5;
  x = -half_range:step:half_range;

  A = 1;
  sigma1 = 1;
  sigma2 = sqrt(2);

  u1 = gauss(x, sigma1, A);
  u2 = gauss(x, sigma2, A);

  v1 = fft(u1);
  v2 = fft(u2);

  error = 0.1;
  xiD = unifrnd(-error, error, 1, length(x));
  xiE = unifrnd(-error, error, 1, length(x));
  delta = std(xiD);
  epsilon = std(xiE);

  filt = tikhonov_filter(v1, v2, step, 2 * half_range, delta, epsilon);
  result = abs(ifft(fft(u2+xiE) .* filt));

  graph(1, 1, 1, x, [u1+xiD; u2+xiE; result],
        [
         "Поврежденный сигнал 1";
         "Поврежденный сигнал 2";
         "Фильтр сигнала"],
        "Регулиризация Тихонова");
end

function u = gauss(x, sigma, A)
  u = A * exp(- (x / sigma).^2);
end

function rho = residual(alpha, v1, v2, deltaX, T, delta, epsilon)
  N = length(v1);
  m = 0:N-1;
  coef = deltaX / N;
  const = 1 + (2 * pi * m / T).^2;
  denominator = (abs(v2).^2 * deltaX^2 + alpha .* const).^2;
  gamma = coef * sum(abs(v2).^2 * deltaX^2 .* abs(v1).^2 .* const ./ denominator, 2);
  beta = coef * sum(alpha.^2 * const .* abs(v1).^2 ./ denominator, 2);
  rho = beta - (delta - epsilon * sqrt(gamma))^2;
end

function h = tikhonov_filter(v1, v2, deltaX, T, delta, epsilon)
  N = length(v1);
  m = 0:N-1;
  coef = deltaX / N;
  const = 1 + (2 * pi * m / T).^2;
  rho = @(alpha) residual(alpha, v1, v2, deltaX, T, delta, epsilon);
  alpha = fzero(rho, [0, 1]);
  h = 0:length(v1)-1;
  for k = 1:length(h)
    h(k) = coef * sum(exp(2 * pi * I * k .* m / length(v1)) .* v1 .* conj(v2)
           ./ (abs(v2).^2 .* deltaX^2 + alpha * const), 2);
  endfor
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
