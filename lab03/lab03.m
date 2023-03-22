function lab03()
  step = 0.01;
  half_range = 5;
  x = -half_range:step:half_range;

  rect1 = rect(x, 1, 1);
  rect2 = rect(x, 2, 2);

  conv_rr = convolution(rect1, rect2, step);
  convolution_graph(2, 4, 1:2,
                    x, [rect1; rect2; conv_rr],
                    ["Прямоугольный 1"; "Прямоугольный 2"; "Свертка"],
                    "Свертка двух прямоугольных сигналов");

  gauss1 = gauss(x, 2, 1);
  gauss2 = gauss(x, 1, 2);

  conv_gg = convolution(gauss1, gauss2, step);
  convolution_graph(2, 4, 3:4,
                    x, [gauss1; gauss2; conv_gg],
                    ["Гаусс 1"; "Гаусс 2"; "Свертка"],
                    "Свертка двух сигналов Гаусса");

  conv_gg = convolution(rect2, gauss2, step);
  convolution_graph(2, 4, 6:7,
                    x, [rect2; gauss2; conv_gg],
                    ["Прямоугольный 2"; "Гаусс 2"; "Свертка"],
                    "Свертка двух прямоугольного сигнала и сигнала Гаусса");
end

function w = convolution(sig1, sig2, k)
  z = zeros(size(sig1));
  n = length(z);
  tilda_sig1 = [sig1 z];
  tilda_sig2 = [sig2 z];

  w = ifft(fft(tilda_sig1).*fft(tilda_sig2)) * k;
  indent = fix(n / 2);
  w = w(n - indent:n + indent);
end

function u = rect(x, L, A)
  u = zeros(size(x));
  u(abs(x / L) <= 1) = A;
end

function u = gauss(x, sigma, A)
  u = A * exp(- (x / sigma).^2);
end

function convolution_graph(rows, columns, num, x, us, lbls, title_string)
  subplot(rows, columns, num);
  title(title_string);
  hold on;
  grid on;
  for signal = us'
    plot(x, signal);
  endfor
  legend(lbls);
end
