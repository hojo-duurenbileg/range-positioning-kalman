function [ ax ] = EllipsePlot(x, a, b, beta)

% compute points corresponding to axis-oriented ellipse
t = -0.01:0.01:(2*pi);

X = x(1) + a * cos(t) * sin(beta) - b * sin(t) * cos(beta);
Y = x(2) + a * cos(t) * cos(beta) + b * sin(t) * sin(beta);

plot(X, Y, '-');

ax = gca;

end

