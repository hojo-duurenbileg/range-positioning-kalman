%% Setup
clear;
close all;
clc;

raw_data = load('data/ranges.txt');    % read data

A(1).x = 0;                     % Anchor point 1
A(1).y = 0;
A(2).x = 100;                   % Anchor point 1
A(2).y = 0;
A(3).x = 100;                   % Anchor point 1
A(3).y = 100;
A(4).x = 0;                     % Anchor point 1
A(4).y = 100;

x_o = [50 50]';                  % POE [x, y]

z_meas = raw_data(:, 2:5);      % range measurements

%% Task 1. a)

x_ee = zeros(150, 2)';              % allocate array for coordinate estimates

for i = 1:size(z_meas)
    
    x_ee(:,i) = ParaLS(x_o, z_meas(i, :), A);      % Parametric least-squares for each epoch
    x_o = x_ee(:,i);                               % Use last estimate as initial estimate
    
    if i == 7
        [x_3, P_3] = ParaLS(x_o, z_meas(i, :), A);
    end
    
end

% Plot
f = figure(1);
f.Name = 'Plot of (x, y) at each epoch';
set(gcf,'color', 'w');

s = scatter(x_ee(1,:),x_ee(2,:));
s.Marker = '.';
title('Plot of (x, y) at each epoch');
xlabel('x (m)');
ylabel('y (m)');

ax = gca;
axis equal
ax.XLim = [-10 110];
ax.YLim = [-10 110];
ax.XGrid = 'on';
ax.YGrid = 'on';

%% Task 1. b)

z_meas = raw_data(1:51, 2:5);                       % Get stationary data

x_o = [50 50]';                                      % POE [x, y]

[x, P, aposteriori, r]= ParaLS(x_o, z_meas, A);        % Batch least-squares

[a,b,beta] = ErrorEllipses( P );

f = figure(2);
f.Name = 'Plot of error ellipse for BATCH estimate of x';

ax = EllipsePlot(x, a, b, beta);
title('Plot of error ellipse for BATCH estimate of x');
xlabel('x (m)');
ylabel('y (m)');
axis equal
ax.XGrid = 'on';
ax.YGrid = 'on';

%% Task 1. c)
f = figure(4);
f.Name = 'Plot of residuals for Anchor pt 2 BATCH';
set(gcf,'color', 'w');

plot(r(2:4:end));
title('Plot of residuals for Anchor pt 2 BATCH');
xlabel('epoch (s)');
ylabel('residuals (m)');
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.XLim = [0 51];

%% Task 1. d)

f = figure(3);
f.Name = 'Plot of epoch 7 error ellipse overlayed with stationary points';
set(gcf,'color', 'w');

s = scatter(x_ee(1,1:51),x_ee(2,1:51));
s.Marker = '.';
title('Plot of epoch 7 error ellipse overlayed with stationary points');
xlabel('x');
ylabel('y');

ax = gca;
xlabel('x (m)');
ylabel('y (m)');
axis equal;
ax.XGrid = 'on';
ax.YGrid = 'on';

[a_3,b_3,beta_3] = ErrorEllipses( P_3 );
hold on
EllipsePlot(x_3,a_3,b_3,beta_3);
hold off

hold on
EllipsePlot(x, a, b, beta);
hold off


%% Task 2. a)

z_meas = raw_data(1:51, 2:5);

x_o = [50 50]';

[x, P, aposteriori, r]= SummationLS(x_o, z_meas, A);

[a,b,beta] = ErrorEllipses( P );

f = figure(5);
f.Name = 'Plot of error ellipse for Summations estimate of x';
set(gcf,'color', 'w');

ax = EllipsePlot(x, a, b, beta);
title('Plot of error ellipse for Summations estimate of x');
xlabel('x (m)');
ylabel('y (m)');
axis equal;
ax.XGrid = 'on';
ax.YGrid = 'on';

%% Task 2. b)

z_meas = raw_data(1:51, 2:5);

x_o = [50 50]';

[x, P, aposteriori, r]= SequentialLS(x_o, z_meas, A);

[a,b,beta] = ErrorEllipses( P );

f = figure(6);
f.Name = 'Plot of error ellipse for Sequential estimate of x';
set(gcf,'color', 'w');

ax = EllipsePlot(x, a, b, beta);
title('Plot of error ellipse for Sequential estimate of x');
xlabel('x (m)');
ylabel('y (m)');
axis equal;
ax.XGrid = 'on';
ax.YGrid = 'on';

%% Task 3. a)

z_meas = raw_data(1:150, 2:5);

x_o = [50 50]';

[x_s, P, aposteriori, r]= IteratedSequentialLS(x_o, z_meas, A);

[a,b,beta] = ErrorEllipses( P );

% Plot
f = figure(7);
f.Name = 'SEQUENTIAL: Plot of (x, y) at each epoch';
set(gcf,'color', 'w');

s = scatter(x_s(1,:),x_s(2,:));
s.Marker = '.';
title('SEQUENTIAL: Plot of (x, y) at each epoch');
xlabel('x (m)');
ylabel('y (m)');

ax = gca;
axis equal;
ax.XLim = [0 1200];
ax.YLim = [0 300];
ax.XGrid = 'on';
ax.YGrid = 'on';

%% Task 3. b)

Sp = 10;
deltaT = 1;

z_meas = raw_data(1:150, 2:5);

x_o = [50 50]';

[x_k1] = NoVelocityKalmanFilterLS(x_o, z_meas, A, Sp, deltaT);

% Plot
f = figure(8);
f.Name = 'KALMAN without Velocity: Plot of (x, y) at each epoch Sp = 10';
set(gcf,'color', 'w');

s = plot(x_k1(1,:),x_k1(2,:));
s.LineWidth = 2;
title('KALMAN without Velocity: Plot of (x, y) at each epoch Sp = 10');
xlabel('x (m)');
ylabel('y (m)');

ax = gca;
axis equal;
ax.XLim = [45 110];
ax.YLim = [45 100];
ax.XGrid = 'on';
ax.YGrid = 'on';

%% Task 3. c)

Sp = 0.5;
Sv = 0.1;
deltaT = 1;

z_meas = raw_data(1:150, 2:5);

x_o = [50 50]';

[x_k] = KalmanFilterLS(x_o, z_meas, A, Sp, Sv, deltaT);

% Plot
f = figure(9);
f.Name  = 'KALMAN with Velocity: Plot of (x, y) at each epoch Sp = 0.5, Sv = 0.1';
set(gcf,'color', 'w');

s = plot(x_k(1,:),x_k(2,:));
s.LineWidth = 2;
title('KALMAN with Velocity: Plot of (x, y) at each epoch Sp = 0.5, Sv = 0.1');
xlabel('x (m)');
ylabel('y (m)');

ax = gca;
axis equal;
ax.XLim = [45 110];
ax.YLim = [45 100];
ax.XGrid = 'on';
ax.YGrid = 'on';

v = sqrt(x_k(3, :).^2 + x_k(4, :).^2);

% Plot
f = figure(10);
f.Name = 'Velocity plot at each epoch Sp = 0.5, Sv = 0.1';
set(gcf,'color', 'w');

s = plot(1:150, v);
s.LineWidth = 2;
title('Velocity plot at each epoch Sp = 0.5, Sv = 0.1');
xlabel('epoch (s)');
ylabel('velocity (m/s)');

ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';

