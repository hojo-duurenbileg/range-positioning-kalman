function [ x ] = NoVelocityKalmanFilterLS( x_o, z_meas, A, Sp, deltaT )

%% x

n = size(z_meas, 1) * 4;                                               % number of rows 
m = 2;                                                                 % number of columns

H_k = zeros(4, 2);                                                     % allocate H matrix

x = zeros(size(z_meas, 1), 2)';

% first run parametric least-squares to get x_1
[x_o, P_k_minus, aposteriori_k, ~, x_k_minus] = ParaLS(x_o, z_meas(1, 1:4), A);

x(:,1) = x_o;                                                       

Q_k = [Sp*deltaT 0; 0 Sp*deltaT];

P_k_minus = P_k_minus + Q_k;

z_meas = reshape(z_meas', [1, n])';                                    % vectorize measurements

i = 5;
t = 2;
while i < n
    r_o = 0;
    r_o(1) = sqrt( (A(1).x - x_o(1))^2 + (A(1).y - x_o(2))^2);         % distance to A1
    r_o(2) = sqrt( (A(2).x - x_o(1))^2 + (A(2).y - x_o(2))^2);         % distance to A2
    r_o(3) = sqrt( (A(3).x - x_o(1))^2 + (A(3).y - x_o(2))^2);         % distance to A3
    r_o(4) = sqrt( (A(4).x - x_o(1))^2 + (A(4).y - x_o(2))^2);         % distance to A4
    r_o = r_o';                                                                    % transpose r_o
    
    H_k(1,1) = -(A(1).x - x_o(1)) / r_o(1);                    % H2 matrix
    H_k(1,2) = -(A(1).y - x_o(2)) / r_o(1);                    % H2 matrix

    H_k(2,1) = -(A(2).x - x_o(1)) / r_o(2);                    % H2 matrix
    H_k(2,2) = -(A(2).y - x_o(2)) / r_o(2);                    % H2 matrix

    H_k(3,1) = -(A(3).x - x_o(1)) / r_o(3);                    % H2 matrix
    H_k(3,2) = -(A(3).y - x_o(2)) / r_o(3);                    % H2 matrix

    H_k(4,1) = -(A(4).x - x_o(1)) / r_o(4);                    % H2 matrix
    H_k(4,2) = -(A(4).y - x_o(2)) / r_o(4);                    % H2 matrix

    R_k = eye(4);                                            % R2 - Covariance of measurements 
    
    dz = z_meas(i:i+3, 1) - r_o;                             % misclosure vector  

    N = R_k + H_k*P_k_minus*H_k';                            % N    
    N_ = mldivide(N, eye(size(N)));                          % N inverse

    k_k = P_k_minus * H_k' *( N_ );                          % Gain matrix

    x_k_plus = x_k_minus + k_k * (dz - H_k * x_k_minus);    % Corrections to the state vector
    
    x_o = x_o + x_k_plus;
    
    P_k_plus = (eye(2) - k_k * H_k) * P_k_minus;             % P2
    
    i = i + 4;                                               % Move to next epoch
    
    x_kplus1_minus = x_k_plus;
    
    P_kplus1_minus = P_k_plus + Q_k;
    
    x_k_minus = x_kplus1_minus;
    
    P_k_minus = P_kplus1_minus;
    
    x(:, t) = x_o;
    
    t = t + 1;
end

% final residuals
% r_o = 0;
% r_o(1) = sqrt( (A(1).x - x(1))^2 + (A(1).y - x(2))^2);         % distance to A1
% r_o(2) = sqrt( (A(2).x - x(1))^2 + (A(2).y - x(2))^2);         % distance to A2
% r_o(3) = sqrt( (A(3).x - x(1))^2 + (A(3).y - x(2))^2);         % distance to A3
% r_o(4) = sqrt( (A(4).x - x(1))^2 + (A(4).y - x(2))^2);         % distance to A4
% r_o = repmat(r_o, 1, (n / 4))';                                % match size of z_meas
% 
% r = r_o - z_meas;                                              % residuals
% 
% aposteriori = r' * eye(n) * r / (n - m);                       % A posteriori variance factor
% 
% N = H_k'*eye(4)*H_k;                                           % N
% N_ = mldivide(N, eye(size(N)));                                % N inverse
% 
% P = aposteriori * N_;

end
