function [ x_final, P, aposteriori, r ] = SequentialLS( x_o, z_meas, A )
% MAKE SURE TO RETURN x_final or x based on which task you are using!!!!!

%% x

n = size(z_meas, 1) * 4;                                               % number of rows 
m = 2;                                                                 % number of columns

H = zeros(4, 2);                                                       % allocate H matrix
x = zeros(size(z_meas, 1), 2)';

% first run parametric least-squares to get x_1
[x_1, P_1, aposteriori_1, ~, dx_1] = ParaLS(x_o, z_meas(1, 1:4), A);

P_1 = (1 / aposteriori_1) * P_1; 

x_o = x_1;

x(:,1) = x_o;                                                       

z_meas = reshape(z_meas', [1, n])';                                    % vectorize measurements

i = 5;
t = 2;
while i < n
    r_o = 0;
    r_o(1) = sqrt( (A(1).x - x_1(1))^2 + (A(1).y - x_1(2))^2);         % distance to A1
    r_o(2) = sqrt( (A(2).x - x_1(1))^2 + (A(2).y - x_1(2))^2);         % distance to A2
    r_o(3) = sqrt( (A(3).x - x_1(1))^2 + (A(3).y - x_1(2))^2);         % distance to A3
    r_o(4) = sqrt( (A(4).x - x_1(1))^2 + (A(4).y - x_1(2))^2);         % distance to A4
    r_o = r_o';                                                        % transpose r_o
    
    H(1,1) = -(A(1).x - x_1(1)) / r_o(1);                    % H2 matrix
    H(1,2) = -(A(1).y - x_1(2)) / r_o(1);                    % H2 matrix

    H(2,1) = -(A(2).x - x_1(1)) / r_o(2);                    % H2 matrix
    H(2,2) = -(A(2).y - x_1(2)) / r_o(2);                    % H2 matrix

    H(3,1) = -(A(3).x - x_1(1)) / r_o(3);                    % H2 matrix
    H(3,2) = -(A(3).y - x_1(2)) / r_o(3);                    % H2 matrix

    H(4,1) = -(A(4).x - x_1(1)) / r_o(4);                    % H2 matrix
    H(4,2) = -(A(4).y - x_1(2)) / r_o(4);                    % H2 matrix

    R = eye(4);                                              % R2 - Covariance of measurements 

    dz = z_meas(i:i+3, 1) - r_o;                             % misclosure vector  

    N = R + H*P_1*H';                                        % N    
    N_ = mldivide(N, eye(size(N)));                          % N inverse

    k = P_1 * H' *( N_ );                                    % Gain matrix

    dx = dx_1 + k * ( dz  - H * dx_1);                       % Corrections to the state vector

    P = (eye(2) - k * H) * P_1;                              % P2
    
    i = i + 4;                                               % Move to next epoch
    
    P_1 = P;                                                 % set next iteration values 
    dx_1 = dx;
    
    x(:, t) = x(:, t-1) + dx;
    
    t = t + 1;
    
end

x_final = x_o + dx;

% final residuals
r_o = 0;
r_o(1) = sqrt( (A(1).x - x_final(1))^2 + (A(1).y - x_final(2))^2);         % distance to A1
r_o(2) = sqrt( (A(2).x - x_final(1))^2 + (A(2).y - x_final(2))^2);         % distance to A2
r_o(3) = sqrt( (A(3).x - x_final(1))^2 + (A(3).y - x_final(2))^2);         % distance to A3
r_o(4) = sqrt( (A(4).x - x_final(1))^2 + (A(4).y - x_final(2))^2);         % distance to A4
r_o = repmat(r_o, 1, (n / 4))';                                            % match size of z_meas

r = r_o - z_meas;                                                          % residuals

aposteriori = r' * eye(n) * r / (n - m);                       % A posteriori variance factor

H(1:4:n, 1) = -(A(1).x - x_final(1)) ./ r_o(1:4:n);                    % H matrix
H(1:4:n, 2) = -(A(1).y - x_final(2)) ./ r_o(1:4:n);                    % H matrix

H(2:4:n, 1) = -(A(2).x - x_final(1)) ./ r_o(2:4:n);                    % H matrix
H(2:4:n, 2) = -(A(2).y - x_final(2)) ./ r_o(2:4:n);                    % H matrix

H(3:4:n, 1) = -(A(3).x - x_final(1)) ./ r_o(3:4:n);                    % H matrix
H(3:4:n, 2) = -(A(3).y - x_final(2)) ./ r_o(3:4:n);                    % H matrix

H(4:4:n, 1) = -(A(4).x - x_final(1)) ./ r_o(4:4:n);                    % H matrix
H(4:4:n, 2) = -(A(4).y - x_final(2)) ./ r_o(4:4:n);                    % H matrix


N = H'*eye(n)*H;                                               % N
N_ = mldivide(N, eye(size(N)));                                % N inverse

P = aposteriori * N_;


% MAKE SURE TO RETURN x_final or x based on which task you are using!!!!!

end
