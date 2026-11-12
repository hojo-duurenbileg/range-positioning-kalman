function [x, P, aposteriori, r] = SummationLS(x_o, z_meas, A)

%% x

initial_poe = x_o;                                                     % initial point of expansion

n = size(z_meas, 1) * 4;                                               % number of rows 
m = 2;                                                                 % number of columns

z_meas = reshape(z_meas', [1, n])';                                    % vectorize measurements

H = zeros(n, 2);                                                       % allocate H matrix

dx = [1 1]';                                                           % set high value to begin iteration
iter = 0;                                                              % count iterations


while (abs(dx(1)) > 0.0001) && (abs(dx(2)) > 0.0001)
    
    iter = iter + 1;                                                   % count iterations

    r_o = 0;                                                           % set r_o zero;    
    
    r_o(1) = sqrt( (A(1).x - x_o(1))^2 + (A(1).y - x_o(2))^2);         % distance to A1
    r_o(2) = sqrt( (A(2).x - x_o(1))^2 + (A(2).y - x_o(2))^2);         % distance to A2
    r_o(3) = sqrt( (A(3).x - x_o(1))^2 + (A(3).y - x_o(2))^2);         % distance to A3
    r_o(4) = sqrt( (A(4).x - x_o(1))^2 + (A(4).y - x_o(2))^2);         % distance to A4
    
    r_o = repmat(r_o, 1, (n / 4))';                                    % match size of r_o
    
    dz = z_meas - r_o;                                                 % misclosure vector  

   	H(1:4:n, 1) = -(A(1).x - x_o(1)) ./ r_o(1:4:n);                    % H matrix
    H(1:4:n, 2) = -(A(1).y - x_o(2)) ./ r_o(1:4:n);                    % H matrix
    
    H(2:4:n, 1) = -(A(2).x - x_o(1)) ./ r_o(2:4:n);                    % H matrix
    H(2:4:n, 2) = -(A(2).y - x_o(2)) ./ r_o(2:4:n);                    % H matrix

    H(3:4:n, 1) = -(A(3).x - x_o(1)) ./ r_o(3:4:n);                    % H matrix
    H(3:4:n, 2) = -(A(3).y - x_o(2)) ./ r_o(3:4:n);                    % H matrix

    H(4:4:n, 1) = -(A(4).x - x_o(1)) ./ r_o(4:4:n);                    % H matrix
    H(4:4:n, 2) = -(A(4).y - x_o(2)) ./ r_o(4:4:n);                    % H matrix

    R = eye(size(H, 1));                                               % Covariance of measurements 
    R_ = mldivide(R, eye(size(R)));                                    % Inverse R matrix 
    
    N = H(1:4,:)' * R_(1:4,1:4) * H(1:4,:);                            % First sum of N 
    b = H(1:4,:)' * R_(1:4,1:4) * dz(1:4);                             % First sum of B 
    
    % Sum all N and b
    for i = 5:4:n  
        N = N + H(i:(i+3),:)'*R_(i:(i+3),i:(i+3))*H(i:(i+3),:);        
        b = b + H(i:(i+3),:)'*R_(i:(i+3),i:(i+3))*dz(i:(i+3));        
    end
        
    N_ = mldivide(N, eye(size(N)));
    
    dx = N_ * b;

    x_o = x_o + dx;   
    
end

initial_r_o(1) = sqrt( (A(1).x - initial_poe(1))^2 + (A(1).y - initial_poe(2))^2);
initial_r_o(2) = sqrt( (A(2).x - initial_poe(1))^2 + (A(2).y - initial_poe(2))^2);
initial_r_o(3) = sqrt( (A(3).x - initial_poe(1))^2 + (A(3).y - initial_poe(2))^2);
initial_r_o(4) = sqrt( (A(4).x - initial_poe(1))^2 + (A(4).y - initial_poe(2))^2);

initial_r_o = repmat(initial_r_o, 1, (n / 4))';

initial_dz = z_meas - initial_r_o;

r = initial_dz;

x = x_o;

%% P

aposteriori = r' * R * r / (n - m);                                    % A posteriori variance factor

P = aposteriori * N_;                                                  % Covariance of state


end

