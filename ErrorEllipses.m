function [ a, b, beta ] = ErrorEllipses( P )
%ERRORELLIPSES Summary of this function goes here
%   Detailed explanation goes here


v_X = P(1,1);
v_Y = P(2,2);

v_XY = P(1,2);

a = sqrt( 0.5 * ( v_X + v_Y + sqrt( (v_X - v_Y)^2 + (4 * v_XY^2) ) ) );

b = sqrt( 0.5 * ( v_X + v_Y - sqrt( (v_X - v_Y)^2 + (4 * v_XY^2) ) ) );

beta = (pi/2) - atan2(2 * v_XY, (v_X - v_Y)) / 2; 

end

