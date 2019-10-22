function [x, P] = mu_acc(x, P, yacc, Ra, g0)
% Measurement update when accelerometer measurements are available

q = x(1:4);
R = Qq(q);

% The innovation
e = yacc - R'*g0; 

% The matrix H in the linearized model y = H x + e
[Q1, Q2, Q3, Q4] = dQqdq(q);
H = cat(2, Q1'*g0, Q2'*g0, Q3'*g0, Q4'*g0);

% The innovation covariance
S = H*P*H' + Ra;

% The Kalman gain
K = (S\(H*P'))'; % or
% K = P*H'*inv(S);

% The estimate update
x = x + K*e;

% The covariance matrix update
P = (eye(4) - K*H)*P;
