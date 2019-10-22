function [x, P] = mu_mag(x, P, ymag, Rm, m0)
% Measurement update when accelerometer measurements are available

% Sanity checks
if (any(isnan(x)) | any(isnan(P)) | any(isnan(ymag)))
    x
    P
    ymag
    %keyboard
end

q = x(1:4);
R = Qq(q);

% The innovation
e = ymag - R'*m0; 

% The matrix H in the linearized model y = H x + e
[Q1, Q2, Q3, Q4] = dQqdq(q);
H = cat(2, Q1'*m0, Q2'*m0, Q3'*m0, Q4'*m0);

% The innovation covariance
S = H*P*H' + Rm;

% The Kalman gain
K = (S\(H*P'))'; % or
% K = P*H'*inv(S);

% The estimate update
x = x + K*e;

% The covariance matrix update
P = (eye(4) - K*H)*P;
