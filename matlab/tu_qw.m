function [x, P] = tu_qw(x, P, omega, T, Rw)

% Use nonlinear model to update x, linear model to update P

Sw = Somega(omega);
%Sqq = Sq(x);
wThalf = 0.5*norm(omega)*T;
x = (cos(wThalf)*eye(4) + T/2*(sin(wThalf)/wThalf)*Sw)*x;

F = eye(4) + 0.5*T*Sw;

P = F*PF + Rw;
