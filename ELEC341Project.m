% Question #1

% Define constants
J1=10/9; J2=10; c=0.1; k=1; kI=1;

% Define State Space
A = [0 0 1 0; 0 0 0 1; -k/J1 k/J1 -c/J1 c/J1; k/J2 -k/J2 c/J2 -c/J2]
B = [0; 0; kI/J1; 0]
C = [0 1 0 0]
D = [0]
F = [0; 0; 0; 1/J2]         % disturbance torque vector


% Question #2: Compute Open-Loop Eigenvalues
eigs = eig(A)


% Question #3: Simulate the unit step changes input I and in
% disturbance torque Td
sys1 = ss(A, B, C, D);      % define sys1 as the ss with input I (Td = 0)
sys2 = ss(A, F, C, D);      % define sys2 as the ss with input Td (I = 0)

% Apply unit step to both models and plot
step(sys1, sys2, 25000)
title('Unit step changes in the input I and disturbance torque Td');

% Controllability of the open-loop system
cont = [B A*B A^2*B A^3*B];
det_cont = det(cont)


% Question #4:
syms K1 K2 K3 K4 s;                 % declare controller variables
K = [K1, K2, K3, K4];               % gain vector K

% We want closed-loop poles of -1, -2, -1+-j. This polynomial provides
% those poles
PEQ1 = (s + 2)*(s + 1)*(s^2 + 2*s + 2);

% This is the polynomial provided by the closed loop polynomial with feedback
% det(sI - (A - BK))
PEQ2 = det(s.*eye(4) - (A - B*K));

% In order to find the values of K that will provide the closed-loop poles
% we want given our state-space model, we must equate the two above
% equations and solve for K
TEMP = coeffs(PEQ1 - PEQ2, s);
[K1, K2, K3, K4] = solve(TEMP(1)==0, TEMP(2)==0, TEMP(3)==0, TEMP(4)==0);
K = double([K1, K2, K3, K4])        % ensure that K matrix is of type double


% Question #5:
Kr = double(-1 / (C*inv(A-B*K)*B))  % ensure that Kr matrix is of type double


% Question #6:
sys3 = ss(A-B*K, Kr.*B, C, D);
step(sys3, 20);
ylim([0 1.5])

% Question #7:
syms L1 L2 L3 L4;
L = [L1; L2; L3; L4];

PEQ1 = (s+4)*(s+2)*(s^2 + 4*s + 8); % desired closed loop polynomial
PEQ2 = det(s.*eye(4) - (A - L*C));    % closed loop polynomial with feedback

TEMP = coeffs(PEQ1-PEQ2, s);        % coefficients to solve observer vals
[L1, L2, L3, L4] = solve(TEMP(1)==0, TEMP(2)==0, TEMP(3)==0, TEMP(4)==0);
L = [L1; L2; L3; L4]                % observer gain vector L
