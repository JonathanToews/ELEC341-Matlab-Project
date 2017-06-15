% Question #1
% Define constants
J1=10/9; J2=10; c=0.1; k=1; kI=1;

% Define State Space
A = [0 0 1 0; 0 0 0 1; -k/J1 k/J1 -c/J1 c/J1; k/J2 -k/J2 c/J2 -c/J2];
B = [0; 0; kI/J1; 0];
C = [0 1 0 0];
D = [0];
F = [0; 0; 0; 1/J2];         % disturbance torque vector


% Question #2: Compute Open-Loop Eigenvalues
eigs = eig(A)


% Question #3: Simulate the unit step changes in I and Td
%figure(1);


sys1 = ss(A, B, C, D);      % define sys1 as the ss with input I (Td = 0)
sys2 = ss(A, F, C, D);      % define sys2 as the ss with input Td (I = 0)
%step(sys1, sys2, 25000)     % Apply unit step to both models and plot
%title('#3 Unit step changes in the input I and disturbance torque Td');
%legend('Input response', 'Disturbance response');

% Controllability of the open-loop system
Sc = ctrb(A, B);
assert(length(A) == rank(Sc));

% Question #4: Design state feedback controller gain K
CCLP = [-1 -2 -1-1i -1+i];       % closed loop poles
K = acker(A, B, CCLP)            % controller gain K vector

% Question #5: Compute Kr
Kr = (-1 / (C*inv(A-B*K)*B))    % ensure that Kr matrix is of type double

% Question #6: Simulate closed loop response to unit step changes
%figure();

A_fb = A-B*K;
B_fb = Kr.*B;
C_fb = C;
D_fb = D;

A_fb_td = A-B*K;
B_fb_td = [Kr.*B F];
C_fb_td = C;
D_fb_td = [D 0];

% Find settling time for setpoint
sys_fb = ss(A_fb, B_fb, C_fb, D_fb);
step_info = stepinfo(sys_fb);
settling_time = step_info.SettlingTime

% Ensure that the graph is large enough to see settling time and then some
Tstop = round(10 * settling_time)
delT = 0.01;
T = 0:delT:Tstop;

% Step function for setpoint from time = 0
step_setpoint = ones(length(T), 1);

% Step function for Td after steady-state. Start by ensuring we are in the
% steady state before incorporating Td
step_td_start_index = 5 * round(settling_time / delT)

% Vector of 0s until well after settling time (then becomes 1 to
% incorporate Td
step_td = [zeros(step_td_start_index, 1) ; ones(length(T) - step_td_start_index, 1)];

% 2 columned matrix for the input (1st column is setpoint input; 2nd column
% is Td input)
input = [step_setpoint step_td];

sys_cl = ss(A_fb_td, B_fb_td, C_fb_td, D_fb_td);    % construct a system model
[Y, Tsim, X] = lsim(sys3, input, T);                % simulate
plot(Tsim,Y)                                        % plot the output vs. time
title('Step Response with Zero Initial Conditions')

% sys_cl = ss(A_fb, B_fb, C_fb, D_fb);
% Ts = 0.1;
% sys_cl_d = c2d(sys_cl, Ts);
% sysd = setmpcsignals(sys_cl_d, 'MV', 1, 'MO', 1);
% MPCobj = mpc(sysd);
% Tstop = 30;
% 
% num_sim_steps = round(Tstop/Ts);
% r = ones(num_sim_steps, 1);
% v = [zeros(2*num_sim_steps/3, 1); ones(num_sim_steps/3, 1)];
% 
% sim(MPCobj, num_sim_steps, r)

% sys_cl = ss(A_fb_td, B_fb_td, C_fb_td, D_fb_td);
% Ts = 0.1;
% sys_cl_d = c2d(sys_cl, Ts);
% sysd = setmpcsignals(sys_cl_d, 'MD', 1, 'MD', 2, 'MO', 1);
% MPCobj = mpc(sysd);
% Tstop = 30;
% 
% num_sim_steps = round(Tstop/Ts);
% r = ones(num_sim_steps, 1);
% v = [zeros(2*num_sim_steps/3, 1); ones(num_sim_steps/3, 1)];
% 
% sim(MPCobj, num_sim_steps, r, v)

% Question #7:
OCLP = [-4 -2 -2-2i -2+2i];
L = acker(A', C', OCLP)';

% Question #8:
figure();
sys4 = ss(A-L*C, B, C, D);
step(sys4);

[y,t,x] = step(sys4, 8);
plot(t,x)
hold on
legend('x0', 'x1', 'x2', 'x3')
title('#8 Open Loop Observer Simulation')

% Question #9:

ACL = [A-B*K B*K; zeros(length(A)) A-L*C];
BCL = [B*Kr; zeros(length(A), 1)];
CCL = [C zeros(1, length(A))];
figure();
sys6 = ss(ACL, BCL, CCL, D)
step(sys6, 8);
ylim([0 1.5])
title('#9 Combined Observer Feedback Simulation')


% Appendix

% Sanity Check 

%syms K1 K2 K3 K4 s;                 % declare controller variables
%K = [K1, K2, K3, K4];               % gain vector K

% We want closed-loop poles of -1, -2, -1+-j. This polynomial provides
% those poles
%PEQ1 = (s + 2)*(s + 1)*(s^2 + 2*s + 2);

% This is the polynomial provided by the closed loop polynomial with feedback
% det(sI - (A - BK))
%PEQ2 = det(s.*eye(4) - (A - B*K));

% In order to find the values of K that will provide the closed-loop poles
% we want given our state-space model, we must equate the two above
% equations and solve for K
%TEMP = coeffs(PEQ1 - PEQ2, s);
%[K1, K2, K3, K4] = solve(TEMP(1)==0, TEMP(2)==0, TEMP(3)==0, TEMP(4)==0);
%K = double([K1, K2, K3, K4])        % ensure that K matrix is of type double


%syms L1 L2 L3 L4;
%L = [L1; L2; L3; L4];

%PEQ1 = (s+4)*(s+2)*(s^2 + 4*s + 8); % desired closed loop polynomial
%PEQ2 = det(s.*eye(4) - (A - L*C));    % closed loop polynomial with feedback

%TEMP = coeffs(PEQ1-PEQ2, s);        % coefficients to solve observer vals
%[L1, L2, L3, L4] = solve(TEMP(1)==0, TEMP(2)==0, TEMP(3)==0, TEMP(4)==0);
%L = [L1; L2; L3; L4]                % observer gain vector L