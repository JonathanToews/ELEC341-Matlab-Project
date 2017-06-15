% Define constants
J1=10/9; J2=10; c=0.1; k=1; kI=1;
x0 = [1 2 3 4];                     % Initial States

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question #1
% Determine state space representation
A = [0 0 1 0; 0 0 0 1; -k/J1 k/J1 -c/J1 c/J1; k/J2 -k/J2 c/J2 -c/J2];
B = [0; 0; kI/J1; 0];
C = [0 1 0 0];
D = [0];
F = [0; 0; 0; 1/J2];         % disturbance torque vector

A_2 = A;
B_2 = [B F];
C_2 = C;
D_2 = [D 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question #2
% Compute Open-Loop Eigenvalues

eigs = eig(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question #3
% Simulate the unit step changes in I and Td

t = [0:0.01:30];
input = ones(length(t), 1);

% define sys_fb_I as the ss with input I (Td = 0)
sys_I = ss(A, B, C, D);
response_I = lsim(sys_I, input, t, x0);

figure(1);
subplot(311)
plot(t, response_I)
xlabel('Time (s)')
ylabel('\phi_{2} (rads)')
title('Output Response to Unit Step in I (Open-Loop)')

subplot(312)
plot(t, input)
xlabel('Time (s)')
ylabel('Input I')
title('Input Over Time')
ylim([0 2])

subplot(313)
plot(t, zeros(length(t), 1))
xlabel('Time (s)')
ylabel('Disturbance T_d')
title('Disturbance Over Time')
ylim([0 2])

% define sys_fb_td as the ss with input Td (I = 0)
sys_td = ss(A, F, C, D);
response_td = lsim(sys_td, input, t, x0);

figure(2);
subplot(311)
plot(t, response_td)
xlabel('Time (s)')
ylabel('\phi_{2} (rads)')
title('Output Response to Unit Step in T_d')

subplot(312)
plot(t, zeros(length(t), 1))
xlabel('Time (s)')
ylabel('Input I')
title('Input Over Time')
ylim([0 2])

subplot(313)
plot(t, input)
xlabel('Time (s)')
ylabel('Disturbance T_d')
title('Disturbance Over Time')
ylim([0 2])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question #4
% Design state feedback controller gain K

% Check controllability of the open-loop system (can only arbitrarily
% choose eigenvalues of A-BK if system is controllable)
Sc = ctrb(A, B);
assert(length(A) == rank(Sc));

% Find gain vector K given the closed loop poles below
CCLP = [-1 -2 -1-1i -1+i];
K = acker(A, B, CCLP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question #5: Compute Kr
Kr = (-1 / (C*inv(A-B*K)*B))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question #6
% Simulate closed loop response to unit step changes

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
t_stop = round(10 * settling_time)
delta_t = 0.01;
t = 0:delta_t:t_stop;

% Step function for setpoint from time = 0
step_setpoint = ones(length(t), 1);

% Step function for Td after steady-state. Start by ensuring we are in the
% steady state before incorporating Td
step_td_start_index = 5 * round(settling_time / delta_t)

% Vector of 0s until well after settling time (then becomes 1 to
% incorporate Td
step_td = [zeros(step_td_start_index, 1) ; ones(length(t) - step_td_start_index, 1)];

% 2 columned matrix for the input (1st column is setpoint input; 2nd column
% is Td input)
input = [step_setpoint step_td];

sys_cl = ss(A_fb_td, B_fb_td, C_fb_td, D_fb_td);    % construct a system model
response_cl = lsim(sys_cl, input, t, x0);

figure(3);
subplot(311)
plot(t, response_cl)
xlabel('Time (s)')
ylabel('\phi_{2} (rads')
title('Output Response to Unit Step in I and T_d (With Feedback)')

subplot(312)
plot(t, input(:,1))
xlabel('Time (s)')
ylabel('Input I')
title('Input Over Time')
ylim([0 2])

subplot(313)
plot(t, input(:,2))
xlabel('Time (s)')
ylabel('Disturbance T_d')
title('Disturbance Over Time')
ylim([0 2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question #7
% Design observer

% Can only arbitrarily choose observer poles if system is observable
So = obsv(A, C)
assert(length(A) == rank(So));

% Find gain vector L given the closed loop poles below
OCLP = [-4 -2 -2-2i -2+2i];
L = acker(A', C', OCLP)'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question #8
% Simulate open-loop system

t = [0:0.01:30];
input = ones(length(t), 1);

% define sys_obs_I as the open-loop sustem with input I (Td = 0)
sys_obs_I = ss(A-L*C, B, C, D);      
response_obs_I = lsim(sys_obs_I, input, t, x0);

figure(4);
subplot(311)
plot(t, response_obs_I)
xlabel('Time (s)')
ylabel('\phi_{2} (rads')
title('Output Response to Unit Step in I (Observer)')

subplot(312)
plot(t, input)
xlabel('Time (s)')
ylabel('Input I')
title('Input Over Time')
ylim([0 2])

subplot(313)
plot(t, zeros(length(t), 1))
xlabel('Time (s)')
ylabel('Disturbance T_d')
title('Disturbance Over Time')
ylim([0 2])

% define sys_obs_I as the open-loop sustem with input Td (I = 0)
sys_obs_td = ss(A-L*C, F, C, D);      
response_obs_td = lsim(sys_obs_td, input, t, x0);

figure(5);
subplot(311)
plot(t, response_obs_td)
xlabel('Time (s)')
ylabel('\phi_{2} (rads')
title('Output Response to Unit Step in Td (Observer)')

subplot(312)
plot(t, input)
xlabel('Time (s)')
ylabel('Input I')
title('Input Over Time')
ylim([0 2])

subplot(313)
plot(t, zeros(length(t), 1))
xlabel('Time (s)')
ylabel('Disturbance T_d')
title('Disturbance Over Time')
ylim([0 2])

%%%% DON"T KNOW WHAT TO DO HERE
% Plot x vs x^
[y,t,x_hat] = lsim(sys_obs_I, zeros(length(t), 1), t, x0)
figure(6);
subplot(111);
plot(t, x_hat)


% 
% n = length(A);
% x = x_original_sim;
% 
% x1_label = '\phi_{1}';
% x2_label = '\phi_{2}';
% x3_label = '\omega_{1}';
% x4_label = '\omega_{2}';
% 
% figure(6);
% subplot(211);
% plot(t,x_hat)
% 
% subplot(212);
% plot(t, x_original_sim);
% 
% plot(t,x(:,1),'-g', t,x_hat(:,1),':g', t,x(:,2),'-b', t,x_hat(:,2),':b', t,x(:,3),'-r', t,x_hat(:,3),':r', t,x(:,4),'-o', t,x_hat(:,4),':o')
% xlabel('Time (s)')
% ylabel('Amplitude')
% legend(x1_label, strcat(x1_label,'^'), strcat(x2_label, x2_label,'^'), x3_label, strcat(x3_label,'^'), x4_label, strcat(x4_label,'^'))
% title('#8 Open Loop Observer Simulation')


% Question #9:

A_cl = [A-B*K B*K; zeros(size(A)) A-L*C];
B_cl = [B*Kr; zeros(size(B))];
C_cl = [C zeros(size(C))];
D_cl = [D];

A_cl_td = A_cl;
B_cl_td = [B_cl [F ; zeros(size(B))]];
C_cl_td = C_cl;
D_cl_td = D_cl;

% Find settling time for setpoint
sys_cl = ss(A_cl, B_cl, C_cl, D_cl);
step_info = stepinfo(sys_cl);
settling_time = step_info.SettlingTime

% Ensure that the graph is large enough to see settling time and then some
t_stop = round(10 * settling_time)
delta_t = 0.01;
t = 0:delta_t:t_stop;

% Step function for setpoint from time = 0
step_setpoint = ones(length(t), 1);

% Step function for Td after steady-state. Start by ensuring we are in the
% steady state before incorporating Td
step_td_start_index = 5 * round(settling_time / delta_t)

% Vector of 0s until well after settling time (then becomes 1 to
% incorporate Td
step_td = [zeros(step_td_start_index, 1) ; ones(length(t) - step_td_start_index, 1)];

% 2 columned matrix for the input (1st column is setpoint input; 2nd column
% is Td input)
input = [step_setpoint step_td];

sys_cl_td = ss(A_cl_td, B_cl_td, C_cl_td, D_cl_td);    % construct a system model
response_cl_td = lsim(sys_cl_td, input, t, [x0 x0]);

figure(7);
subplot(311)
plot(t, response_cl_td)
xlabel('Time (s)')
ylabel('\phi_{2} (rads')
title('Output Response to Unit Step in I and T_d (With Feedback And Observer)')

subplot(312)
plot(t, input(:,1))
xlabel('Time (s)')
ylabel('Input I')
title('Input Over Time')
ylim([0 2])

subplot(313)
plot(t, input(:,2))
xlabel('Time (s)')
ylabel('Disturbance T_d')
title('Disturbance Over Time')
ylim([0 2])

% figure(8);
% subplot(311)
% sys_combined = ss(A_CL, B_CL, C_CL, D_CL)
% step(sys_combined, 8)
% ylim([0 1.5])
% title('#9 Combined Observer Feedback Simulation')


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




%For question 6
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