% Define constants
J1=10/9; J2=10; c=0.1; k=1; kI=1;
x_plant = [0 0 0 0];                     % Initial States
x_obs = [0 0 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question #1
% Determine state space representation
A = [0 0 1 0; 0 0 0 1; -k/J1 k/J1 -c/J1 c/J1; k/J2 -k/J2 c/J2 -c/J2];
BI = [0; 0; kI/J1; 0];
BTd = [0; 0; 0; 1/J2];         % disturbance torque vector
C = [0 1 0 0];
D = [0];

Sc = ctrb(A, BI);
assert(length(A) == rank(Sc));

% Find gain vector K given the closed loop poles below
CCLP = [-1 -2 -1-1i -1+i];
K = place(A, BI, CCLP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question #5: Compute Kr
Kr = (-1 / (C*inv(A-BI*K)*BI))



So = obsv(A, C)
assert(length(A) == rank(So));

% Find gain vector L given the closed loop poles below
OCLP = [-4 -2 -2-2i -2+2i];
L = place(A', C', OCLP)'

% get_param('ClosedLoop/A', 'Gain')
% get_param('ClosedLoop/BI', 'Gain')
% get_param('ClosedLoop/K', 'Gain')
% set_param('ClosedLoop/A', 'Gain', mat2str(A))
% set_param('ClosedLoop/BI', 'Gain', mat2str(B_I))
% set_param('ClosedLoop/C', 'Gain', mat2str(C))
% set_param('ClosedLoop/K', 'Gain', mat2str(K))
% set_param('ClosedLoop/Kr', 'Gain', num2str(Kr))

simOut = sim('ClosedLoop');
% simOut = sim('OpenLoop_1')