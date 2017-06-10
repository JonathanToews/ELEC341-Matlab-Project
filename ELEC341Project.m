%te4sting

% Define Matrices for System
A = [3 0 1 0; 0 0 0 1; -9/10 9/10 -9/100 9/100; 1/10 -1/10 1/100 -1/100]
B = [0 0; 9/10; 0]
C = [0 1 0 0]
D = [0]

% Question #2: Compute Open-Loop Eigenvalues
eigs = eig(A)

% Question #3: Not sure if for matlab
% This is what I'm thinking
% For the responses to unit step changes in the input I
    % We can get the response (Y(s)) using 
    % Y(s) = (C(sI-A)^-1 B + D)U(s) (assuming x(0) = 0)
    % We then find (sI-A)^-1 and calculate out Y(s) by setting Td = 0 (also
    % set I = 1/s because unit step)
    % Then maybe we plot this?
% For the responses to unit step changes in the input Td
    % We use the same formula as above, but this time we set I = 0 (also Td
    % = 1/s because unit step)