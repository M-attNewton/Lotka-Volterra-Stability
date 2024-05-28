%% Implementing the Whirling Pendulum example in IEEE CDC 2002 by AP

%% Set up variables
clear all
syms x1 x2 u1 u2 %u1 = ln(x1), u2 = ln(x2)
vars = [x1; x2; u1; u2];

% Define parameters
X_star = [1;0.5];

A =[ -0.7027    0.6120
     -0.8975    0.5415];
 
b = [0.3967;
     0.6267];

% Modified vector field
x1_dot = b(1)*(x1+(2/3)) + A(1,1)*(x1+(2/3))^2 + A(1,2)*(x1+(2/3))*(x2+(1/2));
x2_dot = b(2)*(x2+(1/2)) + A(2,1)*(x2+(1/2))^2 + A(2,2)*(x1+(2/3))*(x2+(1/2));
f = [x1_dot;
     x2_dot;
     ((x1+(2/3))^-1)*x1_dot;
     ((x2+(1/2))^-1)*x2_dot];
 
% Initialize the sum of squares program
prog = sosprogram(vars);

% The Lyapunov function V(x): 
Z1 = monomials(vars,0:1); %Z1 = [x1; x2; u1; u2; 1];
Z2 = monomials(vars,0:2);
[prog,V] = sospolyvar(prog,Z1,'wscoeff')
%[prog,p1] = sospolyvar(prog,Z2,'wscoeff')
%[prog,p2] = sospolyvar(prog,Z2,'wscoeff')

%% Define SOSP constraints
% Constraint 1: -dV/dx*f >= 0
% Region D
%D1 = -x1 - (2/3); D2 = -x2 - (1/2);
expr2 = -(diff(V,x1)*f(1) + diff(V,x2)*f(2) + diff(V,u1)*f(3) + diff(V,u2)*f(4));% + ... 
     %p1*D1 + p2*D2;
prog = sosineq(prog,expr2);

% Constraint 2: p >= 0
%prog = sosineq(prog,p1);
%prog = sosineq(prog,p2);

% Constraint 3: coeff_4 <= 0, coeff_5 <= 0
prog = sosineq(prog,-coeff_4);
prog = sosineq(prog,-coeff_5);

% Constraint 4: Fix parameters to match true Lyapunov function
%expr3 = coeff_1 + coeff_2*(2/3) + coeff_3*(1/2) + coeff_4*log(2/3) + coeff_5*log(1/2);
%expr3 = coeff_1 + coeff_2*(0) + coeff_3*(0) + coeff_4*log(0) + coeff_5*log(0);
%prog = soseq(prog,expr3);
prog = soseq(prog,coeff_2 - 1);
prog = soseq(prog,coeff_3 - 1);

%% Call solver
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);

%% Get solution
SOLV = sosgetsol(prog,V)

% Subsitute in log terms and reduce significant figures
SOLV2 = subs(SOLV, u1, log(x1));
SOLV3 = subs(SOLV2, u2, log(x2))
SOLV3 = vpa(SOLV3,4)

% True Lyapunov function
K = -1.7836;
%K = 0;
trueLyap = x1 - (2/3)*log(x1) + x2 - (1/2)*log(x2) + K;
%trueLyap = x1+(2/3) - (2/3)*log(x1+(2/3)) + x2+(1/2) - (1/2)*log(x2+(1/2)) + K;

% Plot solution and true function
axis1 = -1; axis2 = 3;
subplot(2,2,1)
fsurf(trueLyap,[axis1 axis2])
title('True')

subplot(2,2,2)
fsurf(SOLV3,[axis1 axis2])
title('Solution')

subplot(2,2,3)
fcontour(trueLyap,[axis1 axis2])
title('True')

subplot(2,2,4)
fcontour(SOLV3,[axis1 axis2])
title('Solution')
