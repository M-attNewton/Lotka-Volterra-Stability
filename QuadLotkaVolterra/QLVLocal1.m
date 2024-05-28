QuqQQuad%% Predator prey dynamics with quadratic Lyapunov function

%% Set up variables
clear all; %echo on;
syms x1 x2;
vars = [x1; x2];

% Define parameters
a = 3; b = 2; c = 4; d = 3;

% Construct the vector field dx/dt = f
f = [a*x1 - a*x1^2 - b*x1*x2;
     c*x2 - c*x2^2 - d*x1*x2];

% Initialize the sum of squares program
prog = sosprogram(vars);

% The Lyapunov function V(x): 
Z1 = monomials(vars,0:2);
Z2 = monomials(vars,0:2);
[prog,V] = sospolyvar(prog,Z1,'wscoeff')
[prog,p] = sospolyvar(prog,Z2,'wscoeff')

%% Define SOSP constraints
% Constraint 1 : V(x) - (x1^2 + x2^2) >= 0
prog = sosineq(prog,V - (x1^2 + x2^2));

% Constraint 2: -dV/dx*f >= 0
expr2 = -(diff(V,x1)*f(1) + diff(V,x2)*f(2)) - p*(x1^2 + x2^2 - 2);
prog = sosineq(prog,expr2);


%% Call solver
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);

%% Get solution
SOLV = sosgetsol(prog,V)
%echo off;

% True Lyapunov function
K = 0;
%trueLyap = d*x1 - c*log(x1) + b*x2 - a*log(x2) + K
trueLyap = x1 - (2/3)*log(x1) + x2 - (1/2)*log(x2);


%% Junk


% Constraint 3: 0 <= x1 <= 1, 0 <= x2 <= 1, p(x) = x1^2 + x2^2
%expr3 = (x1^2 + x2^2)*(-x1 + x1 - 1 - x2 + x2 - 1);
%expr3 = x1^2*(-x1-1) + x1^4*(x1-1) + x2^2*(-x2-1) + x2^4*(x2-1);
%expr3 = (x1^2 + x2^2)*(-x1-1-x2-1) + (x1^4+x2^4)*(x1-1+x2-1);

%prog = soseq(prog, -x1^3); prog = soseq(prog, x1^2*(x1-1)); 
%prog = soseq(prog, -x2^3); prog = soseq(prog, x2^2*(x2-1));  


%prog = sosmatrixineq(prog, [x1-0.2;x2-0.2]);
%prog = sosmatrixineq(prog, [-x1+0.8;-x2+0.8]);

% prog = sosmatrixineq(prog, [-x1+0.2;-x2+0.2]);
% prog = sosmatrixineq(prog, [x1-0.8;x2-0.8]);