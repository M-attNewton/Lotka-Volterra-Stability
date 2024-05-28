
%% Set up variables
clear all
syms x1 x2 x3
vars = [x1; x2; x3];

% Define parameters, values are not D-stable
A = [   -2.6755    1.3400   -0.2294;
   -0.1028   -0.0963    0.2026;
     0.2621   -1.9300   -1.6367];
A = [    -0.3792    0.0457    0.4302;
    -1.1050   -0.6037    0.3413;
    -0.1886    0.0455   -0.0847];

X_star = [1;0.5;0.25];
B = -A*X_star;

% Analytic solution for equilbrium
%sequil = [(c*(a-b))/(a*c-b*d), (a*(c-d))/(a*c-b*d)];
sequil = X_star;

% % Shifted vector field
X = [x1 + sequil(1); x2 + sequil(2); x3 + sequil(3)];
f = diag(X)*B + diag(X)*(A*X);
%x1_dot = a*(x1+(sequil(1))) - a*(x1+(sequil(1)))^2 - b*(x1+(sequil(1)))*(x2+(sequil(2)));
%x2_dot = c*(x2+(sequil(2))) - c*(x2+(sequil(2)))^2 - d*(x1+(sequil(1)))*(x2+(sequil(2)));
% f2 = [x1_dot; x2_dot];

% Initialize the sum of squares program
prog = sosprogram(vars);

% The Lyapunov function V(x): 

% Log-linear coefficients
[prog,coeff_1] = sospolyvar(prog,sym(1));
[prog,coeff_2] = sospolyvar(prog,sym(1));
[prog,coeff_3] = sospolyvar(prog,sym(1));
[prog,coeff_4] = sospolyvar(prog,sym(1));
[prog,coeff_5] = sospolyvar(prog,sym(1));
[prog,coeff_6] = sospolyvar(prog,sym(1));
[prog,coeff_7] = sospolyvar(prog,sym(1));

% Quadratic terms
[prog,V_quad] = sossosvar(prog,vars);
%V_quad = 0;

V = V_quad + coeff_1*x1 + coeff_2*x2 + coeff_3*x3 + coeff_4*log(x1 + sequil(1)) + ...
    coeff_5*log(x2 + sequil(2)) + coeff_6*log(x3 + sequil(3)) + coeff_7;

% Set up polynomials
Z1 = monomials(vars,0:1); 
Z2 = monomials(vars,0:2);

%[prog,V] = sospolyvar(prog,Z1,'wscoeff')
[prog,p1] = sospolyvar(prog,Z2,'wscoeff');
[prog,p2] = sospolyvar(prog,Z2,'wscoeff');
[prog,p3] = sospolyvar(prog,Z2,'wscoeff');

%% Define SOSP constraints
% Constraint 1: -dV/dx*f >= 0
% Region D
D1 = -x1 - sequil(1); D2 = -x2 - sequil(2);  D3 = -x3 - sequil(3);
expr2 = -(diff(V,x1)*f(1) + diff(V,x2)*f(2) + diff(V,x3)*f(3)) + p1*D1 + p2*D2 + p3*D3;
prog = sosineq(prog,expr2);

% Constraint 2: p >= 0
prog = sosineq(prog,p1);
prog = sosineq(prog,p2);
prog = sosineq(prog,p3);

% Constraint 3: (d^2V/dx^2 >= 0)
% Q = subs(jacobian(jacobian(V,[x1,x2]),[x1,x2]),[x1,x2],[sequil(1),sequil(2)]);
% prog = sosineq(prog,[x1,x2]*Q*[x1;x2] - 0.1*x1^2 - 0.1*x2^2);
Q = subs(jacobian(jacobian(V,[x1,x2,x3]),[x1,x2,x3]),[x1,x2,x3],[0,0,0]);
prog = sosineq(prog,[x1,x2,x3]*Q*[x1;x2;x3] - 0.1*x1^2 - 0.1*x2^2 - 0.1*x3^2);
% prog = sosineq(prog,-coeff_3);
% prog = sosineq(prog,-coeff_4);

% Constraint 4: Fix parameters to match true Lyapunov function
% expr3 = coeff_5 + coeff_3*log(sequil(1)) + coeff_4*log(sequil(2));
% prog = soseq(prog,expr3);
%
prog = soseq(prog,coeff_1 - 1);
prog = soseq(prog,coeff_2 - 1);
prog = soseq(prog,coeff_3 - 1);

%% Call solver
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);

%% Get solution
SOLV = sosgetsol(prog,V)

% True Lyapunov function
K = sequil(1)*log(sequil(1)) + sequil(2)*log(sequil(2));
trueLyap = x1 - sequil(1)*log(x1 + sequil(1)) + x2 - sequil(2)*log(x2 + sequil(2)) + K;

% Plot solution and true function
axis1 = -1; axis2 = 3;
subplot(2,2,1)
fsurf(trueLyap,[axis1 axis2])
title('True')

subplot(2,2,2)
fsurf(SOLV,[axis1 axis2])
title('Solution')

subplot(2,2,3)
fcontour(trueLyap,[axis1 axis2])
title('True')

subplot(2,2,4)
fcontour(SOLV,[axis1 axis2])
title('Solution')


