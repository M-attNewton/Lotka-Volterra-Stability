
%% Set up variables
clear all
syms x1 x2 x3
vars = [x1; x2; x3];

% Define parameters, values are not D-stable
A = [-1.0460    1.2383    0.0411;
    -2.6477    0.7553    1.4488;
    -2.9632    0.3359   -0.2762];

X_star = [1;0.5;0.25];
B = -A*X_star;

% Analytic solution for equilbrium
%sequil = [(c*(a-b))/(a*c-b*d), (a*(c-d))/(a*c-b*d)];
sequil = X_star;

% % Shifted vector field
X = [x1 + sequil(1); x2 + sequil(2); x3 + sequil(3)];
f = diag(X)*B + diag(X)*(A*X);
f2 = f - 10*[(x1 + sequil(1))*x2*x3; x1*(x2 + sequil(2))*x3; -2*x1*x2*(x3 + sequil(3))];

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
%vars2 = [x1 - sequil(1); x2 - sequil(2); x3 - sequil(3)];
[prog,V_quad] = sossosvar(prog,vars);
V_quad = subs(V_quad, [x1;x2;x3], [x1 - sequil(1); x2 - sequil(2); x3 - sequil(3)]);
%V_quad = 0;

V = coeff_1*x1 + coeff_2*x2 + coeff_3*x3 + coeff_4*log(x1 + sequil(1)) + ...
    coeff_5*log(x2 + sequil(2)) + coeff_6*log(x3 + sequil(3)) + coeff_7 + V_quad;

%V = V_quad + coeff_4*log(x1) + coeff_5*log(x2) + coeff_6*log(x3);

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
% expr2 = -(diff(V_quad,x1)*f2(1) + diff(V_quad,x2)*f2(2) + diff(V_quad,x3)*f2(3))...
%     -(diff(V,x1)*f(1) + diff(V,x2)*f(2) + diff(V,x3)*f(3)) ...
%     + (coeff_1 + coeff_2 + coeff_3)*x1*x2*x3 ...
%    ... %+ (coeff_4/(x1 + sequil(1)))...
%     + p1*D1 + p2*D2 + p3*D3;...
%     
% % Multiply by (x1 + sequil(1))*(x2 + sequil(2))*(x3 + sequil(3)) to deal
% % with SOS errors
% expr2 = -(x1 + sequil(1))*(x2 + sequil(2))*(x3 + sequil(3))* ...
%     ( diff(V_quad,x1)*f2(1) + diff(V_quad,x2)*f2(2) + diff(V_quad,x3)*f2(3) ...
%     + diff(V,x1)*f(1) + diff(V,x2)*f(2) + diff(V,x3)*f(3) ...
%     - (coeff_1 + coeff_2 + coeff_3)*x1*x2*x3) ...
%     - x1*x2*x3*(coeff_4*(x2 + sequil(2))*(x3 + sequil(3))...
%     + coeff_5*(x1 + sequil(1))*(x3 + sequil(3))...
%     + coeff_6*(x2 + sequil(2))*(x1 + sequil(1)))...
%     + p1*D1 + p2*D2 + p3*D3;

expr2 = -(diff(V,x1)*f2(1) + diff(V,x2)*f2(2) + diff(V,x3)*f2(3))...
        + p1*D1 + p2*D2 + p3*D3;
    %-(diff(V,x1)*x1*x2*x3 + diff(V,x2)*x1*x2*x3 + diff(V,x3)*x1*x2*x3);
prog = sosineq(prog,expr2);

% Constraint 2: p >= 0
prog = sosineq(prog,p1);
prog = sosineq(prog,p2);
prog = sosineq(prog,p3);

% Constraint 3: (d^2V/dx^2 >= 0)
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
solver_opt.solver = 'sdpt3'; %sedumi
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


