%% Implementing the Whirling Pendulum example in IEEE CDC 2002 by AP

%% Set up variables
clear all
syms x1 x2 x3
vars = [x1; x2; x3];

% Define parameters, values are D-stable
A = [    -1.1499   -0.5502   -1.0848;
    -1.0589   -1.7808    1.8048;
    -0.3987   -1.6770   -0.9149];

A = [-1.1863    0.6078    0.1033;
    -0.6940   -0.6834    0.2419;
    -0.1298   -1.6213   -1.0654];


A = [   -1.9513   -0.7378    1.2711;
   -0.1717   -1.0671    1.0499;
   -0.2224   -0.6060   -1.1032];

X_star = [1;0.5;0.25];
B = -A*X_star;   
sequil = X_star;

% Analytic solution for equilbrium
%sequil = [(c*(a-b))/(a*c-b*d), (a*(c-d))/(a*c-b*d)];
%sequil = A\B; %inv(A)*B;

% % Shifted vector field
X = [x1 + sequil(1); x2 + sequil(2); x3 + sequil(3)];
f = diag(X)*B + diag(X)*(A*X);

% Initialize the sum of squares program
prog = sosprogram(vars);

% The Lyapunov function V(x): 
[prog,coeff_1] = sospolyvar(prog,sym(1));
[prog,coeff_2] = sospolyvar(prog,sym(1));
[prog,coeff_3] = sospolyvar(prog,sym(1));
[prog,coeff_4] = sospolyvar(prog,sym(1));
[prog,coeff_5] = sospolyvar(prog,sym(1));
[prog,coeff_6] = sospolyvar(prog,sym(1));
[prog,coeff_7] = sospolyvar(prog,sym(1));

V = coeff_1*x1 + coeff_2*x2 + coeff_3*x3 + coeff_4*log(x1 + sequil(1)) + ...
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
D1 = -x1 - sequil(1); D2 = -x2 - sequil(2); D3 = -x3 - sequil(3); 
expr2 = -(diff(V,x1)*f(1) + diff(V,x2)*f(2) + diff(V,x3)*f(3)) + p1*D1 + p2*D2 + p3*D3;
prog = sosineq(prog,expr2);

% Constraint 2: p >= 0
prog = sosineq(prog,p1);
prog = sosineq(prog,p2);
prog = sosineq(prog,p3);

% Constraint 3: coeff_4 <= 0, coeff_5 <= 0 (d^2V/dx^2 <= 0)
prog = sosineq(prog,-coeff_4);
prog = sosineq(prog,-coeff_5);
prog = sosineq(prog,-coeff_6);

% Constraint 4: Fix parameters to match true Lyapunov function
%expr3 = coeff_5 + coeff_3*log(sequil(1)) + coeff_4*log(sequil(2));
%prog = soseq(prog,expr3);
prog = soseq(prog,coeff_1 - 1);
prog = soseq(prog,coeff_2 - 1);
prog = soseq(prog,coeff_3 - 1);

%% Call solver
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);

%% Get solution
SOLV = sosgetsol(prog,V)

% True Lyapunov function
%K = sequil(1)*log(sequil(1)) + sequil(2)*log(sequil(2)) + sequil(3)*log(sequil(3));
%trueLyap = x1 - sequil(1)*log(x1 + sequil(1)) + x2 - sequil(2)*log(x2 + sequil(2)) + K;

K = subs(SOLV,[x1,x2,x3],[sequil(1),sequil(2),sequil(3)]);
SOLV = SOLV-K;
SOLV = subs(SOLV,[x1, x2, x3],[x1-sequil(1), x2-sequil(2), x3-sequil(3)])
% Subsitute in log terms and reduce significant figures
SOLV3 = vpa(SOLV,4)

%% Plot Solution
% Calculate cross sections to plot 
SOLVXY = subs(SOLV3,x3,sequil(3));
SOLVXZ = subs(SOLV3,x2,sequil(2));
SOLVYZ = subs(SOLV3,x1,sequil(1));

% Plot solution and true function
axis1 = 0.1; axis2 = 3;

subplot(1,3,1)
fsurf(SOLVXY,[axis1 axis2])
title('x_{3} = q_{3}')

subplot(1,3,2)
fsurf(SOLVXZ,[axis1 axis2])
title('x_{2} = q_{2}')

subplot(1,3,3)
fsurf(SOLVYZ,[axis1 axis2])
title('x_{1} = q_{1}')

% 
% subplot(2,2,1)
% fsurf(trueLyap,[axis1 axis2])
% title('True')
% 
% subplot(2,2,2)
% fsurf(SOLV,[axis1 axis2])
% title('Solution')
% 
% subplot(2,2,3)
% fcontour(trueLyap,[axis1 axis2])
% title('True')
% 
% subplot(2,2,4)
% fcontour(SOLV,[axis1 axis2])
% title('Solution')
% 

% A_Dstable(:,:,1) =
% 
%    -0.8213   -0.6662   -0.8043
%    -1.1916   -0.7542    0.6389
%     0.6673   -0.9483   -0.8150
% 
% 
% A_Dstable(:,:,2) =
% 
%    -0.2852   -0.1610    0.6937
%     0.3070   -0.7816   -0.6416
%    -0.8487   -0.4971   -0.4969
% 
% 
% A_Dstable(:,:,3) =
% 
%    -1.1499   -0.5502   -1.0848
%    -1.0589   -1.7808    1.8048
%    -0.3987   -1.6770   -0.9149
% 
% 
% A_Dstable(:,:,4) =
% 
%    -1.1863    0.6078    0.1033
%    -0.6940   -0.6834    0.2419
%    -0.1298   -1.6213   -1.0654
% 
% 
% A_Dstable(:,:,5) =
% 
%    -0.3327   -1.2961    1.5898
%    -0.2433   -0.0110   -0.3184
%    -0.9175    1.4335   -1.4970
% 
% 
% A_Dstable(:,:,6) =
% 
%    -1.5353   -0.9973   -1.9327
%     0.2041   -1.7103    2.0088
%    -0.4879   -0.4088   -0.9142
% 
% 
% A_Dstable(:,:,7) =
% 
%    -0.6352    0.1582    0.8868
%    -0.2967   -2.3241   -0.1806
%    -1.0363   -0.0044   -0.6594
% 
% 
% A_Dstable(:,:,8) =
% 
%    -0.6923    0.8832   -0.6990
%    -0.2189   -1.0131    0.6556
%    -0.1391   -1.3959   -0.9278
% 
% 
% A_Dstable(:,:,9) =
% 
%    -0.7502    0.6120   -1.3384
%    -1.3215   -0.7652   -0.3953
%    -1.0110   -1.1910   -0.1395
%    
%    A_Dunstable(:,:,1) =
% 
%    -0.8939    1.2219   -0.9357
%    -1.0469    0.2974   -1.1144
%    -0.4512    0.1056   -0.6957
% 
% 
% A_Dunstable(:,:,2) =
% 
%    -1.9978   -1.1943    0.0574
%    -1.3830   -0.9422    1.2849
%    -0.2195   -0.2373    0.0097
% 
% 
% A_Dunstable(:,:,3) =
% 
%    -0.6919    1.3339   -0.4263
%    -0.7247    1.1175   -2.0412
%     0.1684   -0.1490   -1.8653
% 
% 
% A_Dunstable(:,:,4) =
% 
%    -0.9629    1.0583    1.0496
%     0.2524   -2.5791    0.0603
%    -0.8328    1.4834   -0.4918
% 
% 
% A_Dunstable(:,:,5) =
% 
%    -1.0460    1.2383    0.0411
%    -2.6477    0.7553    1.4488
%    -2.9632    0.3359   -0.2762

