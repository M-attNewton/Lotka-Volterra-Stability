
%% Set up variables
clear all
syms x1 x2 
vars = [x1; x2];

% Define parameters, values are not D-stable
A = [-0.6825    0.6237;
    -0.7868    0.2040];

X_star = [1;0.5];
B = -A*X_star;

% Analytic solution for equilbrium
%sequil = [(c*(a-b))/(a*c-b*d), (a*(c-d))/(a*c-b*d)];
sequil = X_star;

% % Shifted vector field
X = [x1 + sequil(1); x2 + sequil(2)];
f = diag(X)*B + diag(X)*(A*X);
%x1_dot = a*(x1+(sequil(1))) - a*(x1+(sequil(1)))^2 - b*(x1+(sequil(1)))*(x2+(sequil(2)));
%x2_dot = c*(x2+(sequil(2))) - c*(x2+(sequil(2)))^2 - d*(x1+(sequil(1)))*(x2+(sequil(2)));
% f2 = [x1_dot; x2_dot];

% Initialize the sum of squares program
prog = sosprogram(vars);

% The Lyapunov function V(x): 
[prog,coeff_1] = sospolyvar(prog,sym(1));
[prog,coeff_2] = sospolyvar(prog,sym(1));
[prog,coeff_3] = sospolyvar(prog,sym(1));
[prog,coeff_4] = sospolyvar(prog,sym(1));
[prog,coeff_5] = sospolyvar(prog,sym(1));

V = coeff_1*x1 + coeff_2*x2 + coeff_3*log(x1 + sequil(1)) + ...
    coeff_4*log(x2 + sequil(2)) + coeff_5;

% Set up polynomials
Z1 = monomials(vars,0:1); 
Z2 = monomials(vars,0:2);

%[prog,V] = sospolyvar(prog,Z1,'wscoeff')
[prog,p1] = sospolyvar(prog,Z2,'wscoeff');
[prog,p2] = sospolyvar(prog,Z2,'wscoeff');

%% Define SOSP constraints
% Constraint 1: -dV/dx*f >= 0
% Region D
D1 = -x1 - sequil(1); D2 = -x2 - sequil(2);
expr2 = -(diff(V,x1)*f(1) + diff(V,x2)*f(2)) + p1*D1 + p2*D2;
prog = sosineq(prog,expr2);

% Constraint 2: p >= 0
prog = sosineq(prog,p1);
prog = sosineq(prog,p2);

% Constraint 3: coeff_4 <= 0, coeff_5 <= 0 (d^2V/dx^2 >= 0)
prog = sosineq(prog,-coeff_3);
prog = sosineq(prog,-coeff_4);

% Constraint 4: Fix parameters to match true Lyapunov function
% expr3 = coeff_5 + coeff_3*log(sequil(1)) + coeff_4*log(sequil(2));
% prog = soseq(prog,expr3);
 prog = soseq(prog,coeff_1 - 1);
 prog = soseq(prog,coeff_2 - 1);

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

%%
% A_Dstable(:,:,25) =
% 
%    -1.0544   -0.4945
%    -0.3795   -0.4934
% 
% 
% A_Dstable(:,:,26) =
% 
%    -2.3213   -0.0037
%    -0.0977   -1.9589
% 
% 
% A_Dstable(:,:,27) =
% 
%    -0.4142    0.3384
%    -2.1983   -0.1562
% 
% 
% A_Dstable(:,:,28) =
% 
%    -1.1838   -1.8295
%     0.3984   -2.0730
% 
% 
% A_Dstable(:,:,29) =
% 
%    -0.2859   -1.2384
%     0.1932   -0.9297
% 
% 
% A_Dstable(:,:,30) =
% 
%    -1.5149   -0.8619
%    -0.0824   -1.0650
% 
% 
% A_Dstable(:,:,31) =
% 
%    -1.2868    0.4337
%    -0.9881   -0.5247
% 
% 
% A_Dstable(:,:,32) =
% 
%    -0.6019    1.2034
%    -0.9900   -0.1599
% 
% 
% A_Dstable(:,:,33) =
% 
%    -0.5694   -0.4956
%    -0.0803   -1.1126
% 
% 
% A_Dstable(:,:,34) =
% 
%    -0.3120    0.2361
%    -0.4026   -0.6295
% 
% 
% A_Dstable(:,:,35) =
% 
%    -1.5161   -1.0751
%    -1.1640   -1.2036
% 
% 
% A_Dstable(:,:,36) =
% 
%    -0.5466   -0.1738
%    -1.5704   -2.2032
% 
% 
% A_Dstable(:,:,37) =
% 
%    -1.7197    0.8070
%    -0.1986   -1.5311
% 
% A_Dunstable(:,:,1) =
% 
%    -1.2513    0.2948
%    -0.3687    0.0462
% 
% 
% A_Dunstable(:,:,2) =
% 
%    -0.8948    0.0283
%     0.4067   -2.2471
% 
% 
% A_Dunstable(:,:,3) =
% 
%    -1.0765    1.0410
%    -0.7334    0.6742
% 
% 
% A_Dunstable(:,:,4) =
% 
%    -1.3286    0.9334
%    -1.2502    0.8687
% 
% 
% A_Dunstable(:,:,5) =
% 
%    -1.9338    0.0596
%    -0.5555    0.0128
% 
% 
% A_Dunstable(:,:,6) =
% 
%     0.2730   -1.4828
%     0.3023   -1.2307
% 
% 
% A_Dunstable(:,:,7) =
% 
%    -0.2272    0.2366
%     0.3135   -0.7517
% 
% 
% A_Dunstable(:,:,8) =
% 
%    -0.6097    0.3042
%     0.0197   -0.7878
% 
% 
% A_Dunstable(:,:,9) =
% 
%    -1.1542    1.8291
%     0.0596   -0.8947
% 
% 
% A_Dunstable(:,:,10) =
% 
%    -2.0852    0.3182
%     0.3156   -0.6826
% 
% 
% A_Dunstable(:,:,11) =
% 
%    -1.0057    1.4934
%     0.0222   -0.3391
% 
% 
% A_Dunstable(:,:,12) =
% 
%    -1.1411    0.9544
%     0.1863   -2.1928
% 
% 
% A_Dunstable(:,:,13) =
% 
%    -0.6825    0.6237
%    -0.7868    0.2040
% 
% 
% A_Dunstable(:,:,14) =
% 
%    -1.7063    0.3066
%     0.5488   -1.1223
% 
% 
% A_Dunstable(:,:,15) =
% 
%    -1.1654    0.2788
%     0.2138   -1.1990
