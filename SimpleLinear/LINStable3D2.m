%% Implementing the Whirling Pendulum example in IEEE CDC 2002 by AP

%% Set up variables
clear all
syms x1 x2 x3
vars = [x1; x2; x3];

% Define parameters, values are D-stable
A = [   -0.3509    0.2008   -0.5486;
    -1.0364   -1.5921   -1.3065;
     0.9168    1.5111   -0.4038];

% Works for both D-stable and D-unstable,
% Just needs to be linearly stable
A = [     0.3232   -0.8915   -1.3291;
   -0.7888   -0.6541    1.1016;
    0.9632   -2.3251   -1.7300];

A = [-0.4637   -0.9426    1.4830;
    1.7063   -2.9739   -0.7103;
   -0.1755   -0.7583    0.3918];

X_star = [1;0.5;0.25];
%X_star = [0;0;0];
B = -A*X_star;   
sequil = X_star;

% Analytic solution for equilbrium
%sequil = [(c*(a-b))/(a*c-b*d), (a*(c-d))/(a*c-b*d)];
%sequil = A\B; %inv(A)*B;

% % Shifted vector field
X = [x1 + sequil(1); x2 + sequil(2); x3 + sequil(3)];
%f2 = diag(X)*B + diag(X)*(A*X);
f = B + A*X;
%f = A*X;

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

[prog,V_quad] = sossosvar(prog,vars);
[prog,V_quar] = sossosvar(prog,vars);
V_quar = subs(V_quar, [x1 x2 x3], [x1^1 x2^2 x3^3]);

V = coeff_1*x1 + coeff_2*x2 + coeff_3*x3 + coeff_4*log(x1 + sequil(1)) + ...
    coeff_5*log(x2 + sequil(2)) + coeff_6*log(x3 + sequil(3)) + coeff_7;
V = V_quad;%  + V_quar;

% Set up polynomials
Z1 = monomials(vars,0:1); 
Z2 = monomials(vars,0:2);

%[prog,V] = sospolyvar(prog,Z1,'wscoeff')
%[prog,p1] = sospolyvar(prog,Z2,'wscoeff');
%[prog,p2] = sospolyvar(prog,Z2,'wscoeff');
%[prog,p3] = sospolyvar(prog,Z2,'wscoeff');

%% Define SOSP constraints
% Constraint 1: -dV/dx*f >= 0
% Region D
%D1 = -x1 - sequil(1); D2 = -x2 - sequil(2); D3 = -x3 - sequil(3); 
expr2 = -(diff(V,x1)*f(1) + diff(V,x2)*f(2) + diff(V,x3)*f(3)); %+ p1*D1 + p2*D2 + p3*D3;
prog = sosineq(prog,expr2);

% Constraint 2: p >= 0
%prog = sosineq(prog,p1);
%prog = sosineq(prog,p2);
%prog = sosineq(prog,p3);

% Constraint 3: coeff_4 <= 0, coeff_5 <= 0 (d^2V/dx^2 <= 0)
Q = subs(jacobian(jacobian(V,[x1,x2,x3]),[x1,x2,x3]),[x1,x2,x3],[sequil(1),sequil(2),sequil(3)]);
prog = sosineq(prog,[x1,x2,x3]*Q*[x1;x2;x3] - 0.1*x1^2 - 0.1*x2^2 - 0.1*x3^2);

% Constraint 4: Fix parameters to match true Lyapunov function
%expr3 = coeff_5 + coeff_3*log(sequil(1)) + coeff_4*log(sequil(2));
%prog = soseq(prog,expr3);
%prog = soseq(prog,coeff_1 - 1);
%prog = soseq(prog,coeff_2 - 1);
%prog = soseq(prog,coeff_3 - 1);

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
% A_Dstable(:,:,1) =
% 
%    -0.3509    0.2008   -0.5486
%    -1.0364   -1.5921   -1.3065
%     0.9168    1.5111   -0.4038
% 
% 
% A_Dstable(:,:,2) =
% 
%    -1.2583   -1.6456   -0.3879
%     1.0313   -0.4606    0.6377
%     0.9409   -1.4811   -0.5230
% 
% 
% A_Dstable(:,:,3) =
% 
%    -0.2458   -2.0377   -0.5064
%     0.1529   -1.0190    0.2373
%     1.5018   -0.0050   -0.6539
% 
% 
% A_Dstable(:,:,4) =
% 
%    -1.5455   -0.9343   -0.2770
%     0.2595   -1.1344   -1.1950
%     1.1979    1.1093   -1.0448
% 
% 
% A_Dstable(:,:,5) =
% 
%    -0.4725    0.5456   -0.1723
%    -0.2239   -0.8378   -0.4366
%     1.1880    0.5127   -0.5380
% 
% 
% A_Dstable(:,:,6) =
% 
%    -0.2894    0.4036    0.5666
%    -0.7669   -0.4516    0.0761
%    -1.4322   -0.4147   -1.1970
% 
% 
% A_Dstable(:,:,7) =
% 
%    -0.5408   -0.9728   -0.2623
%     0.0210   -1.3699   -0.1865
%    -1.1046    0.7677   -0.6703
% 
% 
% A_Dstable(:,:,8) =
% 
%    -0.0322   -1.1824    0.3160
%     1.4931   -2.8696   -0.2952
%    -1.6067   -1.1268   -0.1625
% 
% 
% A_Dstable(:,:,9) =
% 
%    -1.3466    0.5119    0.9070
%    -0.3835   -1.5626   -1.1037
%    -0.5354    0.1622   -0.7989
% 
% 
% A_Dstable(:,:,10) =
% 
%    -0.5392    0.7638   -0.5480
%    -0.8957   -0.9675    0.2194
%    -0.5943   -0.2500   -1.4261
% 
% 
% A_Dstable(:,:,11) =
% 
%    -0.5197   -0.9593    0.1797
%     1.2413   -0.8550   -0.9292
%    -1.2108    0.2549   -1.3319
% 
% 
% A_Dstable(:,:,12) =
% 
%    -1.1423   -1.1602   -0.3801
%    -1.0968   -1.1245   -0.7114
%    -0.3178    0.0362   -1.7173

%%
% A_Dunstable(:,:,17) =
% 
%     0.3232   -0.8915   -1.3291
%    -0.7888   -0.6541    1.1016
%     0.9632   -2.3251   -1.7300
% 
% 
% A_Dunstable(:,:,18) =
% 
%     0.0581    0.2242   -0.8577
%     0.4172   -0.5938   -2.5517
%     0.9143    1.1478    0.0624
% 
% 
% A_Dunstable(:,:,19) =
% 
%    -1.3013    1.5938   -1.1706
%     0.2308   -0.2513   -0.6292
%     1.8415    0.9270   -0.5342
% 
% 
% A_Dunstable(:,:,20) =
% 
%     0.1994    0.0093   -1.7216
%     0.8374   -0.5164    0.0812
%     1.2339   -0.4224   -0.6816
% 
% 
% A_Dunstable(:,:,21) =
% 
%    -0.5723   -0.2955    0.8436
%    -0.7451   -0.4705    0.5108
%     0.1300    0.3421   -0.4251
% 
% 
% A_Dunstable(:,:,22) =
% 
%    -0.0603   -0.7064    0.5918
%     0.4263   -0.1703   -0.4820
%     0.3683    0.8732   -1.2133
% 
% 
% A_Dunstable(:,:,23) =
% 
%    -0.7639    0.1176    0.2406
%     0.4560   -0.5943   -0.0609
%    -1.1776    0.6062   -1.0172
% 
% 
% A_Dunstable(:,:,24) =
% 
%    -1.0920   -0.1558    0.6585
%     0.6363   -0.8388   -0.4628
%    -2.0440   -0.1499    0.3804
% 
% 
% A_Dunstable(:,:,25) =
% 
%    -1.2564    0.3312   -2.3020
%     0.5219   -1.7017    0.1512
%     0.1519    1.1068    0.5030
% 
% 
% A_Dunstable(:,:,26) =
% 
%    -1.1212    0.5406    0.2993
%     0.1903   -1.3742   -1.6730
%     0.6937    0.5414    0.7667
% 
% 
% A_Dunstable(:,:,27) =
% 
%    -1.2049    1.7661   -0.8899
%    -1.0664    0.2580    1.0001
%    -0.5991   -0.5064   -1.5395
% 
% 
% A_Dunstable(:,:,28) =
% 
%     0.2146    0.4156    0.1287
%    -0.8217   -0.2682   -0.6076
%     0.1918    0.8417   -0.6184
% 
% 
% A_Dunstable(:,:,29) =
% 
%    -1.5388   -1.3771    0.4074
%     1.5493    0.8331    0.5363
%    -1.1512   -0.3417   -3.4167
% 
% 
% A_Dunstable(:,:,30) =
% 
%     0.5263    0.5177    1.3294
%    -0.2722   -0.0779   -0.0783
%    -2.1815    0.5829   -1.0671
% 
% 
% A_Dunstable(:,:,31) =
% 
%    -1.6088   -1.3075   -1.6265
%    -0.6196   -0.8078    1.0413
%     1.9921    1.6468   -0.0557
% 
% 
% A_Dunstable(:,:,32) =
% 
%    -0.2514    0.8935   -1.4229
%    -0.4769    0.6176    0.0080
%     0.9471   -1.3986   -1.5126
% 
% 
% A_Dunstable(:,:,33) =
% 
%    -1.1163    0.6233    0.0984
%    -1.0471   -1.1641   -0.7008
%     1.2449    0.9288    0.1322
% 
% 
% A_Dunstable(:,:,34) =
% 
%    -1.4511    1.3166    0.0280
%    -0.2837   -0.1480   -1.1869
%     0.8883   -0.3388   -1.6157
% 
% 
% A_Dunstable(:,:,35) =
% 
%    -0.5602    2.4078    0.0589
%    -0.3690   -0.1843    0.3661
%     1.7606   -1.9412   -1.4257
% 
% 
% A_Dunstable(:,:,36) =
% 
%    -0.8386   -1.4047    2.0560
%     0.3196   -0.6901   -0.3820
%    -0.4506    0.3362    0.5249
% 
% 
% A_Dunstable(:,:,37) =
% 
%    -0.6774    0.3000   -0.8142
%     0.3838   -1.0950   -0.5934
%     0.9841    0.5427   -0.5796
% 
% 
% A_Dunstable(:,:,38) =
% 
%    -0.2446   -1.6127   -0.9360
%     0.3680   -0.8491   -0.1410
%     1.9897   -0.4385    0.8457
% 
% 
% A_Dunstable(:,:,39) =
% 
%     0.2447    1.1277   -0.1204
%    -1.1749   -0.7639    0.9828
%     0.8340    0.1109   -1.5063
% 
% 
% A_Dunstable(:,:,40) =
% 
%    -2.7657    1.8207   -0.1071
%     1.8663   -1.9597    1.1053
%    -0.9081   -0.5743   -1.0065
% 
% 
% A_Dunstable(:,:,41) =
% 
%    -0.9312    1.0108    0.3428
%     0.4398   -1.6851   -1.0467
%    -0.8468    0.4362   -0.6309


