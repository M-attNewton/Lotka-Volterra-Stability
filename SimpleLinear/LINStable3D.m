%% Implementing the Whirling Pendulum example in IEEE CDC 2002 by AP

%% Set up variables
clear all
syms x1 x2 x3
vars = [x1; x2; x3];

% Define parameters, values are D-stable
% A = [    -1.1499   -0.5502   -1.0848;
%     -1.0589   -1.7808    1.8048;
%     -0.3987   -1.6770   -0.9149];
% 
% A = [-1.1863    0.6078    0.1033;
%     -0.6940   -0.6834    0.2419;
%     -0.1298   -1.6213   -1.0654];

A = [-1.7648    1.8522   -0.7684;
   -0.5165   -1.1294   -1.0110;
   -0.8923    0.9811   -1.4073];

A = [-0.6203   -1.6592   -0.0660;
     1.7346    0.1874    0.4172;
    -1.1939    2.4415   -0.4887];

A = [    -0.3313   -0.4280   -0.1974;
     0.9467   -0.7081    1.1106;
     0.1764    0.2788   -1.1591];

X_star = [0;0;0];
%B = -A*X_star;   
sequil = X_star;

% Analytic solution for equilbrium
%sequil = [(c*(a-b))/(a*c-b*d), (a*(c-d))/(a*c-b*d)];
%sequil = A\B; %inv(A)*B;

% % Shifted vector field
X = [x1 + sequil(1); x2 + sequil(2); x3 + sequil(3)];
%f2 = diag(X)*B + diag(X)*(A*X);
f = A*X;

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
V = V_quad %+ V_quar;

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
% expr2 = ...
% -(diff(V,x1)*(A(1,1)*(x1+sequil(1)) + A(1,2)*(x2+sequil(2)) + A(1,3)*(x3+sequil(3))) + ...
% diff(V,x2)*(A(2,1)*(x1+sequil(1)) + A(2,2)*(x2+sequil(2)) + A(2,3)*(x3+sequil(3)))  + ...
% diff(V,x3)*(A(3,1)*(x1+sequil(1)) + A(3,2)*(x2+sequil(2)) + A(3,3)*(x3+sequil(3))) );
% expr2 = ...
% -(coeff_1*A(1,1)*(x1+sequil(1)) + coeff_4*A(1,1) ...
% +(coeff_1 + coeff_4/(x1 + sequil(1)))*(A(1,2)*(x2+sequil(2)) + A(1,3)*(x3+sequil(3))));%+...
% (coeff_2 + coeff_5/(x2 + sequil(2)))*(A(2,1)*(x1+sequil(1)) + A(2,2)*(x2+sequil(2)) + A(2,3)*(x3+sequil(3)))  + ...
% (coeff_3 + coeff_6/(x3 + sequil(3)))*(A(3,1)*(x1+sequil(1)) + A(3,2)*(x2+sequil(2)) + A(3,3)*(x3+sequil(3))) );
prog = sosineq(prog,expr2);

% Constraint 2: p >= 0
%prog = sosineq(prog,p1);
%prog = sosineq(prog,p2);
%prog = sosineq(prog,p3);

% Constraint 3: coeff_4 <= 0, coeff_5 <= 0 (d^2V/dx^2 <= 0)
Q = subs(jacobian(jacobian(V,[x1,x2,x3]),[x1,x2,x3]),[x1,x2,x3],[0,0,0]);
prog = sosineq(prog,[x1,x2,x3]*Q*[x1;x2;x3] - 0.1*x1^2 - 0.1*x2^2 - 0.1*x3^2);
%prog = sosineq(prog,-coeff_4);
%prog = sosineq(prog,-coeff_5);
%prog = sosineq(prog,-coeff_6);

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
%    -1.7648    1.8522   -0.7684
%    -0.5165   -1.1294   -1.0110
%    -0.8923    0.9811   -1.4073
% 
% 
% A_Dstable(:,:,2) =
% 
%    -0.8151   -0.1243   -0.1246
%     0.5654   -1.7749    0.5027
%     0.4656   -0.3116   -1.9001
% 
% 
% A_Dstable(:,:,3) =
% 
%    -0.0370    2.3223   -0.7280
%    -0.4268   -0.9593    0.1406
%    -1.4538   -0.1347   -1.1736
% 
% 
% A_Dstable(:,:,4) =
% 
%    -0.2778    0.0494    0.4446
%    -0.2194   -0.8299   -0.1308
%    -0.9720    1.9691   -0.0434
% 
% 
% A_Dstable(:,:,5) =
% 
%    -3.3909    0.2465    1.2731
%    -1.4234   -0.1581    1.7882
%    -0.4887   -0.3087   -0.6143
% 
% 
% A_Dstable(:,:,6) =
% 
%    -1.3407   -0.0404   -0.3278
%    -1.6239   -0.6972   -0.2779
%     0.1462   -0.3640   -1.2061
% 
% 
% A_Dstable(:,:,7) =
% 
%    -0.5099   -2.3851   -0.0319
%     1.7136   -0.3591   -0.7876
%    -0.6633   -0.7912   -1.1896
% 
% 
% A_Dstable(:,:,8) =
% 
%    -0.2109   -0.4029   -1.2545
%     0.4166   -0.0351   -0.6324
%     0.6221    0.6392   -1.4184
% 
% 
% A_Dstable(:,:,9) =
% 
%    -0.0378    2.1527    0.3084
%    -0.8081   -0.5239    0.4314
%    -0.5052   -0.1320   -0.5348
% 
% 
% A_Dstable(:,:,10) =
% 
%    -0.0751   -0.8640   -1.1832
%     0.7227   -0.3986    0.2825
%     1.1678   -2.0464   -1.6315
% 
% 
% A_Dstable(:,:,11) =
% 
%    -1.6320   -1.1523    0.1716
%    -0.7921   -0.8471   -1.8590
%    -0.7552    0.4734   -0.9055
% 
% 
% A_Dstable(:,:,12) =
% 
%    -1.8380    0.3135   -0.1136
%    -0.8367   -0.9488   -0.2124
%     1.2564   -0.4468   -0.7636
% 
% 
% A_Dstable(:,:,13) =
% 
%    -0.0634   -0.2584   -0.1880
%    -0.1477   -0.7728    2.6208
%     0.6411   -1.2527   -1.8164
% 
% 
% A_Dstable(:,:,14) =
% 
%    -0.6382    0.2633    0.5816
%    -1.7191   -0.5461   -0.7973
%    -0.1667   -0.0204   -0.8889
% 
% 
% A_Dstable(:,:,15) =
% 
%    -1.3220   -0.1235    0.0426
%     0.4438   -1.0664   -0.3994
%    -1.4314    0.9999   -0.2313
% 
% 
% A_Dstable(:,:,16) =
% 
%    -1.5448   -0.3645    0.6172
%     0.5898   -0.1458    1.3128
%    -1.7147   -0.7979   -0.6780
% 
% 
% A_Dstable(:,:,17) =
% 
%    -0.4052    1.3407    1.6136
%    -1.6458   -0.4223    2.1548
%    -1.1563   -0.6362   -0.0695
% 
% 
% A_Dstable(:,:,18) =
% 
%    -1.0517   -0.0951   -0.8621
%    -0.7545   -0.0814   -0.2389
%    -0.4510   -0.1578   -1.6840
% 
% 
% A_Dstable(:,:,19) =
% 
%    -1.5218   -1.6963   -0.4720
%    -0.2746   -0.2387    0.8038
%    -1.2686   -1.8449   -0.8880
%    
%    
% A_Dunstable(:,:,19) =
% 
%    -2.1473    1.2121    0.7003
%    -0.3437   -0.8147    1.2308
%     0.9813   -0.1023   -1.2101
% 
% 
% A_Dunstable(:,:,20) =
% 
%    -0.3313   -0.4280   -0.1974
%     0.9467   -0.7081    1.1106
%     0.1764    0.2788   -1.1591
% 
% 
% A_Dunstable(:,:,21) =
% 
%    -0.6203   -1.6592   -0.0660
%     1.7346    0.1874    0.4172
%    -1.1939    2.4415   -0.4887
% 
% 
% A_Dunstable(:,:,22) =
% 
%     0.4142   -1.1827    1.2928
%     0.6329   -1.7690    0.6219
%    -0.2809   -1.0881   -0.7493
% 
% 
% A_Dunstable(:,:,23) =
% 
%    -0.9778    0.6928    0.0721
%     0.1116    0.1859   -1.2077
%     0.4234    0.6442   -0.9371
% 
% 
% A_Dunstable(:,:,24) =
% 
%     0.1876    0.7340   -1.0286
%    -2.7001   -1.0431   -1.0387
%    -0.3644   -0.2882   -1.1295
% 
% 
% A_Dunstable(:,:,25) =
% 
%    -0.2400    1.3505   -0.2132
%    -0.4394    0.1435    0.3282
%    -0.0974    0.1881   -0.1341
% 
% 
% A_Dunstable(:,:,26) =
% 
%    -0.8463   -0.2827   -1.0921
%     0.5672   -1.2259   -0.4975
%     1.0109    0.7173    0.4580
% 
% 
% A_Dunstable(:,:,27) =
% 
%     0.0526   -0.4961   -0.4094
%     2.0112   -1.1423   -1.8756
%    -0.6444    1.2852   -0.3233
% 
% 
% A_Dunstable(:,:,28) =
% 
%    -1.5731    1.4770    1.3411
%     0.7713   -0.7686    0.1342
%     0.4483   -0.4347   -1.2922
% 
% 
% A_Dunstable(:,:,29) =
% 
%     0.1693   -1.0101   -1.5731
%     1.5973   -2.2046   -0.9963
%    -1.4994    1.1948   -1.0953
% 
% 
% A_Dunstable(:,:,30) =
% 
%    -0.7192    0.2379    0.9629
%     1.2128   -0.4250   -0.9694
%    -0.6851    0.2586   -1.5043
% 
% 
% A_Dunstable(:,:,31) =
% 
%     0.3944    0.6024    1.6538
%    -0.5245   -0.6479   -1.0292
%    -1.3494    0.7883   -0.6849
% 
% 
% A_Dunstable(:,:,32) =
% 
%    -0.7015   -0.4542    0.0589
%     1.7860    0.6509    0.4431
%    -2.1122   -0.6980   -0.7750
% 
% 
% A_Dunstable(:,:,33) =
% 
%     0.0479    0.5375   -0.2291
%    -0.4971   -0.9455   -0.5182
%     0.5110   -0.6849   -1.2680


%%

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

