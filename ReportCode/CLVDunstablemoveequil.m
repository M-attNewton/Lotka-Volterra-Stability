
%% Set up variables
clear all
k = 1;
%for c1 = -2:1:2
%for c2 = -2:1:2
%for c3 = -2:1:2
c1 = -0.1; c2 = -0.1; c3 = -0.1;
c1 = 0; c2 = 0; c3 = 0;
%clearvars -except c1 c2 c3 k test
syms x1 x2 x3
vars = [x1; x2; x3];

% Define parameters, values are not D-stable
A = [    -0.3792    0.0457    0.4302;
    -1.1050   -0.6037    0.3413;
    -0.1886    0.0455   -0.0847];

X_star = [1;0.5;0.25];
B = -A*X_star;

% Analytic solution for equilbrium
%sequil = [(c*(a-b))/(a*c-b*d), (a*(c-d))/(a*c-b*d)];
%sequil = X_star;

% % Shifted vector field
X = [x1; x2; x3];
%X = [x1 + sequil(1); x2 + sequil(2); x3 + sequil(3)];
f_quad = diag(X)*B + diag(X)*(A*X);
f = f_quad + [c1*x1*x2*x3; c2*x1*x2*x3; c3*x1*x2*x3];
f2 = B + A*X + [c1*x2*x3; c2*x1*x3; c3*x1*x2];

eqns = [f2(1) == 0; f2(2) == 0; f2(3) == 0]; 
equil = solve(eqns, [x1; x2; x3]); 

% Test linear stability
stability = zeros(size(equil.x1));
instability = stability;
for i = 1:length(equil.x1)
    
    % Calculate Jacobian for equilibrium
    jacob = jacobian([f2(1); f2(2); f2(3)], [x1; x2; x3]);
    jacobSubs(:,:,i) = subs(jacob, [x1; x2; x3], [equil.x1(i); equil.x2(i); equil.x3(i)]);
    
    % Find eigenvalues of Jacobian
    eigen(:,i) = eig(jacobSubs(:,:,i));
    
    % Check sign of eigenvalues to determine stability of equilibruim
    if real(double(subs(eigen(1,i)))) >= 0 || ...
       real(double(subs(eigen(2,i)))) >= 0 || ...
       real(double(subs(eigen(3,i)))) >= 0
   
        stability(i) = 0;
        
    else
        if vpa(equil.x1(i)) > 0 && vpa(equil.x2(i)) > 0 && vpa(equil.x3(i)) > 0
            stability(i) = 1;
        else
            stability(i) = 0;
        end
    end
    
    if stability(i) == 0 && vpa(equil.x1(i)) > 0 && vpa(equil.x2(i)) > 0 && vpa(equil.x3(i)) > 0
        instability(i) = 1;
    end
    
end

% Find relevant equilibria - assume only one for now
ind = find(stability);
if isempty(ind)
    disp('No stable points found');
    return
end

if size(ind,1) >= 2
    disp('Two stable points in positive orthant found')
    return
end

sequil = vpa([equil.x1(ind); equil.x2(ind); equil.x3(ind)]);

% Shift to equilibrium
f = subs(f, [x1; x2; x3], [x1+sequil(1); x2+sequil(2); x3+sequil(3)]);
f2 = subs(f2, [x1; x2; x3], [x1+sequil(1); x2+sequil(2); x3+sequil(3)]);

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
%V_quad = subs(V_quad, [x1;x2;x3], [x1 - sequil(1); x2 - sequil(2); x3 - sequil(3)]);
%V_quad = 0;

V = coeff_1*x1 + coeff_2*x2 + coeff_3*x3 + coeff_4*log(x1 + sequil(1)) + ...
    coeff_5*log(x2 + sequil(2)) + coeff_6*log(x3 + sequil(3)) + coeff_7 + V_quad;

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

expr2 = -(coeff_1*f(1) + coeff_2*f(2) + coeff_3*f(3) + coeff_4*f2(1) + coeff_5*f2(2) + coeff_6*f2(3)) ...
        -(diff(V_quad,x1)*f(1) + diff(V_quad,x2)*f(2) + diff(V_quad,x3)*f(3))... 
        + p1*D1 + p2*D2 + p3*D3;
%expr2 = -(diff(V,x1)*f2(1) + diff(V,x2)*f2(2) + diff(V,x3)*f2(3))...
      %  + p1*D1 + p2*D2 + p3*D3;
    %-(diff(V,x1)*x1*x2*x3 + diff(V,x2)*x1*x2*x3 + diff(V,x3)*x1*x2*x3);
    expr3 = -(diff(V,x1)*f(1) + diff(V,x2)*f(2) + diff(V,x3)*f(3)) + p1*D1 + p2*D2 + p3*D3;
    
test2 = -(coeff_1*f(1) + coeff_2*f(2) + coeff_3*f(3) + coeff_4*f2(1) + coeff_5*f2(2) + coeff_6*f2(3)) ...
        -(diff(V_quad,x1)*f(1) + diff(V_quad,x2)*f(2) + diff(V_quad,x3)*f(3));
test3 = -(diff(V,x1)*f(1) + diff(V,x2)*f(2) + diff(V,x3)*f(3));
test4 = -(coeff_1*f(1) + coeff_2*f(2) + coeff_3*f(3) + coeff_4*f2(1) + coeff_5*f2(2) + coeff_6*f2(3)) ...
        -(diff(V_quad,x1)*f(1) + diff(V_quad,x2)*f(2) + diff(V_quad,x3)*f(3))...
        +(diff(V,x1)*f(1) + diff(V,x2)*f(2) + diff(V,x3)*f(3));
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
solver_opt.solver = 'sedumi'; %sdpt3
prog = sossolve(prog,solver_opt);

%% Get solution
SOLV = sosgetsol(prog,V)
SOLVXY = subs(SOLV, x3, 0);
SOLVXZ = subs(SOLV, x2, 0);
SOLVYZ = subs(SOLV, x1, 0);
% test(c1+3,c2+3,c3+3) = SOLV;
% k = k + 1;
% 
% %end
% %end
% %end
% 
% g = 1;
% for i = 1:size(test,1)
% for j = 1:size(test,2)
% for k = 1:size(test,3)
% 
% temp = subs(test(i,j,k), [x1 x2 x3], [1 1 1]);
% 
% if abs(temp) > 0.00001
%     
%     result(g,1) = i-3;
%     result(g,2) = j-3;
%     result(g,3) = k-3;
%     g = g + 1;
%     
% end
% 
% end
% end
% end
% 
% % True Lyapunov function
% K = sequil(1)*log(sequil(1)) + sequil(2)*log(sequil(2));
% trueLyap = x1 - sequil(1)*log(x1 + sequil(1)) + x2 - sequil(2)*log(x2 + sequil(2)) + K;
% 
% % Plot solution and true function
% axis1 = -1; axis2 = 3;
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
% A_Dunstable(:,:,1) =
% 
%    -0.3792    0.0457    0.4302
%    -1.1050   -0.6037    0.3413
%    -0.1886    0.0455   -0.0847
% 
% 
% A_Dunstable(:,:,2) =
% 
%    -1.7925   -0.2878    0.5348
%    -1.0799   -1.2913    0.6905
%     0.4116   -0.8837   -0.5960
% 
% 
% A_Dunstable(:,:,3) =
% 
%    -0.4716    0.0370    0.7681
%    -0.3846   -1.3484    0.4667
%     0.1077   -0.0944   -2.0150
% 
% 
% A_Dunstable(:,:,4) =
% 
%    -2.6755    1.3400   -0.2294
%    -0.1028   -0.0963    0.2026
%     0.2621   -1.9300   -1.6367
% 
% 
% A_Dunstable(:,:,5) =
% 
%    -0.6757   -1.1762    0.7103
%    -0.3065   -0.5929    0.2673
%    -0.6708    0.1365   -1.5993
% 
% 
% A_Dunstable(:,:,6) =
% 
%    -1.1269    0.6280    1.0653
%    -0.4507    0.4956   -0.0922
%    -1.3070    1.6625    0.2304
% 
% 
% A_Dunstable(:,:,7) =
% 
%    -1.1764   -0.8197    0.1252
%     0.2041    0.1402   -1.7213
%    -0.1171    0.0227   -0.9846
% 
% 
% A_Dunstable(:,:,8) =
% 
%    -0.8699    0.2945   -0.4968
%     0.3507   -1.2428   -0.3141
%    -1.1285    0.2973   -1.2938
% 
% 
% A_Dunstable(:,:,9) =
% 
%    -0.5363    0.9657   -0.9589
%    -0.4468    0.3763   -0.0403
%    -0.5499    0.8848   -1.2158
% 
% 
% A_Dunstable(:,:,10) =
% 
%    -2.6550   -0.5580    1.3238
%    -0.7804   -0.5295    0.7652
%     0.4228   -1.1412   -0.9161
% 
% 
% A_Dunstable(:,:,11) =
% 
%    -1.4939    0.8202    2.0261
%    -0.1347   -1.8552    0.4454
%    -1.5987   -0.9512    0.5220
% 
% 
% A_Dunstable(:,:,12) =
% 
%    -0.5040    0.0866    0.3481
%    -0.7698   -0.0451    0.0459
%    -1.8794    1.6805   -1.2689
% 
% 
% A_Dunstable(:,:,13) =
% 
%    -1.3504    1.4942   -0.4342
%     0.3646   -1.8142    1.0863
%    -1.1681   -0.5653    0.8145
% 
% 
% A_Dunstable(:,:,14) =
% 
%    -1.4829    2.3111    0.8588
%    -0.6386   -1.5994    1.5277
%    -2.2144    0.8474    1.8957
% 
% 
% A_Dunstable(:,:,15) =
% 
%    -2.7248    1.0353    0.0529
%    -0.9496    0.4139    0.6629
%    -1.0442   -0.8136   -0.5480
% 
% 
% A_Dunstable(:,:,16) =
% 
%    -0.4756   -0.5604    0.8135
%    -0.4742   -1.1877    0.3165
%    -0.7708    1.8776   -1.0114