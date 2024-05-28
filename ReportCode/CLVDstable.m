
%% Set up variables
clear all
k = 1;
for c1 = -2:1:2
for c2 = -2:1:2
for c3 = -2:1:2
clearvars -except c1 c2 c3 k test
syms x1 x2 x3
vars = [x1; x2; x3];

% Define parameters, values are not D-stable
A = [    -0.7186    0.4736    0.1267;
    -1.0358    0.3034   -0.0877;
    -0.1814   -0.7840   -2.3431];

X_star = [1;0.5;0.25];
B = -A*X_star;

% Analytic solution for equilbrium
%sequil = [(c*(a-b))/(a*c-b*d), (a*(c-d))/(a*c-b*d)];
sequil = X_star;

% % Shifted vector field
X = [x1 + sequil(1); x2 + sequil(2); x3 + sequil(3)];
f = diag(X)*B + diag(X)*(A*X);
f2 = f - 1*[c1*(x1 + sequil(1))*x2*x3; c2*x1*(x2 + sequil(2))*x3; c3*x1*x2*(x3 + sequil(3))];

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
V_quad = 0;

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
solver_opt.solver = 'sedumi'; %sdpt3
prog = sossolve(prog,solver_opt);

%% Get solution
SOLV = sosgetsol(prog,V)

test(c1+3,c2+3,c3+3) = SOLV;
k = k + 1;

end
end
end

g = 1;
for i = 1:size(test,1)
for j = 1:size(test,2)
for k = 1:size(test,3)

temp = subs(test(i,j,k), [x1 x2 x3], [1 1 1]);

if abs(temp) > 0.00001
    
    result(g,1) = i-3;
    result(g,2) = j-3;
    result(g,3) = k-3;
    g = g + 1;
    
end

end
end
end

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

% 
% A_Dunstable(:,:,1) =
% 
%    -0.7186    0.4736    0.1267
%    -1.0358    0.3034   -0.0877
%    -0.1814   -0.7840   -2.3431
% 
% 
% A_Dunstable(:,:,2) =
% 
%    -0.2764    0.3599   -0.1036
%    -0.0406   -1.3129    0.2466
%    -1.7054    1.1447   -0.8646
% 
% 
% A_Dunstable(:,:,3) =
% 
%    -0.3103    0.6238   -0.1895
%    -0.0548   -1.0824    0.6363
%    -0.2763   -1.1093    0.4247
% 
% 
% A_Dunstable(:,:,4) =
% 
%    -1.8705    1.8351   -0.4970
%    -1.0047    0.3671   -0.1909
%    -1.0818    0.3502   -0.5425