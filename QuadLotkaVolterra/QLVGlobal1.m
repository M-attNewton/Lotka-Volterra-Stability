%% Implementing the Whirling Pendulum example in IEEE CDC 2002 by AP
clear all
for epsilon = 0:0.01:0.2
%% Set up variables
clearvars -except epsilon
syms x1 x2 u1 u2 %u1 = ln(x1), u2 = ln(x2)
vars = [x1; x2; u1; u2];

% Define parameters
a = 3; b = 2; c = 4; d = 3;

% Modified vector field
x1_dot = a*(x1+(2/3)) - a*(x1+(2/3))^2 - b*(x1+(2/3))*(x2+(1/2));
x2_dot = c*(x2+(1/2)) - c*(x2+(1/2))^2 - d*(x1+(2/3))*(x2+(1/2));
f = [x1_dot;
     x2_dot;
     ((x1+(2/3))^-1)*x1_dot;
     ((x2+(1/2))^-1)*x2_dot];

% Sub in exp(xy)
% x1_dot = a*(x1+(2/3)) - a*(x1+(2/3))^2 - b*(1 + (u1+u2) + 0.5*(u1+u2)^2 + (1/6)*(u1+u2)^3);
% x2_dot = c*(x2+(1/2)) - c*(x2+(1/2))^2 - d*(1 + (u1+u2) + 0.5*(u1+u2)^2 + (1/6)*(u1+u2)^3);
% f = [x1_dot;
%      x2_dot;
%      ((x1+(2/3))^-1)*x1_dot;
%      ((x2+(1/2))^-1)*x2_dot];

 
% Initialize the sum of squares program
prog = sosprogram(vars);

% The Lyapunov function V(x): 
Z1 = monomials(vars,0:1);
%Z1 = [x1; x2; u1; u2; 1];
Z2 = monomials(vars,0:2);
[prog,V] = sospolyvar(prog,Z1,'wscoeff')
[prog,p1] = sospolyvar(prog,Z2,'wscoeff')
[prog,p2] = sospolyvar(prog,Z2,'wscoeff')

%% Define SOSP constraints
% Constraint 1 : V(x) - (x1^2 + x2^2) >= 0
%prog = sosineq(prog,V - (x1^2 + x2^2));
%prog = sosineq(prog,V - epsilon*(x1^2 + x2^2 + u1^2 + u2^2)); % + u1^2 + u2^2

% Constraint 2: -dV/dx*f >= 0
% Region D
D1 = -x1 - (2/3); D2 = -x2 - (1/2);
%D = (x1^2 + x2^2 - 0.3^2);
%D = -x1 - x2;
%((x1-(2/3))^2/(0.6^2) + (x2-(1/2))^2/(0.45^2));
%((x1-(2/3))^2 + (x2-(1/2))^2 - 0.4^2); %p*(x1^2 + x2^2 - 0.4^2);

expr2 = -(diff(V,x1)*f(1) + diff(V,x2)*f(2) + diff(V,u1)*f(3) + diff(V,u2)*f(4)) +... p1*D;
     p1*D1 + p2*D2;
prog = sosineq(prog,expr2);

% Constraint 3: p >= 0
prog = sosineq(prog,p1);% - epsilon*(x1^2 + x2^2 + u1^2 + u2^2)); %+ u1^2 + u2^2
prog = sosineq(prog,p2);% - epsilon*(x1^2 + x2^2 + u1^2 + u2^2));
prog = sosineq(prog,-coeff_4);
prog = sosineq(prog,-coeff_5);
%prog = soseq(prog,coeff_1 + coeff_4*log(2/3) + coeff_5*log(1/2));


%% Call solver
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);

%% Get solution
SOLV = sosgetsol(prog,V)
%echo off;

if SOLV ~= 0
    break
end

end

% subs(x1,log(x1));
% subs(x2,log(x2));

SOLV2 = subs(SOLV, u1, log(x1));
SOLV3 = subs(SOLV2, u2, log(x2))
SOLV3 = vpa(SOLV3,4)

% True Lyapunov function
K = -1.7836;
K = 0;
%trueLyap = d*x1 - c*log(x1) + b*x2 - a*log(x2) + K
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


