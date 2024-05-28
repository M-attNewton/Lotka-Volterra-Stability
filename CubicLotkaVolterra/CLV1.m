%% 3D vector field with cubic terms 

%% Set up variables
clear all
syms x y z u v w % where: u = ln(x), v = ln(y), w = ln(z)
vars = [x; y; z; u; v; w];

% Define parameters
A = ones(3); b = ones(3,1); c = ones(3,1); % Set all as 1s for now

% True vector field
x_dot = x*(b(1) - A(1,1)*x - A(1,2)*y - A(1,3)*z + c(1)*y*z);
y_dot = y*(b(2) - A(2,1)*x - A(2,2)*y - A(2,3)*z + c(2)*x*z);
z_dot = z*(b(3) - A(3,1)*x - A(3,2)*y - A(3,3)*z + c(3)*x*y);

% Remove zero terms to speed up simultanious equation solver
x_dot = (b(1) - A(1,1)*x - A(1,2)*y - A(1,3)*z + c(1)*y*z);
y_dot = (b(2) - A(2,1)*x - A(2,2)*y - A(2,3)*z + c(2)*x*z);
z_dot = (b(3) - A(3,1)*x - A(3,2)*y - A(3,3)*z + c(3)*x*y);

% Calculate equilibrium points
eqns = [x_dot == 0; y_dot == 0; z_dot == 0];
equil = solve(eqns);

% Test linear stability
for i = 1:length(equil.x)
    
    % Repeat state vector field to avoid errors
    x_dot = x*(b(1) - A(1,1)*x - A(1,2)*y - A(1,3)*z + c(1)*y*z);
    y_dot = y*(b(2) - A(2,1)*x - A(2,2)*y - A(2,3)*z + c(2)*x*z);
    z_dot = z*(b(3) - A(3,1)*x - A(3,2)*y - A(3,3)*z + c(3)*x*y);
    
    % Calculate Jacobian for equilibrium
    jacob = jacobian([x_dot; y_dot; z_dot], [x; y; z]);
    jacobSubs(:,:,i) = subs(jacob, [x; y; z], [equil.x(i); equil.y(i); equil.z(i)]);
    
    % Find eigenvalues of Jacobian
    eigen(:,i) = eig(jacobSubs(:,:,i));
    
    % Check sign of eigenvalues to determine stability of equilibruim
    if real(double(subs(eigen(1,i)))) >= 0 || ...
       real(double(subs(eigen(2,i)))) >= 0 || ...
       real(double(subs(eigen(3,i)))) >= 0
        stability(i) = 0;
    else
        stability(i) = 1;
    end
    
end

% Find relevant equilibria - assume only one for now
ind = find(stability);
sequil = [equil.x(ind); equil.y(ind); equil.z(ind)];

% Shift vector field so equilibrium is at zero
x_dot = (x+sequil(1))*(b(1) - A(1,1)*(x+sequil(1)) - A(1,2)*(y+sequil(2)) ...
    - A(1,3)*(z+sequil(3)) + c(1)*(y+sequil(2))*(z+sequil(3)));
y_dot = (y+sequil(2))*(b(2) - A(2,1)*(x+sequil(1)) - A(2,2)*(y+sequil(2)) ...
    - A(2,3)*(z+sequil(3)) + c(2)*(x+sequil(1))*(z+sequil(3)));
z_dot = (z+sequil(3))*(b(3) - A(3,1)*(x+sequil(1)) - A(3,2)*(y+sequil(2)) ...
    - A(3,3)*(z+sequil(3)) + c(3)*(x+sequil(1))*(y+sequil(2)));

f = [x_dot; y_dot; z_dot; ...
     ((x+sequil(1))^-1)*x_dot; ((y+sequil(2))^-1)*y_dot; ((z+sequil(3))^-1)*z_dot];
 
% Initialize the sum of squares program
prog = sosprogram(vars);

% The Lyapunov function V(x): 
Z1 = monomials(vars,0:2); 
%Z1 = [x*y; x*z; y*z; x; y; z; u; v; w; 1];
Z1 = [x^2; y^2; z^2; x*y; x*z; y*z; x; y; z; u; v; w; 1];
Z2 = monomials(vars,0:2);
[prog,V] = sospolyvar(prog,Z1,'wscoeff')
[prog,p1] = sospolyvar(prog,Z2,'wscoeff')
[prog,p2] = sospolyvar(prog,Z2,'wscoeff')
[prog,p3] = sospolyvar(prog,Z2,'wscoeff')
%[prog,Q] = sospolymatrixvar(prog,monomials([x;y;z],0),[3 3],'symmetric')

%% Define SOSP constraints
%prog = sosineq(prog,V - (x^2 + y^2 + z^2));
% Constraint 1: -dV/dx*f >= 0
% Region D
D1 = -x - sequil(1); D2 = -y - sequil(2); D3 = -z - sequil(3);

expr2 = -(diff(V,x)*f(1) + diff(V,y)*f(2) + diff(V,z)*f(3) + ...
          diff(V,u)*f(4) + diff(V,v)*f(5) + diff(V,w)*f(6)) + ...
          p1*D1 + p2*D2 + p3*D3;
prog = sosineq(prog,expr2);

% Constraint 2: p >= 0
prog = sosineq(prog,p1);
prog = sosineq(prog,p2);
prog = sosineq(prog,p3);

% % Constraint 3: coeff_10 + 2*coeff_1*(sequil(1))^2 <= 0, coeff_11 <= 0, coeff_12 <= 0
%prog = sosineq(prog,-coeff_10 + 2*coeff_1*(sequil(1) + x)^2);
%prog = sosineq(prog,-coeff_11 + 2*coeff_2*(sequil(2) + y)^2);
%prog = sosineq(prog,-coeff_12 + 2*coeff_3*(sequil(3) + z)^2);
Q = [coeff_1, 0.5*coeff_4, 0.5*coeff_5; ...
     0.5*coeff_4, coeff_2, 0.5*coeff_6; ...
     0.5*coeff_5, 0.5*coeff_6, coeff_3];
matexpr1 = 2*Q + [(-coeff_10/(sequil(1))^2), 0, 0;...
                0, (-coeff_11/(sequil(2))^2), 0;...
                0, 0, (-coeff_12/(sequil(3))^2)];
prog = sosmatrixineq(prog,matexpr1,'quadraticMineq');

% prog = sosineq(prog,-coeff_10);
% prog = sosineq(prog,-coeff_11);
% prog = sosineq(prog,-coeff_12);

% Constraint 4: Fix parameters to match true Lyapunov function
% At equilbrium Lyaponuv function must be zero
expr3 = coeff_1*sequil(1)^2 + coeff_4*sequil(1)*sequil(2) + coeff_5*sequil(1)*sequil(3) + ...
    coeff_7*sequil(1) + coeff_2*sequil(2)^2 + coeff_6*sequil(2)*sequil(3) + ...
    coeff_8*sequil(2) + coeff_3*sequil(3)^2 + coeff_9*sequil(3) + coeff_13 + ...
    coeff_10*log(sequil(1)) + coeff_11*log(sequil(2)) + coeff_12*log(sequil(3));
prog = soseq(prog,expr3);

expr4 = coeff_7 + coeff_10/sequil(1) + 2*coeff_1*sequil(1) + coeff_4*sequil(2) ...
    + coeff_5*sequil(3);
prog = soseq(prog,expr4);

expr5 = coeff_8 + coeff_11/sequil(2) + 2*coeff_2*sequil(2) + coeff_4*sequil(1) ...
    + coeff_6*sequil(3);
prog = soseq(prog,expr5);

expr6 = coeff_9 + coeff_12/sequil(3) + 2*coeff_3*sequil(3) + coeff_5*sequil(1) ...
    + coeff_6*sequil(2);
prog = soseq(prog,expr6);

% prog = soseq(prog,coeff_7 - 1);
% prog = soseq(prog,coeff_8 - 1);
% prog = soseq(prog,coeff_9 - 1);

%% Call solver
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);

%% Get solution
SOLV = sosgetsol(prog,V)

% Subsitute in log terms and reduce significant figures
SOLV2 = subs(SOLV, u, log(x));
SOLV2 = subs(SOLV2, v, log(y));
SOLV2 = subs(SOLV2, w, log(z));
SOLV3 = vpa(SOLV2,4)

%% Plot Solution
% Calculate cross sections to plot 
SOLVXY = subs(SOLV3,z,sequil(3));
SOLVXZ = subs(SOLV3,y,sequil(2));
SOLVYZ = subs(SOLV3,x,sequil(1));

% Define axis
axis1 = 0.1; axis2 = 5;
figure
fsurf(SOLVXY,[axis1 axis2])

figure
fcontour(SOLVXY)

% subplot(1,2,1)
% fsurf(SOLVXY,[axis1 axis2])
% 
% subplot(1,2,2)
% fcontour(SOLVXY,[axis1 axis2])
% 
% subplot(2,3,1)
% fsurf(SOLVXY,[axis1 axis2])
% title('X-Y Plane, Z at Equilbruim')
% 
% subplot(2,3,2)
% fsurf(SOLVXZ,[axis1 axis2])
% title('X-Z Plane, Y at Equilbruim')
% 
% subplot(2,3,3)
% fsurf(SOLVYZ,[axis1 axis2])
% title('Y-Z Plane, X at Equilbruim')
% 
% subplot(2,3,4)
% fcontour(SOLVXY,[axis1 axis2])
% title('X-Y Plane, Z at Equilbruim')
% 
% subplot(2,3,5)
% fcontour(SOLVXZ,[axis1 axis2])
% title('X-Z Plane, Y at Equilbruim')
% 
% subplot(2,3,6)
% fcontour(SOLVYZ,[axis1 axis2])
% title('Y-Z Plane, X at Equilbruim')
% 


%% JUNK
    
%     x_dot = x*(a1 - a2*x - a3*y - a4*z + a5*y*z);
%     y_dot = x*(b1 - b2*x - b3*y - b4*z + b5*x*z);
%     z_dot = x*(c1 - c2*x - c3*y - c4*z + c5*x*y);

%a1 a2 a3 a4 a5 b1 b2 b3 b4 b5 c1 c2 c3 c4 c5


% % Take out relevant equilibria
% equilbria = equil.x(equil.x ~= 0)
% for i = 1:length(equil.x)
%     if equil.x(i) ~= 0
%         equilbria(i,1) = equil.x(i);
%     end
% end