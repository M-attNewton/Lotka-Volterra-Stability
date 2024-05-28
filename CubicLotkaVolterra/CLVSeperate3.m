%% 3D vector field with cubic terms, two Lyapunov functions solved seperately
% V = V_nonpoly + V_poly
%% Set up variables
clear all
syms x y z 
vars = [x; y; z;];

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

x_dot_f = (b(1) - A(1,1)*(x+sequil(1)) - A(1,2)*(y+sequil(2)) ...
    - A(1,3)*(z+sequil(3)) + c(1)*(y+sequil(2))*(z+sequil(3)));
y_dot_f = (b(2) - A(2,1)*(x+sequil(1)) - A(2,2)*(y+sequil(2)) ...
    - A(2,3)*(z+sequil(3)) + c(2)*(x+sequil(1))*(z+sequil(3)));
z_dot_f = (b(3) - A(3,1)*(x+sequil(1)) - A(3,2)*(y+sequil(2)) ...
    - A(3,3)*(z+sequil(3)) + c(3)*(x+sequil(1))*(y+sequil(2)));

f = [x_dot; y_dot; z_dot]; 
 
% Initialize the sum of squares program
prog = sosprogram(vars);

% The Lyapunov function V_poly(x): 
Z1 = [x^2; y^2; z^2; x*y; x*z; y*z];
Z2 = monomials(vars,0:1);
[prog,V_poly] = sospolyvar(prog,Z1,'wscoeff');
[prog,p] = sossosvar(prog,Z2,'wscoeff');

% Parameters for Lyapunov function V_nonpoly(x):
[prog,c1] = sossosvar(prog,sym(1));
[prog,c2] = sossosvar(prog,sym(1));
[prog,c3] = sossosvar(prog,sym(1));
prog = sosineq(prog,c1-0.1);
prog = sosineq(prog,c2-0.1);
prog = sosineq(prog,c3-0.1);

%% Define SOSP constraints
% Constraint 1: -dV/dx*f + D*p >= 0
% Region D
D = x^2 + y^2 + z^2 - 0.2^2;

% V_nonpoly defined from previous result
V_nonpoly = c1*(x - sequil(1)*log(x+sequil(1)))+c2*(y-sequil(2)*log(y+sequil(2)))+...
    c3*(z-sequil(3)*log(z+sequil(3)));
%V_nonpoly = c1*x - sequil(1)*log(x+sequil(1)) + c2*y-sequil(2)*log(y+sequil(2))+...
 %   c3*z-sequil(3)*log(z+sequil(3));
%V_nonpoly = c1*(x - sequil(1)*log(x)) + c2*(y-sequil(2)*log(y)) + c3*(z-sequil(3)*log(z));

% V_full is the total Lyapunov function
V_full = V_nonpoly + V_poly;

expr2 = -(diff(V_poly,x)*f(1) + diff(V_poly,y)*f(2) + diff(V_poly,z)*f(3)) + p*D + ...
    - c1*(1-sequil(1))*x_dot_f - c2*(1-sequil(2))*y_dot_f - c3*(1-sequil(3))*z_dot_f;

prog = sosineq(prog,expr2);

% Constraint 2: V_poly >= 0
prog = sosineq(prog,V_poly);

% % Constraint 3: V_full(0) == 0 
% expr3 = c1*(log(sequil(1))*(-sequil(1))) + c2*(log(sequil(2))*(-sequil(2))) +...
% c3*(log(sequil(3))*(-sequil(3)));
% 
% prog = soseq(prog,expr3);

%% Call solver
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);

%% Get solution
SOLV = sosgetsol(prog,V_full,20)

% Shift graph to equilibrium
SOLV2 = subs(SOLV,[x,y,z],[x-sequil(1),y-sequil(2),z-sequil(3)]);

SOLV3 = vpa(SOLV2,4)

%SOLV3 = SOLV

%% Plot Solution
% Calculate cross sections to plot 
SOLVXY = subs(SOLV3,z,sequil(3));
SOLVXZ = subs(SOLV3,y,sequil(2));
SOLVYZ = subs(SOLV3,x,sequil(1));

% Define axis
axis1 = -1; axis2 = 1;
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
