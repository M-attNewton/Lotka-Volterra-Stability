%% 3D vector field with cubic terms 

%% Set up variables
clear all
syms x y z u v w % where: u = ln(x), v = ln(y), w = ln(z)
vars = [x; y; z; u; v; w];

% Define parameters
d = 3; % Diagonal elements
a1 = 2; a2 = 1; a3 = 1; % [xy; xz; yz] terms
c = 3; % Cubic terms
b = 1; % Normalise system so birth rate is 1

% True vector field
x_dot = x*(1 - d*x - a1*y - a2*z - c*y*z);
y_dot = y*(1 - a1*x - d*y - a3*z - c*x*z);
z_dot = z*(1 - a2*x - a3*y - d*z - c*x*y);

% Remove zero terms to speed up simultanious equation solver
x_dot_f = (1 - d*x - a1*y - a2*z - c*y*z);
y_dot_f = (1 - a1*x - d*y - a3*z - c*x*z);
z_dot_f = (1 - a2*x - a3*y - d*z - c*x*y);

% Calculate equilibrium points
jacob = jacobian([x_dot; y_dot; z_dot], [x; y; z]);
%eqns = [x_dot_f == 0; y_dot_f == 0; z_dot_f == 0];

eqns = [x_dot_f == 0; y_dot_f == 0; z_dot_f == 0]; %; ...
        %x > 0; y > 0; z > 0 ]; % real(eig(jacob)) > 0
%equil = solve(eqns, [x y z],  'ReturnConditions', true);
equil = solve(eqns, [x;y;z]); %x+(x==0)*eps; y+(y==0)*eps; z+(z==0)*eps

% Test linear stability
stability = zeros(size(equil.x));
instability = stability;
for i = 1:length(equil.x)
    
    % Repeat state vector field to avoid errors
    x_dot = x*(1 - d*x - a1*y - a2*z - c*y*z);
    y_dot = y*(1 - a1*x - d*y - a3*z - c*x*z);
    z_dot = z*(1 - a2*x - a3*y - d*z - c*x*y);
    
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
        if vpa(equil.x(i)) > 0 && vpa(equil.y(i)) > 0 && vpa(equil.z(i)) > 0
            stability(i) = 1;
        else
            stability(i) = 0;
        end
    end
    
    if stability(i) == 0 && vpa(equil.x(i)) > 0 && vpa(equil.y(i)) > 0 && vpa(equil.z(i)) > 0
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

sequil = vpa([equil.x(ind); equil.y(ind); equil.z(ind)]);

ind2 = find(instability);
isequil = vpa([equil.x(ind2); equil.y(ind2); equil.z(ind2)]);

x_dot = x*(1 - d*x - a1*y - a2*z - c*y*z);
y_dot = y*(1 - a1*x - d*y - a3*z - c*x*z);
z_dot = z*(1 - a2*x - a3*y - d*z - c*x*y);

f = [x_dot; y_dot; z_dot; ...
     (x^-1)*x_dot; (y^-1)*y_dot; (z^-1)*z_dot];
 
% Initialize the sum of squares program
prog = sosprogram(vars);

% The Lyapunov function V(x): 
Z1 = monomials(vars,0:2); 
%Z1 = [x*y; x*z; y*z; x; y; z; u; v; w; 1];
Z1 = [x^2; y^2; z^2; x*y; x*z; y*z; x; y; z; u; v; w; 1];
[prog,V] = sospolyvar(prog,Z1,'wscoeff');

Z2 = monomials(vars,0:2);
%Z2 = monomials([x; y; z],0:1);
[prog,p1] = sospolyvar(prog,Z2,'wscoeff');
%[prog,p1] = sossosvar(prog,Z2,'wscoeff');

%[prog,p2] = sospolyvar(prog,Z2,'wscoeff')
%[prog,p3] = sospolyvar(prog,Z2,'wscoeff')
%[prog,Q] = sospolymatrixvar(prog,monomials([x;y;z],0),[3 3],'symmetric')

%% Define SOSP constraints

% Constraint 1: -dV/dx*f + D*p >= 0
% Region D
D = ((x-sequil(1))^2)/(0.5*sequil(1))^2 + ((y-sequil(2))^2)/(0.5*sequil(2))^2 ...
    + ((z-sequil(3))^2)/(0.5*sequil(3))^2 - 1;
D1 = -x ; D2 = -y; D3 = -z;

expr2 = -(diff(V,x)*f(1) + diff(V,y)*f(2) + diff(V,z)*f(3) + ...
          diff(V,u)*f(4) + diff(V,v)*f(5) + diff(V,w)*f(6))+ p1*D; %...
         % + p1*D1 + p2*D2 + p3*D3;
prog = sosineq(prog,expr2);

% Constraint 2: p >= 0
prog = sosineq(prog,p1);
%prog = sosineq(prog,p2);
%prog = sosineq(prog,p3);

% Constraint 3: d^2V/dx^2 >= 0 at x = star
Q = [coeff_1, 0.5*coeff_4, 0.5*coeff_5; ...
     0.5*coeff_4, coeff_2, 0.5*coeff_6; ...
     0.5*coeff_5, 0.5*coeff_6, coeff_3];

matexpr1 = [(x^2*coeff_1 - coeff_10), x^2*coeff_4, x^2*coeff_5;
            y^2*coeff_4, (y^2*coeff_2 - coeff_11), y^2*coeff_6;
            z^2*coeff_5, z^2*coeff_6, (z^2*coeff_3 - coeff_12)];
prog = sosmatrixineq(prog,matexpr1,'quadraticMineq');

% matexpr1 = 2*Q + [(-coeff_10/(sequil(1))^2), 0, 0;...
%                 0, (-coeff_11/(sequil(2))^2), 0;...
%                 0, 0, (-coeff_12/(sequil(3))^2)];

% matexpr1 = 2*Q + [(-coeff_10/(x)^2), 0, 0;...
%                 0, (-coeff_11/(y)^2), 0;...
%                 0, 0, (-coeff_12/(z)^2)];

% Constraint 4: dV/dx == 0 at x == xstar
expr4 = coeff_7 + coeff_10/sequil(1) + 2*coeff_1*sequil(1) + coeff_4*sequil(2) ...
    + coeff_5*sequil(3);
prog = soseq(prog,expr4);

expr5 = coeff_8 + coeff_11/sequil(2) + 2*coeff_2*sequil(2) + coeff_4*sequil(1) ...
    + coeff_6*sequil(3);
prog = soseq(prog,expr5);

expr6 = coeff_9 + coeff_12/sequil(3) + 2*coeff_3*sequil(3) + coeff_5*sequil(1) ...
    + coeff_6*sequil(2);
prog = soseq(prog,expr6);

% Constraint 5: Fix parameters normalise Lyapunov function
% At equilbrium Lyaponuv function must be zero
expr3 = coeff_1*sequil(1)^2 + coeff_4*sequil(1)*sequil(2) + coeff_5*sequil(1)*sequil(3) + ...
    coeff_7*sequil(1) + coeff_2*sequil(2)^2 + coeff_6*sequil(2)*sequil(3) + ...
    coeff_8*sequil(2) + coeff_3*sequil(3)^2 + coeff_9*sequil(3) + coeff_13 + ...
    coeff_10*log(sequil(1)) + coeff_11*log(sequil(2)) + coeff_12*log(sequil(3));
prog = soseq(prog,expr3);

% x, y and z coefficients must be 1
prog = soseq(prog,coeff_7 - 1);
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
%fsurf(SOLVXY)

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
