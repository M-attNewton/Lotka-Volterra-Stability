%% 3D vector field with cubic terms 

%% Set up variables
clear all; close all;
syms x y z % where: u = ln(x), v = ln(y), w = ln(z)
vars = [x; y; z];

% Define parameters
X_star = [1;0.1;0.1];

A = [-1      1    2
     -0.1    0.2   -0.4
     -0.1    0.3   -0.3];
 
b = [0.7;0.12;0.1]; 
b = -A*X_star;

k = 1;
for c = -2:0.1:2

% True vector field
x_dot = x*(b(1) + A(1,1)*x + A(1,2)*y + A(1,3)*z - c*y*z);
y_dot = y*(b(2) + A(2,1)*x + A(2,2)*y + A(2,3)*z - c*x*z);
z_dot = z*(b(3) + A(3,1)*x + A(3,2)*y + A(3,3)*z - c*x*y);

% Remove zero terms to speed up simultanious equation solver
x_dot_f = (b(1) + A(1,1)*x + A(1,2)*y + A(1,3)*z - c*y*z);
y_dot_f = (b(2) + A(2,1)*x + A(2,2)*y + A(2,3)*z - c*x*z);
z_dot_f = (b(3) + A(3,1)*x + A(3,2)*y + A(3,3)*z - c*x*y);

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
    x_dot = x*(b(1) + A(1,1)*x + A(1,2)*y + A(1,3)*z - c*y*z);
    y_dot = y*(b(2) + A(2,1)*x + A(2,2)*y + A(2,3)*z - c*x*z);
    z_dot = z*(b(3) + A(3,1)*x + A(3,2)*y + A(3,3)*z - c*x*y);
    
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
    
    if stability(i) == 0 && vpa(equil.x(i)) > 0 && vpa(equil.y(i)) > 0 ...
            && vpa(equil.z(i)) > 0
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

x_dot = x*(b(1) + A(1,1)*x + A(1,2)*y + A(1,3)*z - c*y*z);
y_dot = y*(b(2) + A(2,1)*x + A(2,2)*y + A(2,3)*z - c*x*z);
z_dot = z*(b(3) + A(3,1)*x + A(3,2)*y + A(3,3)*z - c*x*y);

x_dot_s = (b(1) + A(1,1)*x + A(1,2)*y + A(1,3)*z - c*y*z);
y_dot_s = (b(2) + A(2,1)*x + A(2,2)*y + A(2,3)*z - c*x*z);
z_dot_s = (b(3) + A(3,1)*x + A(3,2)*y + A(3,3)*z - c*x*y);


f = [subs(x_dot,[x,y,z],[x+sequil(1),y+sequil(2),z+sequil(3)]);
    subs(y_dot,[x,y,z],[x+sequil(1),y+sequil(2),z+sequil(3)]);
    subs(z_dot,[x,y,z],[x+sequil(1),y+sequil(2),z+sequil(3)])];

f_s = [subs(x_dot_s,[x,y,z],[x+sequil(1),y+sequil(2),z+sequil(3)]);
    subs(y_dot_s,[x,y,z],[x+sequil(1),y+sequil(2),z+sequil(3)]);
    subs(z_dot_s,[x,y,z],[x+sequil(1),y+sequil(2),z+sequil(3)])];
 
% Initialize the sum of squares program
prog = sosprogram(vars);

% The Lyapunov function V(x): 
[prog,coeff_1] = sospolyvar(prog,sym(1));
[prog,coeff_2] = sospolyvar(prog,sym(1));
[prog,coeff_3] = sospolyvar(prog,sym(1));
coeff_4 = coeff_1;
coeff_5 = coeff_2;
coeff_6 = coeff_3;
prog = sosineq(prog,coeff_1-0.1);
%prog = sosineq(prog,coeff_2-0.1);
%prog = sosineq(prog,coffe_3-0.1);
%[prog,coeff_4] = sospolyvar(prog,sym(1));
%[prog,coeff_5] = sospolyvar(prog,sym(1));
%[prog,coeff_6] = sospolyvar(prog,sym(1));

[prog,V_quad] = sossosvar(prog,vars);
%V_quad = 0;
V = V_quad+coeff_1*x+coeff_2*y+coeff_3*z-coeff_4*sequil(1)*log(x+sequil(1))...
    -coeff_5*sequil(2)*log(y+sequil(2))-coeff_6*sequil(3)*log(z+sequil(3));
Q = subs(jacobian(jacobian(V,[x,y,z]),[x,y,z]),[x,y,z],[sequil(1),sequil(2),sequil(3)]);
prog = sosineq(prog,[x,y,z]*Q*[x;y;z]-0.1*x^2-0.1*y^2-0.1*z^2);

Z2 = monomials(vars,0:2);
[prog,p1] = sospolyvar(prog,Z2,'wscoeff');
%[prog,p1] = sossosvar(prog,Z2,'wscoeff');

[prog,p2] = sospolyvar(prog,Z2,'wscoeff')
[prog,p3] = sospolyvar(prog,Z2,'wscoeff')

%% Define SOSP constraints

% Constraint 1: -dV/dx*f + D*p >= 0
% Region D
%D = x^2+y^2+z^2-sequil(1)^2;
D1 = -x-sequil(1) ; D2 = -y-sequil(2); D3 = -z-sequil(3);

expr2 =  -(diff(V_quad,x)*f(1) + diff(V_quad,y)*f(2) + diff(V_quad,z)*f(3))...
    -(coeff_1*x+coeff_1*sequil(1)-coeff_4*sequil(1))*f_s(1)-(coeff_2*y+coeff_2*sequil(2)...
    -coeff_5*sequil(2))*f_s(2) -(coeff_3*z+coeff_3*sequil(3)-coeff_6*sequil(3))*f_s(3)...
    +p1*D1+p2*D2+p3*D3;

prog = sosineq(prog,expr2);

% Constraint 2: p >= 0
prog = sosineq(prog,p1);
prog = sosineq(prog,p2);
prog = sosineq(prog,p3);

%% Call solver
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);

%% Get solution
SOLV = sosgetsol(prog,V)
K = subs(SOLV,[x,y,z],[0,0,0]);
SOLV = SOLV-K;
% Subsitute in log terms and reduce significant figures
SOLV3(k) = vpa(SOLV,4)
k = k + 1;

end

%% Plot Solution
% Calculate cross sections to plot 
SOLVXY = subs(SOLV3,z,0);
SOLVXZ = subs(SOLV3,y,0);
SOLVYZ = subs(SOLV3,x,0);

% Define axis
axis1 = double(-sequil(1)+0.1); axis2 = 1;
figure
%fsurf(SOLVXY)
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
