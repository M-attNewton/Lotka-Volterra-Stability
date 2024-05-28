% Work out stability conditions

% Declare variables
clear all
syms x y z u v w b d a1 a2 a3 c % where: u = ln(x), v = ln(y), w = ln(z)
vars = [x; y; z; u; v; w];

% Use descriminate to determine properties
P4 = 2*a1*a3*c^3 + a1^2*a3*c + c^4;
P3 = 6*a1*a2*a3*c^2 + a1^2*a2*a3 - 2*a1*a3*c + 4*a2*c^3 - a1^2*c^2*d ...
    -a1^3*d - a3*c^3 - a1*c^3 - a1^2*c;
P2 = -2*a1^2*a2*c*d + 6*a1*a2^2*a3*c + a1*c^2*d - 3*a2*a3*c^2 - 3*a1*a2*c^2 ...
    -2*a1*a2*a3 + a3*c + 2*a1*c + 6*a2^2*c^2 + 3*a1^2*d - a1^2*a2;
P1 = 2*a1*a2*c*d - 3*a2^2*a3*c - 3*a1*a2^2 + 2*a1*a3^3*a3 + 4*a3^3*c ...
    -a1^2*a2^2*d - 3*a1*d + a2*a3 + 2*a1*a2 - c;
P0 = a1*a2^2*d + d - a2^3*a3 - a1*a2^3 + a2^4 - a2;

delta = 256*P4^3*P0^3 - 192*P4^2*P3*P1*P0^2 - 128*P4^2*P2^2*P0^2 + 144*P4^2*P2*P1^2*P0 ...
    -27*P4^2*P0^4 + 144*P4*P3^2*P2*P0^2 - 6*P4*P3^2*P1^2*P0 - 80*P4*P3*P2^2*P1*P0 ...
    + 18*P4*P3*P2*P1^3 + 16*P4*P2^4*P0 - 4*P4*P2^3*P1^2 - 27*P3^4*P0^2 + 18*P3^3*P2*P1*P0 ...
    -4*P3^3*P1^3 - 4*P3^2+P2^3*P0 + P3^2*P2^2*P1^2;

Q1 = 8*P4*P2 - 3* P3^2;

Q2 = P3^3 + 8*P1*P4^2 - 4*P4*P3*P2;

Q3 = P2^2 - 3*P3*P1 + 12*P4*P0;

Q4 = 64*P4^3*P0 - 16*P4^2*P2^2 + 16*P4*P3^2*P2 - 16*P4^2*P3*P1 - 3*P3^4;



% Sub equilibria into state space to isoluate a equation for z
X = (1/d)*(1 - a2*y - a3*z - c*y*z);
Y = (1/d)*(1 - a2*x - a1*z - c*x*z);
Z = (1/d)*(1 - a3*x - a1*y - x*y*x);

Y = subs(Y,x,X);
Y = solve(Y,y);
Z = subs(Z,x,X);
Z = subs(Z,y,Y);
Z = solve(Z,z);

%equ1 = simplify(Z)

% Calculate jacobian of and sub to find z
x_dot = x*(1 - d*x - a1*y - a2*z - c*y*z);
y_dot = y*(1 - a1*x - d*y - a3*z - c*x*z);
z_dot = z*(1 - a2*x - a3*y - d*z - c*x*y);

jacob = jacobian([x_dot; y_dot; z_dot], [x; y; z]);

eigen = eig(jacob);

equ2 = subs(eigen,x,X);
equ2 = subs(equ2,y,Y);

%equ2 = simplify(equ2)

S = solve([equ1 == 0; real(equ2) > 0])


% % Define state space
% x_dot = x*(1 - d*x - a1*y - a2*z - c*y*z);
% y_dot = y*(1 - a1*x - d*y - a3*z - c*x*z);
% z_dot = z*(1 - a2*x - a3*y - d*z - c*x*y);
% 
% x_dot_f = (1 - d*x - a1*y - a2*z - c*y*z);
% y_dot_f = (1 - a1*x - d*y - a3*z - c*x*z);
% z_dot_f = (1 - a2*x - a3*y - d*z - c*x*y);
% 
% % Calculate equilibrium points
% jacob = jacobian([x_dot; y_dot; z_dot], [x; y; z]);
% eqns = [x_dot_f == 0; y_dot_f == 0; z_dot_f == 0; ...
%         real(eig(jacob)) > 0; x > 0; y > 0; z > 0];
% [equil1,equil2,equil3,params,conditions] = solve(eqns, [x y z],  'ReturnConditions', true);
% 


%equil = solve(eqns, [x;y;z], 'MaxDegree', 5); %x+(x==0)*eps; y+(y==0)*eps; z+(z==0)*eps

% % Test linear stability
% for i = 1:length(equil.x)
%     
%     % Repeat state vector field to avoid errors
%     x_dot = x*(1 - d*x - a1*y - a2*z + c*y*z);
%     y_dot = y*(1 - a1*x - d*y - a3*z + c*x*z);
%     z_dot = z*(1 - a2*x - a3*y - d*z + c*x*y);
%     
%     % Calculate Jacobian for equilibrium
%     jacob = jacobian([x_dot; y_dot; z_dot], [x; y; z]);
%     jacobSubs(:,:,i) = subs(jacob, [x; y; z], [equil.x(i); equil.y(i); equil.z(i)]);
%     
%     % Find eigenvalues of Jacobian
%     eigen(:,i) = eig(jacobSubs(:,:,i));
%     
%     % Check sign of eigenvalues to determine stability of equilibruim
%     if real(double(subs(eigen(1,i)))) >= 0 || ...
%        real(double(subs(eigen(2,i)))) >= 0 || ...
%        real(double(subs(eigen(3,i)))) >= 0
%         stability(i) = 0;
%     else
%         stability(i) = 1;
%     end
%     
% end
% 
% % Find relevant equilibria - assume only one for now
% ind = find(stability);
% sequil = [equil.x(ind); equil.y(ind); equil.z(ind)];

