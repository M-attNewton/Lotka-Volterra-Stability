
syms x y z q1 q2 q3 q4 q5 q6
Q = [q1, q4, q5; q4, q2, q6; q5, q6, q3];
 
vars = [x; y; z];
syms d a1 a2 a3 c
A = [d, a1, a2; a1, d, a3; a2, a3, d];
 
syms x_star y_star z_star
syms alpha1 alpha2 alpha3
 
vars_star = [x_star;y_star;z_star];
f = [c*y*z;c*x*z;c*x*y];
f_star = c*[y_star*z_star; x_star*z_star; x_star*y_star];
alpha = [alpha1, alpha2, alpha3];
V = -2*(vars - vars_star).'*Q*diag(vars)*(A*(vars - vars_star) + f - f_star) - alpha*(f - f_star);