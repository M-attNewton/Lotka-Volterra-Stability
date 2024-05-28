
clear all
n = 3;
k = 0;
m = 0;

X_star = zeros(n,1);
for j = 1:n
    X_star(j) = (0.5)^(j-1);
end

for i = 1:10000

%clearvars -except i k test1 test2 test_mat1 test_mat2
clear D eta1 eta2 A eigen

% Create matrix A

A = randn(n);% - eye(n)*n;
% A = [-0.6825    0.6237;
%    -0.7868    0.2040];

% A = [    -1.1499   -0.5502   -1.0848;
%     -1.0589   -1.7808    1.8048;
%     -0.3987   -1.6770   -0.9149];

% Test D-stability DA + A'D < 0
 A = [    -1.3853    0.1509    0.3106;
    0.0909   -2.3554   -0.0908;
   -0.2140    0.6201   -0.3969];

cvx_begin

variable D(n,n) diagonal
variable eta1 nonnegative
variable eta2 nonnegative

%minimize(eta1 + eta2)

subject to
D*A + A'*D + eta1*eye(n) <= 0
D - eta2*eye(n) >= 0
% D*A + A'*D <= 0
% D >= 0

eta1 == 10^-2
eta2 == 10^-2

% eta1 == 10^-10
% eta2 == 10^-10

cvx_end

% Test eigenvalue stability
%X_star = [1;0.5;0.25];%\A;
eigen = eig(diag(X_star)*A); 
b = -A*X_star;

if cvx_optval == Inf && all(real(eigen) < 0) && all(b > 0) %% && all(X_star > 0)
    %b = -A*X_star;
    disp('Counterexample found')
    return
    k = k + 1;
    A_Dunstable(:,:,k) = A;
end

if cvx_optval ~= Inf && all(real(eigen) < 0) && all(b > 0) %% && all(X_star > 0)
    %b = -A*X_star;
    disp('Counterexample found')
    return
    m = m + 1;
    A_Dstable(:,:,m) = A;
end

end
