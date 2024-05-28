
clear all
n = 3;
k = 0;
test_mat1 = zeros(n);
test_mat2 = test_mat1;

for i = 1:1000

%clearvars -except i k test1 test2 test_mat1 test_mat2
clear D eta1 eta2 A eigen

% Create matrix A

%A = randn(n);% - eye(n)*n;
%A = A-A'-n*eye(n);

%A = [-0.84    0.82    1.61
%   -0.02    0.22   -0.41
%   -0.03    0.31   -0.29];

A = [-1      1    2
     -0.1    0.2   -0.4
     -0.1    0.3   -0.3];

% Test D-stability DA + A'D < 0
cvx_begin

variable D(n,n)
variable eta1 nonnegative
variable eta2 nonnegative

subject to
D*A + A'*D - eta1*eye(n) <= 0
D - eta2*eye(n) >= 0

eta1 == 0.00001
eta2 == 0.00001

cvx_end

% Test eigenvalue stability
X_star = [1;0.1;0.1];%\A;
eigen = eig(diag(X_star)*A); 

if cvx_optval == Inf && all(real(eigen) < 0)%% && all(X_star > 0)
    b = -A*X_star;
    disp('Counterexample found')
    return
    k = k + 1;
    test1(k) = log(det(A));
    test_mat1 = test_mat1 + A;
else
    test2(i) = log(det(A));
    test_mat2 = test_mat2 + A;
end

end
