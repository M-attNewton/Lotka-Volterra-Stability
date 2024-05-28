
clear all
n = 3;
k = 0;
m = 0;

X_star = zeros(n,1);
for j = 1:n
    X_star(j) = (0.5)^(j-1);
end

for i = 1:500

%clearvars -except i k test1 test2 test_mat1 test_mat2
clear D eta1 eta2 A eigen

% Create matrix A
A = randn(n);

% A = [-0.6203   -1.6592   -0.0660;
%      1.7346    0.1874    0.4172;
%     -1.1939    2.4415   -0.4887];

% Test D-stability DA + A'D < 0
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

eta1 == 10^-3
eta2 == 10^-3

% eta1 == 10^-10
% eta2 == 10^-10

cvx_end

% Test eigenvalue stability
%X_star = [1;0.5;0.25];%\A;
eigen = eig(A); 
%b = -A*X_star;

if cvx_optval == Inf && all(real(eigen) < 0)
    B = A;
    return
end

if cvx_optval == Inf && all(real(eigen) < 0) %&& all(b > 0) %% && all(X_star > 0)
    %b = -A*X_star;
    disp('Counterexample found')
    %return
    k = k + 1;
    A_Dunstable(:,:,k) = A;
end

if cvx_optval ~= Inf && all(real(eigen) < 0) %&& all(b > 0) %% && all(X_star > 0)
    %b = -A*X_star;
    disp('Counterexample found')
    %return
    m = m + 1;
    A_Dstable(:,:,m) = A;
end

end
