
clear all
% Define variables
% A =[ -0.7027    0.6120
%      -0.8975    0.5415];
%  
% b = [0.3967;
%      0.6267];
%  
% X_star = [1;0.5];

n = 3;

X_star = [1;0.1;0.1];

A = [-1      1    2
     -0.1    0.2   -0.4
     -0.1    0.3   -0.3];
 
b = [0.7;0.12;0.1]; 
b = -A*X_star;

% Vary alpha until valid one is found
for alpha1 = -1:0.05:1
for alpha2 = -1:0.05:1
for alpha3 = -1:0.05:1

clear Q eta1    
    
cvx_begin

variable Q(n,n,n) 
variable eta1 nonnegative

alph = [alpha1, 0, 0; 0, alpha2, 0; 0, 0, alpha3];

N = diag(X_star)*Q(:,:,1) + diag(X_star)*Q(:,:,2) + diag(X_star)*Q(:,:,3);

M = alph*(-A + N);
M + eta1*eye(n) <= 0 

eta1 == 0.001

for i = 1:n
    for j = 1:n
        for k = 1:n
            
            if i+j+k ~= 6 || i == j || i == k || j == k
                Q(i,j,k) == 0;
            end
            
        end
    end
end

% Q(:,:,1) == Q(:,:,1).';
% Q(:,:,2) == Q(:,:,2).';
% Q(:,:,3) == Q(:,:,3).';
% alpha1*Q(1,2,3) + alpha2*Q(2,1,3) + alpha3*Q(3,1,2) == 0

alpha1*Q(1,2,3) + alpha1*Q(1,3,2) + alpha2*Q(2,1,3) + alpha2*Q(2,3,1) + ...
alpha3*Q(3,1,2) + alpha3*Q(3,2,1) == 0;

cvx_end
if cvx_optval ~= inf || isnan(cvx_optval)
    return
end

end        
end
end



