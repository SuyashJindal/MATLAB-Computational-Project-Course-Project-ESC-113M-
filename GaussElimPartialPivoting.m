% Solve Ax = b using Gauss Elimination with Partial Pivoting
A = [1 1 1 ; 2 1 3;3 4 -2];
b = [4;7;9];
%% Gauss Elimination
% Get Augmented matrix
Ab = [A,b];
n = length(A);
% Row exchange to ensure A(1,1) is the largest column-1
col1 = Ab(:,1);
[dummy,idx] = max(col1);
dummy = Ab(1,:);
Ab(1,:) = Ab(idx,:);
Ab(idx,:) = dummy;
% Computation in the pivot column
for i =2:3
    alpha = Ab(i,1)/Ab(1,1);
    Ab(i,:) = Ab(i,1)/Ab(1,1);
    Ab(i,:) = Ab(i,:) - alpha*Ab(1,:);
end
col2 = Ab(2:end,2);
[dummy,idx] = max(col2);
dummy = Ab(2,:);
Ab(2,:) = Ab(idx,:);
Ab(idx,:) = dummy;
i = 3;
alpha = Ab(i,2)/Ab(2,2);
Ab(i,:) = Ab(i,:) - alpha * Ab(2,:);
%% Back substitutuon
x = zeros(3,1);
for i = 3:-1:1
    x(i) = (Ab(i,end)- (Ab(i,i+1:n)*x(i+1:n)))/Ab(i,i);
end


