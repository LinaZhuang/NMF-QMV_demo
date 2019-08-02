function [P,dist] = align_matrices(M,Mref,criterion) 
%
%   P = align_matrices(M,Mref) 
%   align two matrices per angle of SQUARE ERROR 
%
% ---  input variables ----
%
% Mref  - ref matrix
%
% M - matrix to be aligned  such that 
%    
% P - permutation matrix such that
%
%     P = arg min d(M*Q,Mref)
%            Q
%
% criterion - {'angle', 'MSE'}   
%              angle - mean absolute angle
%              MSE - mean square error
%
%       

% compute distance matrices
[m,n] = size(M);
if strcmp(criterion, 'angle')  %angle - mean absolute angle
    normM = repmat(sqrt(sum(M.^2)),m,1);
    normMref = repmat(sqrt(sum(Mref.^2)),m,1);
    d = (1-abs((M./normM)'*(Mref./normMref)))*100;  % use a [0 100] scale
else  % MSE - mean square error
    normM = repmat(sqrt(sum(M.^2))',1,n);
    normMref = repmat(sqrt(sum(Mref.^2))',1,n);
    d =  sqrt(normM.^2 + normMref'.^2 - 2*M'*Mref); 
end
%save for output
dist = d;

% find permutation
P = zeros(n);


for i =1:n
 % find  min
 [aux, index_row] = min(d);
 [dmin, index_col] = min(aux);
 P(index_row(index_col),index_col) = 1;
 % set column and row in inf
 d(:,index_col) = inf;
 d(index_row(index_col),:) = inf;
    
end

    
    






