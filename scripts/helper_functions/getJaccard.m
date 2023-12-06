function J = getJaccard(A,B)
% 
% call:
% 
%      J = getJaccard(A,B)
%      
% Compute the Jaccard Index, a measure of similarity between two binary (0,1) vector-sets A, B.

% INPUT
% 
%        A  :  Binary vector of set A
%        B  :  Binary vector of set B
%        
% OUTPUT
% 
%        J  :  Jaccard Index
%        
%        

% ----------------------------------------------------------------------------------------------------------------

if nargin<2
    error('The function requires two input-vectors')
end

if (sum(A==0) + sum(A==1)) < length(A)
    error('Input A must be binary (0,1) !!!')
end
  
if (sum(B==0) + sum(B==1)) < length(B)
    error('Input B must be binary (0,1) !!!')
end

J = sum(A.*B)/(sum(A+B) - sum(A.*B));


