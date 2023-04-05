function [i,j,k] = val2indexes(ix, a, b, c)
% val2indexes provides the indexes i,j,k from index ix
%   OUTPUT: 
%       * i - index of the first element
%       * j - index of the second element
%       * k - index of the third element
%   INPUT: 
%       * ix - index of the (i,j,k) combination
%       * a - size of elements containing "i"
%       * b - size of elements containing "j"
%       * c - size of elements containing "k"
%a=4;
%b=1;
%c=4;
%for ix=1:16
    i = mod(ix, a); 
    if i == 0, i = a; end    
    y = mod(ix, (a * b));
    j = ceil(y/a);
    if j == 0, j = b; end 
    k = ceil(ix / (a * b)); 
    if k > c, k = c; end
    %fprintf('%d %d %d\n',i,j,k);
end