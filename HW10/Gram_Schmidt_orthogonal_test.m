clc, clear, close all
%{
  Gram-Schmidt orthogonal decomposition test
%}

u1 = [1;5;9;10];
u2 = [5;8;2;8];
u3 = [9;8;1;3];
u4 = [9;8;1;3];
u5 = [0;-5;19;3];
% almost singular matrix
% singular matrix
X = [u1 u2 u3 u4 u5];
[row, col] = size(X);
Bn = CGS(X, true)
B = CGS(X, false)
V_modified = MGS(X)    % numberic stable
