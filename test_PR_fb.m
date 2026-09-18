
%% test_PR_fb.m
%
%   Test perfect reconstruction property of afb and sfb.
%   The program afb followed by sfb should reconstruct
%   the input signal of afb.
%
% September 15, 2026

N = 16;
N0 = 10;
N1 = 12;

X = complex(rand(1, N), rand(1, N));  % Complex 1D array
[V0, V1] = afb(X, N0, N1);
Y = sfb(V0, V1, N);
max(abs(X - Y))                 % Should be zero

