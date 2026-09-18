function [V0, V1] = afb(X, N0, N1)
% [V0, V1] = afb(X, N0, N1)
% afb: analysis filter bank
% Converts a vector X into two vectors:
% V0 of length N0, V1 of length N1
%
% Need:
%   N0, N1, both even integers
%   2 <= N0 <= length(X)
%   2 <= N1 <= length(X)
%   N0 + N1 > length(X)
%
% % Example (verify perfect reconstruction)
% N = 20;
% X = rand(1, N);
% [V0, V1] = afb(X, 16, 12);
% Y = sfb(V0, V1, N);
% max(abs(X - Y))

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% https://eeweb.engineering.nyu.edu/iselesni/TQWT
% Ivan Selesnick
% selesi@nyu.edu
% NYU Tandon School of Engineering

% November 2010
% Revised: September 06, 2026
% - Simplify code

X = X(:).';                                     % X is row vector
N = length(X);

P = (N-N1)/2;
T = (N0+N1-N)/2 - 1;
S = (N-N0)/2;

% transition-band function
v = (1:T)/(T+1)*pi;
c = cos(v);
trans = (1+c) .* sqrt(2-c)/2;
% i.e., trans = (1+cos(v)) .* sqrt(2-cos(v))/2;

trans_rev = trans(end:-1:1);

% Add 1 to indices because Matlab indexing starts at 1 (not 0)

XT_pos = X(P+2:P+T+1);      % positive-freq transition-band samples
XT_neg = X(N-P-T+1:N-P);    % negative-freq transition-band samples

% low-pass subband
V0 = zeros(1,N0);
V0(1:P+1) = X(1:P+1);                   % pass-band (dc and pos freq)
V0(P+2:P+T+1) = XT_pos.*trans;          % trans-band (pos freq)
V0(N0-P-T+1:N0-P) = XT_neg.*trans_rev;  % trans-band (neg freq)
V0(N0-P+1:N0) = X(N-P+1:N);             % pass-band (neg freq)

% high-pass subband
V1 = zeros(1,N1);
V1(2:T+1) = XT_pos.*trans_rev;      % trans-band (pos freq)
V1(T+2:N1-T) =  X(P+T+2:N-P-T) ;    % pass-pand (pos and neg freq)
V1(N1-T+1:N1) = XT_neg.*trans;      % trans-band (neg freq)

