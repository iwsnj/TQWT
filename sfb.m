function Y = sfb(V0, V1, N)
% Y = sfb(V0, V1, N)
% sfb: synthesis filter bank

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% https://eeweb.engineering.nyu.edu/iselesni/TQWT
% Ivan Selesnick
% selesi@nyu.edu
% NYU Tandon School of Engineering

% November, 2010
% Revised: September 06, 2026
% - Simplify code
% - revised to remove unnecessary zero assignments..
% - revised indexing from a+(b:c) to (a+b):(a+c)

N0 = length(V0);
N1 = length(V1);

S = (N-N0)/2;
P = (N-N1)/2;
T = (N0+N1-N)/2 - 1;

% transition-band function
v = (1:T)/(T+1)*pi;
c = cos(v);
trans = 0.5 * (1+c) .* sqrt(2-c);
% i.e., trans = (1+cos(v)) .* sqrt(2-cos(v))/2;

trans_rev = trans(end:-1:1);

% Add 1 to indices because Matlab indexing starts at 1 (not 0)

% Y0 and Y1 only overlap in the two transition bands.

Y = zeros(1,N);

% low frequencies (dc and positive frequencies)
Y(1:P+1) = V0(1:P+1);                       

% transition band (positive frequencies)
Y(P+2:P+T+1) = V0(P+2:P+T+1).*trans + V1(2:T+1).*trans_rev;                 

% high frequencies
Y(P+T+2:N-P-T) = V1(T+2:N1-T);
% Y(N/2-S+1:N/2+S+1) = V1(T+2:N1-T);
% Y(N/2-S+1:N/2+S+1) = V1(N1/2-S+1:N1/2+S+1);

% transition band (negative frequencies)
Y(N-P-T+1:N-P) = V0(N0-P-T+1:N0-P).*trans_rev + V1(N1-T+1:N1).*trans;

% low frequencies (negative frequencies)
Y(N-P+1:N) = V0(N0-P+1:N0);
