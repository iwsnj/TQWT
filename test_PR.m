%% Tunable Q-factor Wavelet Transform: Verify PR
% Verify perfect reconstruction (PR) property of TQWT

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% https://eeweb.engineering.nyu.edu/iselesni/TQWT
% Ivan Selesnick
% selesi@nyu.edu
% NYU Tandon School of Engineering

% November 2010
% September 06, 2026
% - Updated to include radix-2 TQWT in same file

clear

for k = 1:2

    switch k
        case 1
            Q = 4; r = 3; J = 10;   % High Q-factor wavelet transform
        case 2
            Q = 1; r = 3; J = 5;    % Low Q-factor wavelet transform
    end

    fprintf('\n')
    fprintf('TQWT with Q = %3.2f, r = %3.2f\n', Q, r)

    beta = 2/(Q+1);
    alpha = 1-beta/r;

    for i = 1:2
        switch i
            case 1
                xform = @tqwt;
                inv_xform = @itqwt;
                fprintf('Real TQWT:\n')
            case 2
                xform = @tqwt_radix2;
                inv_xform = @itqwt_radix2;
                fprintf('Radix-2 Real TQWT:\n')
        end

        % Verify PR for various lengths
        for N = [400:2:420 2.^(7:10)]
            x = rand(1,N) + 1j*rand(1,N);    % Make test signal (complex-valued)
            J = floor(log2(beta*N/8)/log2(1/alpha));    % number of levels

            w = xform(x,Q,r,J);             % TQWT
            y = inv_xform(w,Q,r,N);         % Inverse TQWT
            recon_err = max(abs(x - y));    % Reconstruction error
            % print reconstruction error:
            fprintf('N = %4d, J = %3d: recon error = %e \n', N, J, recon_err)
        end

    end

end


