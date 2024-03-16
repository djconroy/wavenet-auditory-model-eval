% This script is a modified version of code by Andrew Hines.
% This function calculates NSIM scores.
% NSIM is based on SSIM. The code for SSIM is available here:
% https://www.cns.nyu.edu/~lcv/ssim/

% neuro_r Reference neurogram to compare against
% neuro_d Degraded neurogram 

function [freq_NSIMs, mean_NSIM] = nsim(neuro_r, neuro_d)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% START OF CODE WRITTEN BY ME %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Approach to implementing asymmetric Gaussian filter
    % taken from https://stackoverflow.com/a/15241151
    vertical = fspecial('gaussian', [5 1], 1);
    horizontal = fspecial('gaussian', [1 3], 0.5);
    window = vertical * horizontal;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% START OF CODE NOT WRITTEN BY ME %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %window = fspecial('gaussian', [3 3], 0.5);
    window = window/sum(sum(window));
    %dynamic range set to max of reference neuro
    L=max(max(neuro_r));
    %C1 and C2 constants
    K=[0.01 0.03];
    C1 = (K(1)*L)^2;
    C2 = ((K(2)*L)^2)/2;
    %Calc mean NSIM(r,d)
    neuro_r = double(neuro_r);
    neuro_d = double(neuro_d);
    mu_r   = filter2(window, neuro_r, 'valid');
    mu_d   = filter2(window, neuro_d, 'valid');
    mu_r_sq = mu_r.*mu_r;
    mu_d_sq = mu_d.*mu_d;
    mu_r_mu_d = mu_r.*mu_d;
    sigma_r_sq = filter2(window, neuro_r.*neuro_r, 'valid') - mu_r_sq;
    sigma_d_sq = filter2(window, neuro_d.*neuro_d, 'valid') - mu_d_sq;
    sigma_r_d = filter2(window, neuro_r.*neuro_d, 'valid') - mu_r_mu_d;
    sigma_r=sign(sigma_r_sq).*sqrt(abs(sigma_r_sq));
    sigma_d=sign(sigma_d_sq).*sqrt(abs(sigma_d_sq));
    L_r_d= (2*mu_r.*mu_d+C1) ./(mu_r_sq+mu_d_sq +C1);
    S_r_d= (sigma_r_d + C2)./(sigma_r.*sigma_d +C2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% START OF CODE WRITTEN BY ME %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Return NSIM scores for windows with specific CFs as well as the mean NSIM score
    freq_NSIMs = mean(L_r_d.*S_r_d, 2);
    mean_NSIM = mean(freq_NSIMs);
end
