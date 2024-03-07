% This script is a slightly modified version of code by Andrew Hines.
% This function calculates NSIM scores. NSIM is based on SSIM.
% The only changes I made to the script:
% 1. Create new variables specifying the window parameters
% 2. Return NSIM scores for windows with specific CFs as well as the mean NSIM score

function [freq_NSIMs, mNSIM] = nsim_calc(neuro_r, neuro_d)
    % Window parameters
    window_CFs = 3; % Must be an odd natural number
    window_bins = 3; % Must be an odd natural number
    gaussian_filter_std_dev = 0.5;

    window = fspecial('gaussian', [window_CFs window_bins], gaussian_filter_std_dev);
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
    freq_NSIMs = mean(L_r_d.*S_r_d, 2);
    mNSIM = mean(freq_NSIMs);
end
