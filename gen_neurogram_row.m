function [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% START OF CODE WRITTEN BY ANDREW HINES %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    psthbinwidth_env = 100e-6;
    env_windowsize = 128;
    env_overlap = 64;
    w_env = hamming(env_windowsize, 'periodic');

    psthbinwidth_tfs = 10e-6;
    tfs_windowsize = 32;
    tfs_overlap = 16;
    w_tfs = hamming(tfs_windowsize, 'periodic');

    Fs = 100e3;
    psthbins_env = round(psthbinwidth_env * Fs);

    psth_env = psth_tfs(1:psthbins_env * floor(length(psth_tfs) / psthbins_env));
    pr_env = sum(reshape(psth_env, psthbins_env, length(psth_env) / psthbins_env)); % pr of spike in each bin

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% START OF CODE WRITTEN BY ME %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if length(pr_env) < env_windowsize
        env_neurogram_row = 0;
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% START OF CODE WRITTEN BY ANDREW HINES %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        psth_env = pr_env / psthbinwidth_env; % psth in units of spikes/s
        b_sp_env = spectrogram(psth_env, w_env, env_overlap, env_windowsize, 1 / psthbinwidth_env);
        synrate_env = abs(b_sp_env / env_windowsize / sqrt(sum(w_env.^2) / length(w_env)));
        env_neurogram_row = synrate_env(1, :); % get the dc value
    end

    psth_tfs = psth_tfs / psthbinwidth_tfs; % psth in units of spikes/s
    b_sp_tfs = spectrogram(psth_tfs, w_tfs, tfs_overlap, tfs_windowsize, 1 / psthbinwidth_tfs);
    synrate_tfs = abs(b_sp_tfs / tfs_windowsize / sqrt(sum(w_tfs.^2) / length(w_tfs)));
    tfs_neurogram_row = synrate_tfs(1, :); % get the dc value
end
