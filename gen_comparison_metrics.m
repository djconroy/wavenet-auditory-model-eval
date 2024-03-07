function gen_comparison_metrics(audio_file_name, spl, snr)
    arguments
        audio_file_name (1, 1)
        spl (1, 1) = 60
        snr (1, 1) = "N/A"
    end

    wavenet_Fs = 16000; % Sampling rate of audio signal inputs the WaveNet model is fitted to
    wavenet_RF = 2048; % Size of the receptive field of the WaveNet model

    addpath(fullfile('models', 'original'))

    [audio_signal, audio_sampling_rate] = audioread(audio_file_name);

    audio_signal = amplify_audio_signal(audio_signal, spl);

    if isStringScalar(snr) && strcmp(snr, "N/A")
        noisy_audio_file = "";

        audio_name = audio_file_name;
    else
        % Add white Gaussian noise to the speech signal in accordance with
        % the specified signal-to-noise ratio (SNR)
        audio_signal = awgn(audio_signal, snr, 'measured');

        % Save the noisified speech signal to a temporary WAVE file for
        % use by the WaveNet model
        noisy_audio_file = "noisy.wav";
        audiowrite(noisy_audio_file, audio_signal, audio_sampling_rate)

        audio_name = noisy_audio_file;
    end

    outputDir = "outputs" + extractBetween(audio_file_name, pwd, '.WAV');
    mkdir(outputDir)

    % Run the WaveNet model
    system("python run_wavenet_model.py --file " + audio_name + " --spl " + num2str(spl));
    disp('Finished running WaveNet model')
    ihcogram_file = "ihcogram.mat";
    % Load the output produced by the WaveNet model, which includes the 'ihcogram' variable
    load(ihcogram_file)

    if noisy_audio_file == "noisy.wav"
        % Delete the temporary WAVE file of the noisified speech signal after
        % it's used by the WaveNet model
        delete(noisy_audio_file)
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% START OF CODE NOT WRITTEN BY ME %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    species = 2; % human
    numCFs = 80;
    CFs = logspace(log10(125), log10(8000), numCFs);

    cohcs  = ones(1, numCFs);  % normal ohc function
    cihcs  = ones(1, numCFs);  % normal ihc function

    numsponts = [10 10 30]; % Number of low-spont, medium-spont, and high-spont fibers at each CF in a healthy AN

    if exist('ANpopulation.mat','file')
        load('ANpopulation.mat');
        disp('Loading existing population of AN fibers saved in ANpopulation.mat')
        if (size(sponts.LS,2)<numsponts(1))||(size(sponts.MS,2)<numsponts(2))||(size(sponts.HS,2)<numsponts(3))||(size(sponts.HS,1)<numCFs||~exist('tabss','var'))
            disp('Saved population of AN fibers in ANpopulation.mat is too small - generating a new population');
            [sponts,tabss,trels] = generateANpopulation(numCFs,numsponts);
        end
    else
        [sponts,tabss,trels] = generateANpopulation(numCFs,numsponts);
        disp('Generating population of AN fibers, saved in ANpopulation.mat')
    end

    sponts_concat = zeros(numCFs,sum(numsponts));
    tabss_concat = zeros(numCFs,sum(numsponts));
    trels_concat = zeros(numCFs,sum(numsponts));

    for CFlp = 1:numCFs
        sponts_concat(CFlp,:) = [sponts.LS(CFlp,1:numsponts(1)) sponts.MS(CFlp,1:numsponts(2)) sponts.HS(CFlp,1:numsponts(3))];
        tabss_concat(CFlp,:) = [tabss.LS(CFlp,1:numsponts(1)) tabss.MS(CFlp,1:numsponts(2)) tabss.HS(CFlp,1:numsponts(3))];
        trels_concat(CFlp,:) = [trels.LS(CFlp,1:numsponts(1)) trels.MS(CFlp,1:numsponts(2)) trels.HS(CFlp,1:numsponts(3))];
    end

    implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
    noiseType = 1;  % 0 for fixed fGn (1 for variable fGn)
    expliketype = 1; % 1 for shifted softplus (preferred); 0 for no expontential-like function; 2 for shifted exponential; 3 for shifted Boltmann

    % stimulus parameters
    Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% START OF CODE WRITTEN BY ANDREW HINES %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Nfir = 30; % proportional to FIR filter length used for resampling: higher Nfir, better accuracy & longer comp time
    stimFs = resample(audio_signal, Fs, audio_sampling_rate, Nfir);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% END OF CODE WRITTEN BY ANDREW HINES %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    stimFs = stimFs';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% START OF CODE WRITTEN BY ME %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    clear audio_signal

    vihcFs_wavenet = resample(ihcogram(1, :), Fs, wavenet_Fs, Nfir);
    vihcFs_wavenet_len = length(vihcFs_wavenet);
    psth_tfs_wavenet = model_Synapse_BEZ2018a(vihcFs_wavenet, CFs(1), 1, 1 / Fs, noiseType, implnt,...
        sponts_concat(1, 1), tabss_concat(1, 1), trels_concat(1, 1), expliketype);
    clear vihcFs_wavenet


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

    psthbins_env = round(psthbinwidth_env * Fs);

    psth_env_wavenet = psth_tfs_wavenet(1:psthbins_env * floor(length(psth_tfs_wavenet) / psthbins_env));
    pr_env = sum(reshape(psth_env_wavenet, psthbins_env, length(psth_env_wavenet) / psthbins_env)); % pr of spike in each bin
    psth_env_wavenet = pr_env / psthbinwidth_env; % psth in units of spikes/s
    clear pr_env
    b_sp_env = spectrogram(psth_env_wavenet, w_env, env_overlap, env_windowsize, 1 / psthbinwidth_env);
    clear psth_env_wavenet
    synrate_env = abs(b_sp_env / env_windowsize / sqrt(sum(w_env.^2) / length(w_env)));
    clear b_sp_env
    cf_env_neurogram = synrate_env(1, :); % get the dc value

    psth_tfs_wavenet = psth_tfs_wavenet / psthbinwidth_tfs; % psth in units of spikes/s
    b_sp_tfs = spectrogram(psth_tfs_wavenet, w_tfs, tfs_overlap, tfs_windowsize, 1 / psthbinwidth_tfs);
    clear psth_tfs_wavenet
    synrate_tfs = abs(b_sp_tfs / tfs_windowsize / sqrt(sum(w_tfs.^2) / length(w_tfs)));
    clear b_sp_tfs
    cf_tfs_neurogram = synrate_tfs(1, :); % get the dc value

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% START OF CODE WRITTEN BY ME %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Memory allocation for neurograms
    env_neurogram = zeros(numCFs, length(cf_env_neurogram));
    env_neurogram_wavenet = zeros(numCFs, length(cf_env_neurogram));
    clear synrate_env cf_env_neurogram
    tfs_neurogram = zeros(numCFs, length(cf_tfs_neurogram));
    tfs_neurogram_wavenet = zeros(numCFs, length(cf_tfs_neurogram));
    clear synrate_tfs cf_tfs_neurogram

    % Memory allocation for SDRs
    SDRs = zeros(1, numCFs);

    stimFs_len = length(stimFs);

    for CFnum = 1:numCFs
        CF = CFs(CFnum);
        vihcFs = model_IHC_BEZ2018a(stimFs, CF, 1, 1 / Fs, stimFs_len,...
            cohcs(CFnum), cihcs(CFnum), species);
        psth_tfs = zeros(1, length(vihcFs));

        parfor spontnum = 1:sum(numsponts)
            single_psth_tfs = model_Synapse_BEZ2018a(vihcFs, CF, 1, 1 / Fs, noiseType,...
                implnt, sponts_concat(CFnum, spontnum), tabss_concat(CFnum, spontnum),...
                trels_concat(CFnum, spontnum), expliketype);
            psth_tfs = psth_tfs + single_psth_tfs;
            % psth_tfs is a reduction variable in this parfor loop
        end
        % Average the PSTH for all the simulated auditory nerves
        psth_tfs = psth_tfs / sum(numsponts);

        % Skip past the number of initial values in vihcFs and psth_tfs corresponding to
        % the duration of the WaveNet model's receptive field
        % This aligns vihcFs with vihcFs_wavenet and psth_tfs with the TFS PSTHs that
        % will be produced from vihcFs_wavenet
        vihcFs = vihcFs(length(vihcFs) - vihcFs_wavenet_len + 1:end);
        psth_tfs = psth_tfs(length(psth_tfs) - vihcFs_wavenet_len + 1:end);

        % Resample the WaveNet model's output to Fs for comparison with the auditory
        % periphery model's output, vihcFs, and to make it suitable for input to
        % the synapse model
        vihcFs_wavenet = resample(ihcogram(CFnum, :), Fs, wavenet_Fs, Nfir);

        SDRs(CFnum) = signal_to_distortion_ratio(vihcFs, vihcFs_wavenet);
        disp("CF: " + num2str(CF) + ', SDR: ' + num2str(SDRs(CFnum)))
        
        clear vihcFs

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% START OF CODE WRITTEN BY ANDREW HINES %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        psth_env = psth_tfs(1:psthbins_env * floor(length(psth_tfs) / psthbins_env));
        pr_env = sum(reshape(psth_env, psthbins_env, length(psth_env) / psthbins_env)); % pr of spike in each bin
        psth_env = pr_env / psthbinwidth_env; % psth in units of spikes/s
        clear pr_env
        b_sp_env = spectrogram(psth_env, w_env, env_overlap, env_windowsize, 1 / psthbinwidth_env);
        clear psth_env
        synrate_env = abs(b_sp_env / env_windowsize / sqrt(sum(w_env.^2) / length(w_env)));
        clear b_sp_env
        cf_env_neurogram = synrate_env(1, :); % get the dc value
        env_neurogram(CFnum, :) = cf_env_neurogram;

        psth_tfs = psth_tfs / psthbinwidth_tfs; % psth in units of spikes/s
        b_sp_tfs = spectrogram(psth_tfs, w_tfs, tfs_overlap, tfs_windowsize, 1 / psthbinwidth_tfs);
        clear psth_tfs
        synrate_tfs = abs(b_sp_tfs / tfs_windowsize / sqrt(sum(w_tfs.^2) / length(w_tfs)));
        clear b_sp_tfs
        cf_tfs_neurogram = synrate_tfs(1, :); % get the dc value
        tfs_neurogram(CFnum, :) = cf_tfs_neurogram;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%% START OF CODE WRITTEN BY ME %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        psth_tfs_wavenet = zeros(1, vihcFs_wavenet_len);

        parfor spontnum = 1:sum(numsponts)
            single_psth_tfs = model_Synapse_BEZ2018a(vihcFs_wavenet, CF, 1, 1 / Fs, noiseType,...
                implnt, sponts_concat(CFnum, spontnum), tabss_concat(CFnum, spontnum),...
                trels_concat(CFnum, spontnum), expliketype);
            psth_tfs_wavenet = psth_tfs_wavenet + single_psth_tfs;
            % psth_tfs_wavenet is a reduction variable in this parfor loop
        end

        clear vihcFs_wavenet

        % Average the PSTH for all the simulated auditory nerves
        psth_tfs_wavenet = psth_tfs_wavenet / sum(numsponts);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% START OF CODE WRITTEN BY ANDREW HINES %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        psth_env_wavenet = psth_tfs_wavenet(1:psthbins_env * floor(length(psth_tfs_wavenet) / psthbins_env));
        pr_env = sum(reshape(psth_env_wavenet, psthbins_env, length(psth_env_wavenet) / psthbins_env)); % pr of spike in each bin
        psth_env_wavenet = pr_env / psthbinwidth_env; % psth in units of spikes/s
        clear pr_env
        b_sp_env = spectrogram(psth_env_wavenet, w_env, env_overlap, env_windowsize, 1 / psthbinwidth_env);
        clear psth_env_wavenet
        synrate_env = abs(b_sp_env / env_windowsize / sqrt(sum(w_env.^2) / length(w_env)));
        clear b_sp_env
        cf_env_neurogram = synrate_env(1, :); % get the dc value
        env_neurogram_wavenet(CFnum, :) = cf_env_neurogram;

        psth_tfs_wavenet = psth_tfs_wavenet / psthbinwidth_tfs; % psth in units of spikes/s
        b_sp_tfs = spectrogram(psth_tfs_wavenet, w_tfs, tfs_overlap, tfs_windowsize, 1 / psthbinwidth_tfs);
        clear psth_tfs_wavenet
        synrate_tfs = abs(b_sp_tfs / tfs_windowsize / sqrt(sum(w_tfs.^2) / length(w_tfs)));
        clear b_sp_tfs
        cf_tfs_neurogram = synrate_tfs(1, :); % get the dc value
        tfs_neurogram_wavenet(CFnum, :) = cf_tfs_neurogram;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%% START OF CODE WRITTEN BY ME %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    clear stimFs ihcogram
    % Delete the MAT-file produced by the run_wavenet_model.py script
    delete(ihcogram_file)

    save(fullfile(outputDir, 'neurograms'), 'env_neurogram', 'env_neurogram_wavenet', 'tfs_neurogram', 'tfs_neurogram_wavenet')

    disp('ENV and TFS NSIM scores for windows of CFs')
    [ENV_NSIMs, mean_ENV_NSIM] = nsim_calc(env_neurogram, env_neurogram_wavenet);
    [TFS_NSIMs, mean_TFS_NSIM] = nsim_calc(tfs_neurogram, tfs_neurogram_wavenet);
    for row_num = 1:length(ENV_NSIMs)
        disp("Row " + num2str(row_num) + " mean ENV NSIM (" + num2str(ENV_NSIMs(row_num))...
            + ") and mean TFS NSIM (" + num2str(TFS_NSIMs(row_num)) + ")")
    end
    disp("Mean ENV NSIM: " + num2str(mean_ENV_NSIM))
    disp("Mean TFS NSIM: " + num2str(mean_TFS_NSIM))
    mean_SDR = mean(SDRs);
    disp("Mean SDR: " + num2str(mean_SDR))

    save(fullfile(outputDir, 'results'), 'SDRs', 'mean_SDR', 'ENV_NSIMs', 'mean_ENV_NSIM', 'TFS_NSIMs', 'mean_TFS_NSIM')
end
