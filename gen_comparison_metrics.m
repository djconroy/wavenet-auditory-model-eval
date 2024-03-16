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

    outputDir = "outputs";
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

    [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs_wavenet);

    % Memory allocation for neurograms
    env_neurogram = zeros(numCFs, length(env_neurogram_row));
    env_neurogram_wavenet = zeros(numCFs, length(env_neurogram_row));
    clear env_neurogram_row
    tfs_neurogram = zeros(numCFs, length(tfs_neurogram_row));
    tfs_neurogram_wavenet = zeros(numCFs, length(tfs_neurogram_row));
    clear tfs_neurogram_row

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

        [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs);
        env_neurogram(CFnum, :) = env_neurogram_row;
        tfs_neurogram(CFnum, :) = tfs_neurogram_row;

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

        [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs_wavenet);
        env_neurogram_wavenet(CFnum, :) = env_neurogram_row;
        tfs_neurogram_wavenet(CFnum, :) = tfs_neurogram_row;
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
