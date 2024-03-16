function gen_comparison_metrics(audio_file_name, spl, snr)
    arguments
        audio_file_name (1, 1)
        spl (1, 1) = 60
        snr (1, 1) = "N/A"
    end

    Fs = 100e3; % Sampling frequency of audio signal inputs to the auditory model
    wavenet_Fs = 16000; % Sampling frequency of audio signal inputs the WaveNet model is fitted to
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
    %%%%%%% START OF CODE WRITTEN BY ANDREW HINES %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Nfir = 30; % proportional to FIR filter length used for resampling:
    % higher Nfir, better accuracy & longer comp time
    stimFs = resample(audio_signal, Fs, audio_sampling_rate, Nfir);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% END OF CODE WRITTEN BY ANDREW HINES %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    stimFs = stimFs';

    clear audio_signal

    [species, numCFs, CFs, cohcs, cihcs, numsponts, sponts_concat, tabss_concat,...
    trels_concat, implnt, noiseType, expliketype] = get_auditory_model_parameters();

    vihcFs_wavenet = resample(ihcogram(1, :), Fs, wavenet_Fs, Nfir);
    vihcFs_wavenet_len = length(vihcFs_wavenet);
    psth_tfs_wavenet = model_Synapse_BEZ2018a(vihcFs_wavenet, CFs(1), 1, 1 / Fs, noiseType,...
        implnt, sponts_concat(1, 1), tabss_concat(1, 1), trels_concat(1, 1), expliketype);
    clear vihcFs_wavenet

    [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs_wavenet);
    clear psth_tfs_wavenet

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

    save(fullfile(outputDir, 'neurograms'), 'env_neurogram', 'env_neurogram_wavenet',...
        'tfs_neurogram', 'tfs_neurogram_wavenet')

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

    save(fullfile(outputDir, 'results'), 'SDRs', 'mean_SDR', 'ENV_NSIMs', 'mean_ENV_NSIM',...
        'TFS_NSIMs', 'mean_TFS_NSIM')
end
