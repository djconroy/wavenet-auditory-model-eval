function timit_pipeline(speech_dataset_dir, spl, snr)
    arguments
        speech_dataset_dir (1, 1)
        spl (1, 1) = 60
        snr (1, 1) = "N/A"
    end

    warning('off', 'MATLAB:MKDIR:DirectoryExists')

    output_dir = "outputs";
    mkdir(output_dir)

    Fs = 100000; % Sampling frequency of audio signal inputs to the auditory model
    wavenet_Fs = 16000; % Sampling frequency of audio signal inputs the WaveNet model is fitted to
    wavenet_RF = 2048; % Size of the receptive field of the WaveNet model

    [species, num_CFs, CFs, Cohcs, Cihcs, num_sponts, sponts_concat,...
    tabss_concat, trels_concat, power_law_implementation, noise_type,...
    explike_type] = get_auditory_model_parameters();

    addpath(fullfile('models', 'original'))

    speech_datastore = audioDatastore(speech_dataset_dir, "IncludeSubfolders", true);

    while hasdata(speech_datastore)
        [speech, speech_info] = read(speech_datastore);

        speech = amplify_audio_signal(speech, spl);

        if isStringScalar(snr) && strcmp(snr, "N/A")
            noisy_speech_file = "";
            speech_file = speech_info.FileName;
        else
            % Add white Gaussian noise to the speech signal in accordance with
            % the specified signal-to-noise ratio (SNR)
            speech = awgn(speech, snr, 'measured');

            % Save the noisified speech signal to a temporary WAVE file for
            % use by the WaveNet model
            noisy_speech_file = "noisy.wav";
            audiowrite(noisy_speech_file, speech, speech_info.SampleRate)

            speech_file = noisy_speech_file;
        end

        % Run the WaveNet model
        system("python run_wavenet_model.py --file " + speech_file + " --spl " + num2str(spl));
        disp('Finished running WaveNet model')
        ihcogram_file = "ihcogram.mat";
        % Load the output produced by the WaveNet model, which includes the 'ihcogram' variable
        load(ihcogram_file)

        if noisy_speech_file == "noisy.wav"
            % Delete the temporary WAVE file of the noisified speech signal after
            % it's used by the WaveNet model
            delete(noisy_speech_file)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% START OF CODE WRITTEN BY ANDREW HINES %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Nfir = 30; % proportional to FIR filter length used for resampling:
        % higher Nfir, better accuracy & longer comp time
        speech_Fs = resample(speech, Fs, speech_info.SampleRate, Nfir);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% END OF CODE WRITTEN BY ANDREW HINES %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        speech_Fs = speech_Fs';

        vihc_Fs_wavenet = resample(ihcogram(1, :), Fs, wavenet_Fs, Nfir);
        vihc_Fs_wavenet_len = length(vihc_Fs_wavenet);
        psth_tfs_wavenet = model_Synapse_BEZ2018a(vihc_Fs_wavenet, CFs(1), 1, 1 / Fs,...
            noise_type, power_law_implementation, sponts_concat(1, 1), tabss_concat(1, 1),...
            trels_concat(1, 1), explike_type);
        clear vihc_Fs_wavenet

        [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs_wavenet);
        clear psth_tfs_wavenet

        % Memory allocation for neurograms
        env_neurogram = zeros(num_CFs, length(env_neurogram_row));
        env_neurogram_wavenet = zeros(num_CFs, length(env_neurogram_row));
        clear env_neurogram_row
        tfs_neurogram = zeros(num_CFs, length(tfs_neurogram_row));
        tfs_neurogram_wavenet = zeros(num_CFs, length(tfs_neurogram_row));
        clear tfs_neurogram_row

        % Memory allocation for SDRs
        SDRs = zeros(1, num_CFs);

        speech_Fs_len = length(speech_Fs);

        for CF_num =1:num_CFs
            CF = CFs(CF_num);
            vihc_Fs = model_IHC_BEZ2018a(speech_Fs, CF, 1, 1 / Fs, speech_Fs_len,...
                Cohcs(CF_num), Cihcs(CF_num), species);
            psth_tfs = zeros(1, speech_Fs_len);

            parfor spont_num = 1:sum(num_sponts)
                single_psth_tfs = model_Synapse_BEZ2018a(vihc_Fs, CF, 1, 1 / Fs,...
                    noise_type, power_law_implementation,...
                    sponts_concat(CF_num, spont_num), tabss_concat(CF_num, spont_num),...
                    trels_concat(CF_num, spont_num), explike_type);
                psth_tfs = psth_tfs + single_psth_tfs;
                % psth_tfs is a reduction variable in this parfor loop
            end
            % Average the PSTH for all the simulated auditory nerves
            psth_tfs = psth_tfs / sum(num_sponts);

            % Skip past the number of initial values in vihc_Fs and psth_tfs corresponding to
            % the duration of the WaveNet model's receptive field
            % This aligns vihc_Fs with vihc_Fs_wavenet and psth_tfs with the TFS PSTHs that
            % will be produced from vihcFs_wavenet
            vihc_Fs = vihc_Fs(speech_Fs_len - vihc_Fs_wavenet_len + 1:end);
            psth_tfs = psth_tfs(speech_Fs_len - vihc_Fs_wavenet_len + 1:end);

            % Resample the WaveNet model's output to Fs for comparison with the auditory
            % periphery model's output, vihc_Fs, and to make it suitable for input to
            % the synapse model
            vihc_Fs_wavenet = resample(ihcogram(CF_num, :), Fs, wavenet_Fs, Nfir);

            SDRs(CF_num) = signal_to_distortion_ratio(vihc_Fs, vihc_Fs_wavenet);

            clear vihc_Fs

            [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs);
            clear psth_tfs
            env_neurogram(CF_num, :) = env_neurogram_row;
            tfs_neurogram(CF_num, :) = tfs_neurogram_row;

            psth_tfs_wavenet = zeros(1, vihc_Fs_wavenet_len);

            parfor spont_num = 1:sum(num_sponts)
                single_psth_tfs = model_Synapse_BEZ2018a(vihc_Fs_wavenet, CF, 1, 1 / Fs,...
                    noise_type, power_law_implementation,...
                    sponts_concat(CF_num, spont_num), tabss_concat(CF_num, spont_num),...
                    trels_concat(CF_num, spont_num), explike_type);
                psth_tfs_wavenet = psth_tfs_wavenet + single_psth_tfs;
                % psth_tfs_wavenet is a reduction variable in this parfor loop
            end
            clear vihc_Fs_wavenet
            % Average the PSTH for all the simulated auditory nerves
            psth_tfs_wavenet = psth_tfs_wavenet / sum(num_sponts);

            [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs_wavenet);
            clear psth_tfs_wavenet
            env_neurogram_wavenet(CF_num, :) = env_neurogram_row;
            tfs_neurogram_wavenet(CF_num, :) = tfs_neurogram_row;
        end

        clear speech_Fs ihcogram

        output_subdir = output_dir + extractBetween(speech_info.FileName, pwd, ".WAV");
        mkdir(output_subdir)

        save(fullfile(output_subdir, 'neurograms'),...
            'env_neurogram', 'env_neurogram_wavenet',...
            'tfs_neurogram', 'tfs_neurogram_wavenet')

        [env_NSIMs, mean_env_NSIM] = nsim(env_neurogram, env_neurogram_wavenet);
        [tfs_NSIMs, mean_tfs_NSIM] = nsim(tfs_neurogram, tfs_neurogram_wavenet);
        mean_SDR = mean(SDRs);

        save(fullfile(output_subdir, 'comparison_metrics'),...
            'SDRs', 'mean_SDR',...
            'env_NSIMs', 'mean_env_NSIM',...
            'tfs_NSIMs', 'mean_tfs_NSIM')

        disp("Progress: " + (100 * progress(speech_datastore)) + "%")
    end

    % Delete the MAT-file produced by the run_wavenet_model.py script
    delete(ihcogram_file)
end
