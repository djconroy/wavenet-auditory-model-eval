function timit_pipeline(speech_dataset_dir, spl, snr)
    arguments
        speech_dataset_dir (1, 1)
        spl (1, 1) = 60
        snr (1, 1) = "N/A"
    end

    warning('off', 'MATLAB:MKDIR:DirectoryExists')

    output_root_dir = "outputs";

    timit_Fs = 16000; % Sampling frequency of TIMIT speech signals
    Fs = 100000; % Sampling frequency of audio signal inputs to the auditory model
    wavenet_Fs = 16000; % Sampling frequency of audio signal inputs the WaveNet model is fitted to
    wavenet_RF = 2048; % Size of the receptive field of the WaveNet model

    resample_ratio = Fs / timit_Fs;

    auditory_periphery_delay = 91;

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
            audiowrite(noisy_speech_file, speech, speech_info.SampleRate, 'BitsPerSample', 64)

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

        [path, sentence] = fileparts(speech_info.FileName);

        % Open the TIMIT word and phoneme transcription files corresponding to the speech signal
        words_file_ID = fopen(fullfile(path, sentence) + ".WRD");
        phonemes_file_ID = fopen(fullfile(path, sentence) + ".PHN");

        words_info = textscan(words_file_ID, "%d%d%s");
        phonemes_info = textscan(phonemes_file_ID, "%d%d%s");

        fclose(words_file_ID);
        fclose(phonemes_file_ID);

        Nfir = 30; % proportional to FIR filter length used for resampling:
        % higher Nfir, better accuracy & longer comp time
        speech_Fs = resample(speech, Fs, speech_info.SampleRate, Nfir);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% END OF CODE WRITTEN BY ANDREW HINES %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        speech_Fs = speech_Fs';

        % Each line of a TIMIT transcription file has the form:
        % (start sample index) (end sample index) (phoneme/word/sentence)

        words_start_indices = words_info{1};
        words_end_indices = words_info{2};
        words = words_info{3};

        phonemes_start_indices = phonemes_info{1};
        phonemes_end_indices = phonemes_info{2};
        phonemes = phonemes_info{3};

        % Here's an example of consecutive lines in a phoneme transcription file:
        % 61538 66230 ae
        % 66230 69560 s
        % The indices of a phoneme's last sample and the next phoneme's first sample are the same

        last_word_num = size(words_start_indices, 1);
        first_word_num = 1;

        % Skip to the first word preceded by at least wavenet_RF - 1 zero-indexed audio samples
        while words_start_indices(first_word_num) < wavenet_RF - 1
            first_word_num = first_word_num + 1;
        end
        num_words = last_word_num - first_word_num + 1;

        words_start_indices = words_start_indices(first_word_num:last_word_num);
        words_end_indices = words_end_indices(first_word_num:last_word_num);
        words = words(first_word_num:last_word_num);

        % Subtract wavenet_RF - 1 from the word transciption file indices and then multiply
        % the indices by the resample ratio to make them align with vihc_Fs_wavenet
        % Then add auditory_periphery_delay to the indices to align them with the delayed
        % responses of the auditory model
        words_start_indices = words_start_indices - wavenet_RF + 1;
        words_start_indices = round(words_start_indices * resample_ratio);
        words_start_indices = words_start_indices + auditory_periphery_delay;
        words_end_indices = words_end_indices - wavenet_RF + 1;
        words_end_indices = round(words_end_indices * resample_ratio);
        words_end_indices = words_end_indices + auditory_periphery_delay;

        num_lines_phonemes_info = size(phonemes_start_indices, 1);
        last_phoneme_num = num_lines_phonemes_info - 1; % Subtract 1 to discount the end marker
        first_phoneme_num = 2; % Start with 2 to discount the start marker

        % Skip to the first phoneme preceded by at least wavenet_RF - 1 zero-indexed audio samples
        while phonemes_start_indices(first_phoneme_num) < wavenet_RF - 1
            first_phoneme_num = first_phoneme_num + 1;
        end
        num_phonemes = last_phoneme_num - first_phoneme_num + 1;

        phonemes_start_indices = phonemes_start_indices(first_phoneme_num:last_phoneme_num);
        phonemes_end_indices = phonemes_end_indices(first_phoneme_num:last_phoneme_num);
        phonemes = phonemes(first_phoneme_num:last_phoneme_num);

        % Subtract wavenet_RF - 1 from the phoneme transciption file indices and then multiply
        % the indices by the resample ratio to make them align with vihc_Fs_wavenet
        % Then add auditory_periphery_delay to the indices to align them with the delayed
        % responses of the auditory model
        phonemes_start_indices = phonemes_start_indices - wavenet_RF + 1;
        phonemes_start_indices = round(phonemes_start_indices * resample_ratio);
        phonemes_start_indices = phonemes_start_indices + auditory_periphery_delay;
        phonemes_end_indices = phonemes_end_indices - wavenet_RF + 1;
        phonemes_end_indices = round(phonemes_end_indices * resample_ratio);
        phonemes_end_indices = phonemes_end_indices + auditory_periphery_delay;

        % Memory allocation for SDRs
        SDRs = zeros(num_CFs, 1);
        words_SDRs = zeros(num_CFs, num_words);
        phonemes_SDRs = zeros(num_CFs, num_phonemes);

        vihc_Fs_wavenet = resample(ihcogram(1, :), Fs, wavenet_Fs, Nfir);
        vihc_Fs_wavenet_len = length(vihc_Fs_wavenet);
        psth_tfs_wavenet = model_Synapse_BEZ2018a(vihc_Fs_wavenet, CFs(1), 1, 1 / Fs,...
            noise_type, power_law_implementation, sponts_concat(1, 1), tabss_concat(1, 1),...
            trels_concat(1, 1), explike_type);
        clear vihc_Fs_wavenet

        % Create cell arrays to store neurograms for words
        words_env_neurograms = cell(1, num_words);
        words_env_neurograms_wavenet = cell(1, num_words);
        words_tfs_neurograms = cell(1, num_words);
        words_tfs_neurograms_wavenet = cell(1, num_words);

        for word_num = 1:num_words
            [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs_wavenet(...
                words_start_indices(word_num):words_end_indices(word_num)));

            % Memory allocation for word neurograms
            words_env_neurograms{word_num} = zeros(num_CFs, length(env_neurogram_row));
            words_env_neurograms_wavenet{word_num} = zeros(num_CFs, length(env_neurogram_row));
            words_tfs_neurograms{word_num} = zeros(num_CFs, length(tfs_neurogram_row));
            words_tfs_neurograms_wavenet{word_num} = zeros(num_CFs, length(tfs_neurogram_row));
        end

        % Create cell arrays to store neurograms for phonemes
        phonemes_env_neurograms = cell(1, num_phonemes);
        phonemes_env_neurograms_wavenet = cell(1, num_phonemes);
        phonemes_tfs_neurograms = cell(1, num_phonemes);
        phonemes_tfs_neurograms_wavenet = cell(1, num_phonemes);

        for phoneme_num = 1:num_phonemes
            [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs_wavenet(...
                phonemes_start_indices(phoneme_num):phonemes_end_indices(phoneme_num)));

            % Memory allocation for phoneme neurograms
            phonemes_env_neurograms{phoneme_num} = zeros(num_CFs, length(env_neurogram_row));
            phonemes_env_neurograms_wavenet{phoneme_num} =...
                zeros(num_CFs, length(env_neurogram_row));
            phonemes_tfs_neurograms{phoneme_num} = zeros(num_CFs, length(tfs_neurogram_row));
            phonemes_tfs_neurograms_wavenet{phoneme_num} =...
                zeros(num_CFs, length(tfs_neurogram_row));
        end

        [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs_wavenet);
        clear psth_tfs_wavenet

        % Memory allocation for sentence neurograms
        env_neurogram = zeros(num_CFs, length(env_neurogram_row));
        env_neurogram_wavenet = zeros(num_CFs, length(env_neurogram_row));
        tfs_neurogram = zeros(num_CFs, length(tfs_neurogram_row));
        tfs_neurogram_wavenet = zeros(num_CFs, length(tfs_neurogram_row));
        clear env_neurogram_row tfs_neurogram_row

        speech_Fs_len = length(speech_Fs);

        for CF_num = 1:num_CFs
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
            % will be produced from vihc_Fs_wavenet
            vihc_Fs = vihc_Fs(speech_Fs_len - vihc_Fs_wavenet_len + 1:end);
            psth_tfs = psth_tfs(speech_Fs_len - vihc_Fs_wavenet_len + 1:end);

            % Resample the WaveNet model's output to Fs for comparison with the auditory
            % periphery model's output, vihc_Fs, and to make it suitable for input to
            % the synapse model
            vihc_Fs_wavenet = resample(ihcogram(CF_num, :), Fs, wavenet_Fs, Nfir);

            SDRs(CF_num) = signal_to_distortion_ratio(vihc_Fs, vihc_Fs_wavenet);

            parfor word_num = 1:num_words
                words_SDRs(CF_num, word_num) = signal_to_distortion_ratio(...
                    vihc_Fs(words_start_indices(word_num):...
                            words_end_indices(word_num)),...
                    vihc_Fs_wavenet(words_start_indices(word_num):...
                                    words_end_indices(word_num)));
            end
            parfor phoneme_num = 1:num_phonemes
                phonemes_SDRs(CF_num, phoneme_num) = signal_to_distortion_ratio(...
                    vihc_Fs(phonemes_start_indices(phoneme_num):...
                            phonemes_end_indices(phoneme_num)),...
                    vihc_Fs_wavenet(phonemes_start_indices(phoneme_num):...
                                    phonemes_end_indices(phoneme_num)));
            end
            clear vihc_Fs

            parfor word_num = 1:num_words
                [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs(...
                    words_start_indices(word_num):words_end_indices(word_num)));
                words_env_neurograms{word_num}(CF_num, :) = env_neurogram_row;
                words_tfs_neurograms{word_num}(CF_num, :) = tfs_neurogram_row;
            end
            parfor phoneme_num = 1:num_phonemes
                [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs(...
                    phonemes_start_indices(phoneme_num):phonemes_end_indices(phoneme_num)));
                phonemes_env_neurograms{phoneme_num}(CF_num, :) = env_neurogram_row;
                phonemes_tfs_neurograms{phoneme_num}(CF_num, :) = tfs_neurogram_row;
            end
            [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs);
            clear psth_tfs
            env_neurogram(CF_num, :) = env_neurogram_row;
            tfs_neurogram(CF_num, :) = tfs_neurogram_row;
            clear env_neurogram_row tfs_neurogram_row

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

            parfor word_num = 1:num_words
                [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs_wavenet(...
                    words_start_indices(word_num):words_end_indices(word_num)));
                words_env_neurograms_wavenet{word_num}(CF_num, :) = env_neurogram_row;
                words_tfs_neurograms_wavenet{word_num}(CF_num, :) = tfs_neurogram_row;
            end
            parfor phoneme_num = 1:num_phonemes
                [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs_wavenet(...
                    phonemes_start_indices(phoneme_num):phonemes_end_indices(phoneme_num)));
                phonemes_env_neurograms_wavenet{phoneme_num}(CF_num, :) = env_neurogram_row;
                phonemes_tfs_neurograms_wavenet{phoneme_num}(CF_num, :) = tfs_neurogram_row;
            end
            [env_neurogram_row, tfs_neurogram_row] = gen_neurogram_row(psth_tfs_wavenet);
            clear psth_tfs_wavenet
            env_neurogram_wavenet(CF_num, :) = env_neurogram_row;
            tfs_neurogram_wavenet(CF_num, :) = tfs_neurogram_row;
            clear env_neurogram_row tfs_neurogram_row
        end

        clear speech_Fs ihcogram

        speech_dir = extractBetween(speech_info.FileName, pwd, ".WAV");
        spl_dir = "SPL" + num2str(spl);
        if isStringScalar(snr) && strcmp(snr, "N/A")
            snr_dir = "SNRNA";
        else
            snr_dir = "SNR" + num2str(snr);
        end
        output_dir = fullfile(output_root_dir + speech_dir, spl_dir, snr_dir);
        mkdir(output_dir)

        save(fullfile(output_dir, 'neurograms'),...
            'env_neurogram', 'env_neurogram_wavenet',...
            'tfs_neurogram', 'tfs_neurogram_wavenet',...
            'words_env_neurograms', 'words_env_neurograms_wavenet',...
            'words_tfs_neurograms', 'words_tfs_neurograms_wavenet',...
            'phonemes_env_neurograms', 'phonemes_env_neurograms_wavenet',...
            'phonemes_tfs_neurograms', 'phonemes_tfs_neurograms_wavenet')

        [env_NSIMs, mean_env_NSIM] = nsim(env_neurogram, env_neurogram_wavenet);
        [tfs_NSIMs, mean_tfs_NSIM] = nsim(tfs_neurogram, tfs_neurogram_wavenet);

        words_env_NSIMs = zeros(length(env_NSIMs), num_words);
        words_tfs_NSIMs = zeros(length(tfs_NSIMs), num_words);
        words_mean_env_NSIMs = zeros(1, num_words);
        words_mean_tfs_NSIMs = zeros(1, num_words);
        parfor word_num = 1:num_words
            if size(words_env_neurograms{word_num}, 2) > 1
                [words_env_NSIMs(:, word_num), words_mean_env_NSIMs(word_num)] =...
                    nsim(words_env_neurograms{word_num},...
                         words_env_neurograms_wavenet{word_num});
            end
            [words_tfs_NSIMs(:, word_num), words_mean_tfs_NSIMs(word_num)] =...
                nsim(words_tfs_neurograms{word_num},...
                     words_tfs_neurograms_wavenet{word_num});
        end

        phonemes_env_NSIMs = zeros(length(env_NSIMs), num_phonemes);
        phonemes_tfs_NSIMs = zeros(length(tfs_NSIMs), num_phonemes);
        phonemes_mean_env_NSIMs = zeros(1, num_phonemes);
        phonemes_mean_tfs_NSIMs = zeros(1, num_phonemes);
        parfor phoneme_num = 1:num_phonemes
            if size(phonemes_env_neurograms{phoneme_num}, 2) > 1
                [phonemes_env_NSIMs(:, phoneme_num), phonemes_mean_env_NSIMs(phoneme_num)] =...
                    nsim(phonemes_env_neurograms{phoneme_num},...
                         phonemes_env_neurograms_wavenet{phoneme_num});
            end
            [phonemes_tfs_NSIMs(:, phoneme_num), phonemes_mean_tfs_NSIMs(phoneme_num)] =...
                nsim(phonemes_tfs_neurograms{phoneme_num},...
                     phonemes_tfs_neurograms_wavenet{phoneme_num});
        end

        mean_SDR = mean(SDRs);
        words_mean_SDRs = mean(words_SDRs, 1);
        phonemes_mean_SDRs = mean(phonemes_SDRs, 1);

        save(fullfile(output_dir, 'comparison_metrics'),...
            'words', 'phonemes',...
            'SDRs', 'mean_SDR',...
            'words_SDRs', 'words_mean_SDRs',...
            'phonemes_SDRs', 'phonemes_mean_SDRs',...
            'env_NSIMs', 'mean_env_NSIM',...
            'tfs_NSIMs', 'mean_tfs_NSIM',...
            'words_env_NSIMs', 'words_mean_env_NSIMs',...
            'words_tfs_NSIMs', 'words_mean_tfs_NSIMs',...
            'phonemes_env_NSIMs', 'phonemes_mean_env_NSIMs',...
            'phonemes_tfs_NSIMs', 'phonemes_mean_tfs_NSIMs')

        disp("Progress: " + (100 * progress(speech_datastore)) + "%")
    end

    % Delete the MAT-file produced by the run_wavenet_model.py script
    delete(ihcogram_file)
end
