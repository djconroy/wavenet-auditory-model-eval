function run_pipeline(spl, snr)
    arguments
        spl (1, 1) = 60
        snr (1, 1) = "N/A"
    end

    addpath(fullfile('models', 'original'))


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% START OF CODE NOT WRITTEN BY ME %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    species = 1; % cat
    numCFs = 80;
    CFs = logspace(log10(125), log10(8000), numCFs);

    cohcs  = ones(1, numCFs);  % normal ohc function
    cihcs  = ones(1, numCFs);  % normal ihc function

    numsponts_healthy = [10 10 30]; % Number of low-spont, medium-spont, and high-spont fibers at each CF in a healthy AN

    if exist('ANpopulation.mat','file')
        load('ANpopulation.mat');
        disp('Loading existing population of AN fibers saved in ANpopulation.mat')
        if (size(sponts.LS,2)<numsponts_healthy(1))||(size(sponts.MS,2)<numsponts_healthy(2))||(size(sponts.HS,2)<numsponts_healthy(3))||(size(sponts.HS,1)<numCFs||~exist('tabss','var'))
            disp('Saved population of AN fibers in ANpopulation.mat is too small - generating a new population');
            [sponts,tabss,trels] = generateANpopulation(numCFs,numsponts_healthy);
        end
    else
        [sponts,tabss,trels] = generateANpopulation(numCFs,numsponts_healthy);
        disp('Generating population of AN fibers, saved in ANpopulation.mat')
    end

    implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
    noiseType = 1;  % 0 for fixed fGn (1 for variable fGn)
    expliketype = 1; % 1 for shifted softplus (preferred); 0 for no expontential-like function; 2 for shifted exponential; 3 for shifted Boltmann

    % PSTH parameters
    nrep = 1;
    psthbinwidth_mr = 100e-6; % mean-rate binwidth in seconds;
    windur_ft=32;
    smw_ft = hamming(windur_ft);
    windur_mr=128;
    smw_mr = hamming(windur_mr);

    % stimulus parameters
    Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% END OF CODE NOT WRITTEN BY ME %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    wavenetFs = 16000; % Sampling rate of audio signal inputs the WaveNet model is fitted to
    wavenetRF = 2048; % Size of the receptive field of the WaveNet model

    timitDatastore = audioDatastore("TIMIT", "IncludeSubfolders", true)

    while hasdata(timitDatastore)
        [speech, speech_info] = read(timitDatastore);

        % Normalize the speech signal to the specified sound pressure level (SPL)
        speech = speech / rms(speech) * 20e-6 * 10 ^ (spl / 20);

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


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% START OF CODE NOT WRITTEN BY ME %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        stim100k = resample(speech, Fs, speech_info.SampleRate).';
        T  = length(stim100k)/Fs;  % stimulus duration in seconds

        simdur = ceil(T*1.2/psthbinwidth_mr)*psthbinwidth_mr;

        % Memory Allocation added by Wissam
        % Run model_IHC_BEZ2018a and model_Synapse_BEZ2018a functions to estimate the
        % size of neurogram_ft and other variables.
        vihc = model_IHC_BEZ2018a(stim100k,CFs(1),nrep,1/Fs,simdur,cohcs(1),cihcs(1),species);
        vihc_mat=zeros(numCFs,length(vihc));
        neurogram_ft=zeros(numCFs,length(vihc));
        psth_ft = model_Synapse_BEZ2018a(vihc,CFs(1),nrep,1/Fs,noiseType,implnt,sponts.LS(1,1),tabss.LS(1,1),trels.LS(1,1),expliketype);
        psthbins = round(psthbinwidth_mr*Fs);  % number of psth_ft bins per psth bin
        psth_mr = sum(reshape(psth_ft,psthbins,length(psth_ft)/psthbins));
        neurogram_mr=zeros(numCFs,length(psth_mr));
        clear vihc psth_ft psthbins psth_mr

        parfor CFlp = 1:numCFs
            CF = CFs(CFlp);
            cohc = cohcs(CFlp);
            cihc = cihcs(CFlp);

            numsponts = round([1 1 1].*numsponts_healthy); % Healthy AN
            sponts_concat = [sponts.LS(CFlp,1:numsponts(1)) sponts.MS(CFlp,1:numsponts(2)) sponts.HS(CFlp,1:numsponts(3))];
            tabss_concat = [tabss.LS(CFlp,1:numsponts(1)) tabss.MS(CFlp,1:numsponts(2)) tabss.HS(CFlp,1:numsponts(3))];
            trels_concat = [trels.LS(CFlp,1:numsponts(1)) trels.MS(CFlp,1:numsponts(2)) trels.HS(CFlp,1:numsponts(3))];

            vihc = model_IHC_BEZ2018a(stim100k,CF,nrep,1/Fs,simdur,cohc,cihc,species);
            vihc_mat(CFlp,:) = vihc;

            for spontlp = 1:sum(numsponts)
                spont = sponts_concat(spontlp);
                tabs = tabss_concat(spontlp);
                trel = trels_concat(spontlp);

                psth_ft = model_Synapse_BEZ2018a(vihc,CF,nrep,1/Fs,noiseType,implnt,spont,tabs,trel,expliketype);
                psthbins = round(psthbinwidth_mr*Fs);  % number of psth_ft bins per psth bin
                psth_mr = sum(reshape(psth_ft,psthbins,length(psth_ft)/psthbins));

                neurogram_ft(CFlp,:) = neurogram_ft(CFlp,:)+filter(smw_ft,1,psth_ft);
                neurogram_mr(CFlp,:) = neurogram_mr(CFlp,:)+filter(smw_mr,1,psth_mr);
            end
        end

        vihc_mat = vihc_mat(:,1:round(T*Fs));
        vihc_mat = resample(vihc_mat', wavenetFs, Fs)';
        vihc_mat = vihc_mat(:, wavenetRF:end); % discard first 2047 samples for better comparion with WaveNet model
        neurogram_ft = neurogram_ft(:,1:ceil(windur_ft/2):end); % 50% overlap in Hamming window
        neurogram_mr = neurogram_mr(:,1:ceil(windur_mr/2):end); % 50% overlap in Hamming window

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% END OF CODE NOT WRITTEN BY ME %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % Run the WaveNet model using the same speech signal input
        system("python run_wavenet_model.py --file " + speech_file + " --spl " + num2str(spl))
        % Load the output produced by the WaveNet model, which includes the 'ihcogram' variable
        load('tmp.mat')

        if noisy_speech_file == "noisy.wav"
            delete noisy_speech_file
        end

        % Open the TIMIT Phone (PHN) file corresponding to the speech signal
        fileID = fopen(speech_info.FileName(1:end - 3) + "PHN");

        % The following line of code is a slightly modified version of code written by Andrew Hines
        phone_info = textscan(fileID, "%d%d%s");
        % phone_info{1}, phone_info{2} and phone_info{3} are column vectors containing
        % start indices, end indices, and phoneme types, respectively

        fclose(fileID);

        % TIMIT Phone data consists of start points, end points and types of phones/phonemes
        % Each line has the form: (start sample index) (end sample index) (phoneme type)
        % Here's an example of consecutive lines:
        % 61538 66230 ae
        % 66230 69560 s
        % The indices of a phone's last sample and the next phone's first sample are the same

        numLines = size(phone_info{1}, 1);
        lastPhoneIndex = numLines - 1; % 2nd last line since the last one contains an end marker
        firstPhoneIndex = 2; % 2nd line since the first one contains a start marker

        % Skip to the first phone preceded by at least (wavenetRF - 1) zero-indexed audio samples
        while phone_info{1}(firstPhoneIndex) < wavenetRF - 1
            firstPhoneIndex = firstPhoneIndex + 1;
        end

        SDRs = zeros(lastPhoneIndex - firstPhoneIndex + 1, numCFs);
        diffs = vihc_mat - ihcogram;

        % Calculate SDRs for each phone/phoneme and each CF
        for phoneIndex = firstPhoneIndex:lastPhoneIndex
            % Subtract (wavenetRF - 1) from the indices of the first and last samples of phones
            phoneStart = phone_info{1}(phoneIndex) - wavenetRF + 1;
            phoneEnd = phone_info{2}(phoneIndex) - wavenetRF + 1;

            for CFnum = 1:numCFs
                SDRs(phoneIndex - firstPhoneIndex + 1, CFnum) = 10 * log10(...
                sum(vihc_mat(CFnum, phoneStart:phoneEnd).^2)...
                /sum(diffs(CFnum, phoneStart:phoneEnd).^2));
            end
        end
    end
end
