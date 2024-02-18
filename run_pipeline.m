function run_pipeline(spl, snr)
    arguments
        spl (1, 1) = 60
        snr (1, 1) = "N/A"
    end

    addpath(fullfile('models', 'original'))

    species = 1; % cat

    % Use the 80 CFs used by the WaveNet model
    numCFs = 80;
    CFs = logspace(log10(125), log10(8000), numCFs);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% START OF CODE NOT WRITTEN BY ME %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


    timitDatastore = audioDatastore("TIMIT", "IncludeSubfolders", true)

    while hasdata(timitDatastore)
        [audio, info] = read(timitDatastore);

        % Normalize the audio signal to the specified sound pressure level (SPL)
        audio = audio / rms(audio) * 20e-6 * 10 ^ (spl / 20);

        if ~(isStringScalar(snr) && strcmp(snr, "N/A"))
            % Add white Gaussian noise to the audio in accordance with the
            % specified signal-to-noise ratio (SNR)
            audio = awgn(audio, snr, 'measured');
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% START OF CODE NOT WRITTEN BY ME %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        stim100k = resample(audio,Fs,info.SampleRate).';
        T  = length(stim100k)/Fs;  % stimulus duration in seconds

        simdur = ceil(T*1.2/psthbinwidth_mr)*psthbinwidth_mr;

        %% Memory Allocation added by Wissam
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
        neurogram_ft = neurogram_ft(:,1:ceil(windur_ft/2):end); % 50% overlap in Hamming window
        neurogram_mr = neurogram_mr(:,1:ceil(windur_mr/2):end); % 50% overlap in Hamming window

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% END OF CODE NOT WRITTEN BY ME %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
