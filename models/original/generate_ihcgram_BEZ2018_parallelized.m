function vihc_mat = generate_ihcgram_BEZ2018_parallelized(stim,Fs_stim,species,numcfs,CFs)

cohcs  = ones(1,numcfs);  % normal ohc function
cihcs  = ones(1,numcfs);  % normal ihc function

% stimulus parameters
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
stim100k = resample(stim,Fs,Fs_stim).';
T  = length(stim100k)/Fs;  % stimulus duration in seconds

% PSTH parameters
nrep = 1;
psthbinwidth_mr = 100e-6; % mean-rate binwidth in seconds;

simdur = ceil(T*1.2/psthbinwidth_mr)*psthbinwidth_mr;

% Run model_IHC_BEZ2018 function to estimate the size of vihc and variable.  
vihc = model_IHC_BEZ2018a(stim100k,CFs(1),nrep,1/Fs,simdur,cohcs(1),cihcs(1),species);
vihc_mat=zeros(numcfs,length(vihc));
clear vihc

% loop over all CFs
parfor CFind = 1:numcfs
    
    CF = CFs(CFind);
    cohc = cohcs(CFind);
    cihc = cihcs(CFind);
        
    vihc = model_IHC_BEZ2018a(stim100k,CF,nrep,1/Fs,simdur,cohc,cihc,species);
    vihc_mat(CFind,:) = vihc;
    
end

vihc_mat = vihc_mat(:,1:round(T*Fs));
