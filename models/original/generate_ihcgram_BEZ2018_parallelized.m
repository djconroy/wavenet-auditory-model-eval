function vihc_mat = generate_ihcgram_BEZ2018_parallelized(stim,Fs_stim,species,numcfs,CFs)

cohcs  = ones(1,numcfs);  % normal ohc function
cihcs  = ones(1,numcfs);  % normal ihc function

% stimulus parameters
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
stim100k = resample(stim,Fs,Fs_stim).';

% PSTH parameters
nrep = 1;

vihc_mat=zeros(numcfs,length(stim100k));

% loop over all CFs
parfor CFind = 1:numcfs
    
    CF = CFs(CFind);
    cohc = cohcs(CFind);
    cihc = cihcs(CFind);
        
    vihc = model_IHC_BEZ2018a(stim100k,CF,nrep,1/Fs,length(stim100k),cohc,cihc,species);
    vihc_mat(CFind,:) = vihc;
    
end
