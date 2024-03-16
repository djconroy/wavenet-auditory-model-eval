function [species, numCFs, CFs, cohcs, cihcs, numsponts, sponts_concat, tabss_concat,...
    trels_concat, implnt, noiseType, expliketype] = get_auditory_model_parameters()
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
end
