clear all
clc

addpath(fullfile('models', 'original'))

species = 1; % cat
numCFs = 80;
CFs = logspace(log10(125), log10(8000), numCFs);

spl = 60;

timitDatastore = audioDatastore("TIMIT", "IncludeSubfolders", true)

while hasdata(timitDatastore)
    [speech, info] = read(timitDatastore);

    speech = speech / rms(speech) * 20e-6 * 10 ^ (spl / 20);

    ihcogram = generate_ihcgram_BEZ2018_parallelized(speech, info.SampleRate, species, numCFs, CFs);

    outputDir = fullfile("outputs" + extractBetween(info.FileName, pwd, '.WAV'), 'original');
    mkdir(outputDir)
    save(fullfile(outputDir, 'ihcogram'), 'ihcogram')

    disp("Progress: " + (100 * progress(timitDatastore)) + "%")
end
