% Change the input arguments as you wish and then run the TIMIT pipeline from the terminal with
% matlab -batch "run_timit_pipeline"

speech_dataset_dir = fullfile(pwd, "datasets", "sample");
%speech_dataset_dir = "";
spl = 60;
snr = "N/A";

timit_pipeline(speech_dataset_dir, spl, snr);
