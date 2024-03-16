% Change the input arguments as you wish and then run the pipeline from the terminal with
% matlab -batch "run_pipeline"

audio_file_name = fullfile(pwd, "datasets", "sample", "LDC93S1.wav");
%audio_file_name = "";
spl = 60;
snr = "N/A";

gen_comparison_metrics(audio_file_name, spl, snr);
