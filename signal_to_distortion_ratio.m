function SDR = signal_to_distortion_ratio(signal, distorted_signal)
    SDR = 20 * log10(rms(signal) / rms(signal - distorted_signal));
end
