function amplified_audio_signal = amplify_audio_signal(audio_signal, spl)
    amplified_audio_signal = audio_signal / rms(audio_signal) * 20e-6 * 10 ^ (spl / 20);
end
