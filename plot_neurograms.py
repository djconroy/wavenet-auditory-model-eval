import argparse
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from math import log10


CF = np.round(np.logspace(log10(125), log10(8000), 80, base=10.0)).astype(int)


def get_time(sample_index, tfs):
    # Returns the time in seconds corresponding to the sample index
    # TFS Neurogram bin width = 16 samples at 100000 Hz
    # ENV Neurogram bin width = 640 samples at 100000 Hz
    return sample_index * 16 / 100000 if tfs else sample_index * 64 / 10000


def plot_neurogram(neurogram, ax: Axes, vmin, vmax, aspect, title, tfs=False):
    image = ax.imshow(neurogram,
                      cmap="RdYlBu_r" if tfs else "YlGnBu_r",
                      aspect=aspect,
                      interpolation="none",
                      vmin=vmin,
                      vmax=vmax,
                      origin="lower")
    xticks = ax.get_xticks()
    ax.set_xticks(xticks, labels=[get_time(xtick, tfs) for xtick in xticks])
    ax.margins(0)
    yticks = [0, 16, 32, 48, 64, 79]
    ax.set_yticks(yticks, labels=[CF[ytick] for ytick in yticks])
    ax.set_ylabel("CFs", rotation="horizontal", verticalalignment="center")
    ax.set_title(title)
    return image


parser = argparse.ArgumentParser(description="Produces plots of neurograms.")
parser.add_argument("neurograms_mat_file",
                    type=str,
                    nargs="?",
                    default="datasets/sample/neurograms_SPL60_SNR25.mat",
                    help="the name of a MAT-file storing the neurograms")
args = parser.parse_args()

neurograms = scipy.io.loadmat(args.neurograms_mat_file, squeeze_me=True)
env_neurogram = neurograms["env_neurogram"]
env_neurogram_wavenet = neurograms["env_neurogram_wavenet"]
words_env_neurograms = neurograms["words_env_neurograms"]
phonemes_env_neurograms = neurograms["phonemes_env_neurograms"]

# Set the max and min values of the colorbar range to the max and min
# values in the ENV neurogram for the original model
vmax = np.max(env_neurogram)
vmin = np.min(env_neurogram)

# Plot ENV neurograms for the original and WaveNet models
fig, axs = plt.subplots(2, 1, sharex=True, sharey=True, layout="compressed")
image1 = plot_neurogram(env_neurogram, axs[0], vmin, vmax,
                        aspect=2, title="ENV Neurogram, Original")
image2 = plot_neurogram(env_neurogram_wavenet, axs[1], vmin, vmax,
                        aspect=2, title="ENV Neurogram, WaveNet")
axs[1].set_xlabel("Time (seconds)")
colorbar = fig.colorbar(image1, ax=axs, pad=0.028, aspect=30)
colorbar.ax.set_ylabel("Auditory Nerve Firing Intensity",
                       rotation=-90, verticalalignment="bottom")
plt.savefig("env_neurograms_plots.png")
#plt.show()


# Set the max and min values to the max and min values in the ENV neurogram for the word wash
vmax = np.max(words_env_neurograms[7])
vmin = np.min(words_env_neurograms[7])

fig = plt.figure(figsize=(6, 6), layout="compressed")
spec = fig.add_gridspec(2, 3)
ax0 = fig.add_subplot(spec[0, :])
ax0.set_xlabel("Time (s)")
plot_neurogram(words_env_neurograms[7], ax0, vmin, vmax, aspect=0.5,
               title="ENV Neurogram for “wash”")
ax10 = fig.add_subplot(spec[1, 0])
ax10.set_xlabel("Time (s)")
plot_neurogram(phonemes_env_neurograms[22], ax10, vmin, vmax, aspect=0.5, title="w")
ax11 = fig.add_subplot(spec[1, 1])
ax11.set_xlabel("Time (s)")
plot_neurogram(phonemes_env_neurograms[23], ax11, vmin, vmax, aspect=0.5, title="ao")
ax12 = fig.add_subplot(spec[1, 2])
ax12.set_xlabel("Time (s)")
plot_neurogram(phonemes_env_neurograms[24], ax12, vmin, vmax, aspect=0.5, title="sh")
plt.savefig("env_neurogram_word_wash.png")
#plt.show()
