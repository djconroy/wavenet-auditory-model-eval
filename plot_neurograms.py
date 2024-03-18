import argparse
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from math import log10


CF = np.logspace(log10(125), log10(8000), 80, base=10.0).astype(int)


def get_time(sample_index, tfs):
    # Returns the time in seconds corresponding to the sample index
    # TFS Neurogram bin width = 16 samples at 100000 Hz
    # ENV Neurogram bin width = 640 samples at 100000 Hz
    return sample_index * 16 / 100000 if tfs else sample_index * 64 / 10000


def plot_neurogram(neurogram, ax: Axes, vmin, vmax, aspect, title, tfs=False):
    image = ax.imshow(neurogram,
                      cmap="YlGnBu_r",
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


parser = argparse.ArgumentParser(
    description="Produces plots of neurograms.")
parser.add_argument("neurograms_mat_file",
                    type=str,
                    nargs="?",
                    default="datasets/sample/neurograms.mat",
                    help="the name of a MAT-file storing the neurograms")
args = parser.parse_args()

neurograms = scipy.io.loadmat(args.neurograms_mat_file)
env_neurogram = neurograms["env_neurogram"]
env_neurogram_wavenet = neurograms["env_neurogram_wavenet"]

# Set the max and min values of the colorbar range to the max and min
# values in the ENV neurogram for the original model
vmax = np.max(env_neurogram)
vmin = np.min(env_neurogram)

# Plot ENV neurograms for the original and WaveNet models
fig, axs = plt.subplots(2, 1, sharex=True, sharey=True, layout="compressed")
image1 = plot_neurogram(env_neurogram, axs[0], vmin, vmax,
                        aspect=1, title="ENV Neurogram, Original")
image2 = plot_neurogram(env_neurogram_wavenet, axs[1], vmin, vmax,
                        aspect=1, title="ENV Neurogram, WaveNet")
axs[1].set_xlabel("Time (seconds)")
colorbar = fig.colorbar(image1, ax=axs, pad=0.028, aspect=30)
colorbar.ax.set_ylabel("Auditory Nerve Firing Intensity",
                       rotation=-90, verticalalignment="bottom")
plt.savefig("env_neurograms_plots.png")
