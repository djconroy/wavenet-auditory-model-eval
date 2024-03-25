import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from math import log10

nsim_CFs = np.logspace(log10(250), log10(8000), 30, base=10.0)
wavenet_CFs = np.logspace(log10(125), log10(8000), 80, base=10.0)

fig, ax = plt.subplots(figsize=(7, 1.6))
ax.set_xscale('log')
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.set_ylim(0, 1)
ax.minorticks_off()
#ax.tick_params(axis='x', labelrotation=90)
ax.scatter(wavenet_CFs, 0.75 * np.ones(80), c="royalblue", marker="|", linewidths=2)
ax.plot(wavenet_CFs[-7:], [0.75] * 7, c="aqua", marker=".")
ax.scatter(nsim_CFs, 0.5 * np.ones(30), c="royalblue", marker="|", linewidths=2)
ax.plot(nsim_CFs[-3:], [0.5] * 3, c="red", marker=".")
ax.scatter(wavenet_CFs, 0.25 * np.ones(80), c="royalblue", marker="|", linewidths=2)
ax.plot(wavenet_CFs[-5:], [0.25] * 5, c="lime", marker=".")
ax.set_xticks([round(cf) for cf in np.logspace(log10(125), log10(8000), 7, base=10.0)])
ax.set_xlabel("Characteristic Frequencies (Hz)")
ax.set_yticks([])
plt.tight_layout()
plt.savefig("nsim_window_scales.png")
#plt.show()
