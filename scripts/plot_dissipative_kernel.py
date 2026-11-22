"""Plot the finite-image Ohmic kernel printed in Appendix C, Eq. (C3)."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


BETA = 256
MIN_SEPARATION = 1
KERNEL_OFFSET = 0.25
FIGURE_SIZE = (10, 4)
FIGURE_DPI = 180
LINE_WIDTH = 1.8
FONT_SIZE = 11
CURVE_COLOR = "#176b91"
GRID_ALPHA = 0.2
OUTPUT_DIRECTORY = Path(__file__).resolve().parents[1] / "figures"


separation = np.arange(MIN_SEPARATION, BETA)
images = np.arange(-BETA // 2, BETA // 2 + 1)
kernel = np.sum(
    1.0 / ((separation[:, None] + BETA * images[None, :]) ** 2 - KERNEL_OFFSET),
    axis=1,
)

assert np.all(np.isfinite(kernel)) and np.all(kernel > 0)

plt.rcParams.update({"font.size": FONT_SIZE})
figure, axes = plt.subplots(1, 2, figsize=FIGURE_SIZE, layout="constrained")
for axis in axes:
    axis.plot(separation, kernel, color=CURVE_COLOR, linewidth=LINE_WIDTH)
    axis.set_xlabel(r"Imaginary-time separation $r=|j-j'|$")
    axis.set_ylabel(r"$D_\beta(r)$")
    axis.set_xlim(MIN_SEPARATION, BETA - MIN_SEPARATION)
    axis.grid(alpha=GRID_ALPHA)
axes[0].set_title("Linear scale")
axes[1].set_yscale("log")
axes[1].set_title("Logarithmic vertical scale")
figure.suptitle(rf"Periodic Ohmic kernel, $\beta={BETA}$ (equal-time term excluded)")
OUTPUT_DIRECTORY.mkdir(exist_ok=True)
for extension in ("png", "pdf"):
    figure.savefig(OUTPUT_DIRECTORY / f"dissipative-kernel-beta{BETA}.{extension}", dpi=FIGURE_DPI)
plt.close(figure)

print(f"D(1) = {kernel[0]:.8g}")
print(f"D(beta/2) = {kernel[BETA // 2 - MIN_SEPARATION]:.8g}")
print(f"Maximum reflection difference from finite image cutoff: {np.max(np.abs(kernel-kernel[::-1])):.3g}")
