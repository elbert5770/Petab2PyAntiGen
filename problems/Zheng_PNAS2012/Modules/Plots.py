import os

import matplotlib.pyplot as plt

from antimony_modules.Histone.methylation_Zheng_PNAS2012 import STATES, state
from .Solver_settings import T_PREEQUIL


def plot_results(paths, results_dict):
    """The 15 H3 states (fraction of total H3) against the measurements."""
    plot_path, MODEL_NAME = paths["plot_path"], paths["MODEL_NAME"]
    fig, axes = plt.subplots(4, 4, figsize=(14, 11), sharex=True)
    for label, item in results_dict.items():
        results, data = item["results"], item["data"][label]
        t = results["time"] - T_PREEQUIL
        for ax, (a, b) in zip(axes.flat, STATES):
            name = f"observable_{state(a, b)}"
            ax.plot(t, results[name], color="tab:red")
            ax.scatter(data["time"] - T_PREEQUIL, data[name], color="k", s=14, zorder=5)
            ax.set_title(state(a, b), fontsize=9)
    for ax in axes.flat[len(STATES):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel("Time [h]")
    fig.supylabel("fraction of total H3")
    tag = paths.get("plot_tag")
    fig.suptitle(MODEL_NAME + (f" - {tag}" if tag else ""))
    fig.tight_layout()
    plot_name = os.path.join(plot_path, MODEL_NAME + (f"_{tag}" if tag else "") + ".png")
    fig.savefig(plot_name, bbox_inches="tight")
    print(f"Plot saved to: {plot_name}")
    plt.close(fig)
