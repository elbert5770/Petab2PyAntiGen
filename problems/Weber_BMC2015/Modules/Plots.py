import os

import matplotlib.pyplot as plt

from .Solver_settings import T_PREEQUIL


def plot_results(paths, results_dict):
    """The eight observables, model against data; the two conditions overlaid (time = hours after the switch)."""
    plot_path, MODEL_NAME = paths["plot_path"], paths["MODEL_NAME"]
    observables = list(dict.fromkeys(c for label, item in results_dict.items() for c in item["data"][label].columns
                                     if c != "time"))
    cols = 4
    rows = (len(observables) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(14, 3.2 * rows), squeeze=False)
    colors = dict(zip(results_dict, ("tab:blue", "tab:orange")))
    for ax, obs in zip(axes.flat, observables):
        for label, item in results_dict.items():
            data = item["data"][label]
            if obs not in data.columns:
                continue
            res = item["results"]
            ax.plot(res["time"] - T_PREEQUIL, res[obs], color=colors[label], label=label.replace("model1_", ""))
            ax.scatter(data["time"] - T_PREEQUIL, data[obs], color=colors[label], s=14, edgecolor="k", linewidth=0.4, zorder=5)
        ax.set_title(obs, fontsize=9)
        ax.set_xlabel("h after the switch")
        ax.legend(fontsize=7)
    for ax in axes.flat[len(observables):]:
        ax.axis("off")
    tag = paths.get("plot_tag")
    fig.suptitle(MODEL_NAME + (f" - {tag}" if tag else ""))
    fig.tight_layout()
    plot_name = os.path.join(plot_path, MODEL_NAME + (f"_{tag}" if tag else "") + ".png")
    fig.savefig(plot_name, bbox_inches="tight")
    print(f"Plot saved to: {plot_name}")
    plt.close(fig)
