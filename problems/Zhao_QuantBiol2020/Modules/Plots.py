import os

import matplotlib.pyplot as plt


def plot_results(paths, results_dict):
    """Cumulative confirmed infected, model against data, one panel per (region, stage)."""
    plot_path, MODEL_NAME = paths["plot_path"], paths["MODEL_NAME"]
    n = len(results_dict)
    cols = 4
    fig, axes = plt.subplots((n + cols - 1) // cols, cols, figsize=(14, 3.4 * ((n + cols - 1) // cols)), squeeze=False)
    for ax, (label, item) in zip(axes.flat, results_dict.items()):
        results, data = item["results"], item["data"][label]
        ax.plot(results["time"], results["Confirmed_Infected_Population"], color="tab:red")
        ax.scatter(data["time"], data["Confirmed_Infected"], color="k", s=14, zorder=5)
        ax.set_title(label.replace("_", " "), fontsize=9)
        ax.set_xlabel("day")
    for ax in axes.flat[n:]:
        ax.axis("off")
    fig.supylabel("cumulative confirmed infected")
    tag = paths.get("plot_tag")
    fig.suptitle(MODEL_NAME + (f" - {tag}" if tag else ""))
    fig.tight_layout()
    plot_name = os.path.join(plot_path, MODEL_NAME + (f"_{tag}" if tag else "") + ".png")
    fig.savefig(plot_name, bbox_inches="tight")
    print(f"Plot saved to: {plot_name}")
    plt.close(fig)
