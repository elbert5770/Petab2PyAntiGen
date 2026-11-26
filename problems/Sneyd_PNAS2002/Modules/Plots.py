import os

import matplotlib.pyplot as plt


def plot_results(paths, results_dict):
    """Open probability against time, model and data, one panel per condition."""
    plot_path, MODEL_NAME = paths["plot_path"], paths["MODEL_NAME"]
    n = len(results_dict)
    cols = 3
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(11, 3.0 * rows), squeeze=False, sharey=True)
    for ax, (label, item) in zip(axes.flat, results_dict.items()):
        results, data = item["results"], item["data"][label]
        ax.plot(results["time"], results["open_probability"], color="tab:red")
        ax.scatter(data["time"], data["open_probability"], color="k", s=12, zorder=5)
        ax.set_title(label.replace("_dose_response__", " dose response ").replace("Ca ", "Ca2+ "), fontsize=9)
        ax.set_xlabel("time (s)")
    for ax in axes.flat[n:]:
        ax.axis("off")
    fig.supylabel("open probability")
    tag = paths.get("plot_tag")
    fig.suptitle(MODEL_NAME + (f" - {tag}" if tag else ""))
    fig.tight_layout()
    plot_name = os.path.join(plot_path, MODEL_NAME + (f"_{tag}" if tag else "") + ".png")
    fig.savefig(plot_name, bbox_inches="tight")
    print(f"Plot saved to: {plot_name}")
    plt.close(fig)
