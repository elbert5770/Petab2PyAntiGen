import os

import matplotlib.pyplot as plt
import numpy as np

from .Data import load_conditions

# the five experiments of Smith et al. (figure panels), by the suffix of the observable ids
FIGURES = {"2B": "dose response at 15 min: insulin dose (molecules)", "2C": "time course after 5e5 insulin: time (min)",
           "2D": "time course after 5e5 insulin: time (min)", "3B": "insulin pulse with and without H2O2 (condition)",
           "3C": "H2O2 outside the cell (molecules) at 960 min"}


def plot_results(paths, results_dict):
    """Model (line or open circle) and data (black) of every measured observable, one panel per figure of the paper."""
    plot_path, MODEL_NAME = paths["plot_path"], paths["MODEL_NAME"]
    cond = load_conditions().set_index("condition")
    panels = {}                                          # (figure, observable) -> list of (x, label, times, model, data)
    for label, item in results_dict.items():
        results, data = item["results"], item["data"][label]
        for obs in [c for c in data.columns if c != "time" and data[c].notna().any()]:
            d = data[data[obs].notna()]
            model = np.interp(d["time"], results["time"], results[obs])
            panels.setdefault((obs.split("__")[-1], obs), []).append((label, d["time"].to_numpy(), model, d[obs].to_numpy()))

    fig, axes = plt.subplots(3, 3, figsize=(12, 9), squeeze=False)
    for ax, ((fig_id, obs), points) in zip(axes.flat, sorted(panels.items())):
        if fig_id in ("2C", "2D"):                       # a time course: time on x
            for label, t, model, y in points:
                order = np.argsort(t)
                ax.plot(t[order], model[order], "o-", color="tab:red", mfc="none", ms=5, lw=1)
                ax.scatter(t, y, color="k", s=14, zorder=5)
        else:
            if fig_id == "2B":
                x = np.array([cond.loc[p[0], "Ins"] for p in points])
                ax.set_xscale("log")
            elif fig_id == "3C":
                x = np.array([cond.loc[p[0], "extracellular_ROS"] for p in points])
            else:
                x = np.arange(len(points))
                ax.set_xticks(x, [p[0].replace("figure3B__", "") for p in points], fontsize=6, rotation=45)
            order = np.argsort(x)
            xs = x[order]
            ax.plot(xs, [points[i][2][0] for i in order], "o", color="tab:red", mfc="none", ms=6)
            ax.scatter(xs, [points[i][3][0] for i in order], color="k", s=14, zorder=5)
        ax.set_title(obs.replace("__", " "), fontsize=9)
        ax.set_xlabel(FIGURES[fig_id], fontsize=7)
    for ax in axes.flat[len(panels):]:
        ax.axis("off")
    tag = paths.get("plot_tag")
    fig.suptitle(MODEL_NAME + (f" - {tag}" if tag else "") + "   (red: model, black: data)")
    fig.tight_layout()
    plot_name = os.path.join(plot_path, MODEL_NAME + (f"_{tag}" if tag else "") + ".png")
    fig.savefig(plot_name, bbox_inches="tight")
    print(f"Plot saved to: {plot_name}")
    plt.close(fig)
