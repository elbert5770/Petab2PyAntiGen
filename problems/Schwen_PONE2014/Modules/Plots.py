import os

import matplotlib.pyplot as plt
import numpy as np

from .Data import load_conditions

PANELS = [("observable_IR1", "IR1 (cytometry)"), ("observable_IR2", "IR2 (cytometry)"), ("observable_IRsum", "IR sum (15 min, by dose)"),
          ("observable_Insulin_nExpID1", "insulin, series 1 (ELISA)"), ("observable_Insulin_nExpID2", "insulin, series 2"),
          ("observable_Insulin_nExpID3", "insulin, series 3"), ("observable_Insulin_nExpID4", "insulin, series 4")]


def plot_results(paths, results_dict):
    """Model (line) and data (dots) of every observable on its log10 scale, one panel each; a time course per condition, except the
    IR sum, which is one point per insulin dose."""
    plot_path, MODEL_NAME = paths["plot_path"], paths["MODEL_NAME"]
    cond = load_conditions().set_index("condition")
    fig, axes = plt.subplots(2, 4, figsize=(15, 7), squeeze=False)
    colors = plt.cm.viridis(np.linspace(0, 0.9, 6))
    for ax, (obs, title) in zip(axes.flat, PANELS):
        k = 0
        doses = []
        for label, item in results_dict.items():
            results, data = item["results"], item["data"][label]
            if obs not in data or data[obs].notna().sum() == 0:
                continue
            d = data[data[obs].notna()]
            if obs == "observable_IRsum":
                doses.append((cond.loc[label, "init_Ins"], float(np.interp(d["time"].iloc[0], results["time"], results[obs])),
                              float(d[obs].iloc[0])))
                continue
            color = colors[k % len(colors)]
            k += 1
            ax.plot(results["time"], results[obs], "-", color=color, lw=1.2,
                    label=f"{cond.loc[label, 'init_Ins']:g} nM")
            ax.scatter(d["time"], d[obs], color=color, s=12, zorder=5)
        if doses:
            doses.sort()
            ax.plot([max(x, 0.03) for x, _, _ in doses], [m for _, m, _ in doses], "o-", color="tab:red", mfc="none", ms=6)
            ax.scatter([max(x, 0.03) for x, _, _ in doses], [y for _, _, y in doses], color="k", s=14, zorder=5)
            ax.set_xscale("log")
            ax.set_xlabel("insulin (nM; 0 drawn at 0.03)", fontsize=7)
        else:
            ax.set_xlabel("time (min)", fontsize=7)
            ax.legend(fontsize=6, ncol=2)
        ax.set_title(title, fontsize=9)
        ax.set_ylabel("log10", fontsize=7)
    for ax in axes.flat[len(PANELS):]:
        ax.axis("off")
    tag = paths.get("plot_tag")
    fig.suptitle(MODEL_NAME + (f" - {tag}" if tag else "") + "   (line: model, dots: data)")
    fig.tight_layout()
    plot_name = os.path.join(plot_path, MODEL_NAME + (f"_{tag}" if tag else "") + ".png")
    fig.savefig(plot_name, bbox_inches="tight")
    print(f"Plot saved to: {plot_name}")
    plt.close(fig)
