#!/usr/bin/env python3

import matplotlib.pyplot as plt
from plot_utils import load_custom_style, get_hue_palette

def preview_palette(n_colors=6, savepath="palette_preview.png"):
    # Load your custom style
    load_custom_style()

    # Get the palette
    palette = get_hue_palette(n_colors)

    # Plot swatches
    fig, ax = plt.subplots(figsize=(n_colors, 2))
    for i, color in enumerate(palette):
        ax.add_patch(plt.Rectangle((i, 0), 1, 1, color=color))
        ax.text(i + 0.5, -0.2, color, ha="center", va="top", fontsize=10)

    ax.set_xlim(0, len(palette))
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Save + show
    plt.savefig(savepath, dpi=150, bbox_inches="tight")
    plt.show()

if __name__ == "__main__":
    preview_palette()
