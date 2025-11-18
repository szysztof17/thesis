import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl

import seaborn as sns
import numpy as np

figures = Path(__file__)

# === Default figure width ===
figure_width = 6.1



METHOD_PALETTE_INDICES = {
    "PIDC": 0,
    "GRNBOOST2": 2,
    "GENIE3": 1,
    # Add other methods here
    # "ANOTHER_METHOD": 3,
}

def get_method_color(method_name):
    """
    Gets the consistent color for a specific method.
    Assumes load_custom_style() has been called.
    """
    if method_name not in METHOD_PALETTE_INDICES:
        raise KeyError(f"Method '{method_name}' not found in METHOD_PALETTE_INDICES. "
                       f"Please add it to plot_utils.py.")
    
    palette_index = METHOD_PALETTE_INDICES[method_name]
    full_palette = get_hue_palette()
    
    if palette_index >= len(full_palette):
        raise IndexError(f"Index {palette_index} for method '{method_name}' is out of bounds "
                         f"for the loaded palette (size {len(full_palette)}).")
        
    return full_palette[palette_index]


# === Matplotlib style ===
def load_custom_style():
    sns.set_style("whitegrid")
    path = Path(__file__).parent.parent / "plt_params.mplstyle"
    plt.style.use(str(path))

    color_cycle = mpl.rcParams['axes.prop_cycle'].by_key()['color']

    sns.set_palette(color_cycle)


def load_custom_style(colours='default'):
    """
    Loads the mplstyle and enforces a specific color palette.
    palette_style: 'default', 'one', or 'two'
    """
    # 1. Reset to base style from file
    sns.set_style("whitegrid")
    path = Path(__file__).parent.parent / "plt_params.mplstyle"
    plt.style.use(str(path))

    # 2. If a specific palette is requested, overwrite the global cycle
    if colours != 'default':
        # Get the specific colors
        new_colors = get_hue_palette(style=colours)
        
        # Set Matplotlib global cycle
        mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=new_colors)
        
        # Set Seaborn global palette
        sns.set_palette(new_colors)
    else:
        # Just ensure seaborn matches the loaded mplstyle
        color_cycle = mpl.rcParams['axes.prop_cycle'].by_key()['color']
        sns.set_palette(color_cycle)



def get_hue_palette(n=None):
    """
    Returns the custom palette defined by the 'axes.prop_cycle' in the mplstyle.
    """
    palette = mpl.rcParams['axes.prop_cycle'].by_key()['color'] 
    if n:
        return palette[:n]
    else:
        return mpl.rcParams['axes.prop_cycle'].by_key()['color'] #

def get_hue_palette(n=None, style='default'):
    """
    Returns a palette based on the style argument.
    style: 'default' (uses rcParams), 'secondary', or 'diverging'
    """
    
    # Define your specific custom lists here
    custom_palettes = {
        'default': mpl.rcParams['axes.prop_cycle'].by_key()['color'],
        'one': ['#fe5f00', '#988f2a', '#6a5837', '#322f20'], 
        'two': ['#750d37', '#6d9f71', '#bc8034', '#d65f5f']      
    }

    # Retrieve the requested palette, fallback to default if not found
    palette = custom_palettes.get(style, custom_palettes['default'])

    if n:
        return palette[:n]
    else:
        return palette



# === Set figure width helper ===
def set_figure_width(width=6.1, ax=None, aspect_ratio=None, n_per_page=None):
    """Set figure width in inches and adjust height proportionally or via aspect ratio."""
    fig = plt.gcf() if ax is None else ax.figure
    if n_per_page:
        width  = width/n_per_page
    current_size = fig.get_size_inches()
    height = width / aspect_ratio if aspect_ratio else current_size[1] * (width / current_size[0])
    fig.set_size_inches(width, height)
    return fig

# === Custom color map ===
def get_custom_cmap():
    base = LinearSegmentedColormap.from_list("cmap", ["#ffffff", sns.color_palette()[0]])
    return LinearSegmentedColormap.from_list("trunc", base(np.linspace(0.05, 1, 256)))

    
def thousands_formatter(x, pos):
    return f'{int(x/1000)}K' if x >= 1000 else int(x)


def make_sequential_cmap(color, low=0.1, high=1, n=256):
    base = LinearSegmentedColormap.from_list("base", ["#ffffff", color])
    return LinearSegmentedColormap.from_list("trunc", base(np.linspace(low, high, n)))

def make_cmap_between(color1, color2, low=0.0, high=1.0, n=256, name="custom"):
    base = LinearSegmentedColormap.from_list(name, [color1, color2])
    return LinearSegmentedColormap.from_list(name + "_trunc", base(np.linspace(low, high, n)))
