import matplotlib.pyplot as plt
from scpp.tools.utils import save_or_show
import seaborn as sns


@save_or_show
def plot_dual_histograms(data1, data2, title1, title2, bins=100, figsize=(8, 4)):
    """
    Plot dual histograms with two subplots.

    Parameters:
    - data1: First set of data for the histogram.
    - data2: Second set of data for the histogram.
    - title1: Title for the first subplot.
    - title2: Title for the second subplot.
    - bins: Number of bins for the histograms. Defaults to 100.
    - figsize: Size of the figure. Defaults to (8, 4).

    Returns:
    None
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    p1 = sns.histplot(data1, bins=bins, kde=False, ax=axes[0], legend=False)
    axes[0].set_title(title1)

    p2 = sns.histplot(data2, bins=bins, kde=False, ax=axes[1], legend=False)
    axes[1].set_title(title2)
    axes[1].set_ylabel("")  # Hide the y-axis label for the second plot
