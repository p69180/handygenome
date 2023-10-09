import colorsys
import collections

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


class LegendHandles(collections.UserList):
    def add_line(self, color='blue', marker='o', linewidth=0, markersize=8, label='label'):
        self.append(
            mpl.lines.Line2D(
                [], [], 
                color=color, 
                marker=marker, 
                markersize=markersize, 
                linewidth=linewidth,
                label=label,
            )
        )
    
    def add_patch(self):
        self.append(
            mpl.patches.Patch(color='red', label='The red data')
        )


def lightness_spectrum(color, n, l_start=0.3, l_end=0.7):
    color_hls = colorsys.rgb_to_hls(*mpl.colors.to_rgb(color))
    return [
        colorsys.hls_to_rgb(color_hls[0], x, color_hls[2]) 
        for x in np.linspace(l_start, l_end, n + 2, endpoint=True)[1:-1]
    ]


def get_yticklabel_size(n_yticks):
    return min((200 / n_yticks), 10)


def get_integer_yticks(ymin, ymax, num_ticks=15):
    ymin = np.rint(ymin)
    ymax = np.rint(ymax)
    step = int(max(1, np.rint((ymax - ymin) / num_ticks)))
    return np.arange(ymin, ymax, step).astype(int)


def draw_suptitle(fig, title, **kwargs):
    kwargs = (
        {'weight': 'bold', 'size': 20, 'y': 1} 
        | kwargs
    )
    fig.suptitle(title, **kwargs)


def get_boxplot_range(data):
    q25, q75 = np.nanquantile(data, [0.25, 0.75])
    iqr = q75 - q25
    offset = 1.5 * iqr
    ymax = q75 + offset
    ymin = q25 - offset
    return ymin, ymax


