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
    return np.arange(ymin, ymax + step, step).astype(int)


def draw_suptitle(fig, title, **kwargs):
    kwargs = (
        {'weight': 'bold', 'size': 20, 'y': 0.97} 
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


def qqplot(data1, data2, q=np.arange(0, 1.01, 0.01), ax=None):
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    quantiles1 = np.quantile(data1, q)
    quantiles2 = np.quantile(data2, q)
    ax.plot(quantiles1, quantiles2, linestyle='', marker='o', markersize=4, zorder=1)
    ax.set_title('QQ plot')

    new_min = min(ax.get_xlim()[0], ax.get_ylim()[0])
    new_max = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.set_xlim(new_min, new_max)
    ax.set_ylim(new_min, new_max)
    ax.set_aspect('equal')
    ax.axline((new_min, new_min), slope=1, color='black', zorder=0, linewidth=0.5)

    return fig, ax


def clear_ax(ax):
    ax.clear()
    ax.spines[:].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


def data_into_point(ax):
    x0, x1, y0, y1 = 0, 1, 0, 1
    disp_x0, disp_y0 = ax.transData.transform((x0, y0))
    disp_x1, disp_y1 = ax.transData.transform((x1, y1))
    inch_x0, inch_y0 = ax.figure.dpi_scale_trans.inverted().transform((disp_x0, disp_y0))
    inch_x1, inch_y1 = ax.figure.dpi_scale_trans.inverted().transform((disp_x1, disp_y1))

    inch_x = inch_x1 - inch_x0
    inch_y = inch_y1 - inch_y0
    point_x = inch_x * 72
    point_y = inch_y * 72

    return inch_x, inch_y, point_x, point_y


