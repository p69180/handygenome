import colorsys
import collections

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


class LegendHandles(collections.UserList):
    #def add_line(self, color='blue', marker='o', linewidth=0, markersize=8, label=None, alpha=1):
    def add_line(self, label=None, **kwargs):
        kwargs = (
            dict(
                color='blue', marker='o', linewidth=0, markersize=8, alpha=1,
            ) 
            | dict(label=label)
            | kwargs
        )
        self.append(
            mpl.lines.Line2D([], [], **kwargs)
        )
    
    def add_patch(self, color='red', label=None):
        self.append(
            mpl.patches.Patch(color=color, label=label)
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


def vstacked_barplot(
    data, 
    colors=None, 
    xlabels=None, 
    legend_labels=None, 
    ax=None, 
    xs=None, 
    bar_width=None, 
    rotation=0, 
    align='center',
    zorder=3,
    dont_set_xticks=False,
):
    """Args:
        data: 
            - ndim == 2 
            - axis 0 means different kinds of data, axis 1 means samples
                (same spatial arrangement as plot itself)
            - upper row (rows with smaller index) comes on top in the figure
    """
    data = np.asarray(data)
    assert data.ndim == 2

    # arg handling
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    if xlabels is None:
        xlabels = np.arange(data.shape[1])
    else:
        assert len(xlabels) == data.shape[1]

    if legend_labels is None:
        legend_labels = np.repeat(None, data.shape[0])
    else:
        assert len(legend_labels) == data.shape[0]

    if colors is None:
        colors = mpl.cm.viridis(np.linspace(0, 1, data.shape[0], endpoint=False))
    else:
        assert len(colors) == data.shape[0]

    # main
    bottoms = np.concatenate(
        [
            np.repeat(0, data.shape[1])[np.newaxis, :],
            np.cumsum(data[:-1, :], axis=0),
        ],
        axis=0,
    )
    if xs is None:
        xs = np.arange(data.shape[1])

    bars = list()
    for val, bot, col, leglabel in zip(data, bottoms, colors, legend_labels):
        kwargs = dict(
            x=xs, 
            height=val, 
            bottom=bot, 
            facecolor=col, 
            label=leglabel, 
            align=align,
            zorder=zorder,
        )
        if bar_width is not None:
            kwargs['width'] = bar_width
        bar = ax.bar(**kwargs)
        textxs = xs
        textys = np.asarray(bot) + 0.5 * np.asarray(val)
        texts = np.asarray(val)
        for x, y, t in zip(textxs, textys, texts):
            if t != 0:
                ax.text(x, y, str(t), ha='center', va='center')

        
        bars.append(bar)

    if not dont_set_xticks:
        ax.set_xticks(xs, labels=xlabels, rotation=rotation)

    if not all((x is None) for x in legend_labels):
        ax.legend()

    return fig, ax, bars


stacked_barplot = vstacked_barplot


def hstacked_barplot(data, group_width, xlabels, legend_labels=None, ax=None):
    data = np.asarray(data)
    assert data.ndim == 2

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    bar_width = group_width / data.shape[0]
    xs_common = np.arange(data.shape[1])
    xs_offsets = -0.5 * group_width + np.linspace(0, group_width, data.shape[0], endpoint=False)

    if legend_labels is None:
        legend_labels = np.repeat(None, data.shape[0])
    else:
        assert len(legend_labels) == data.shape[0]

    bars = list()
    for subdata, offset, label in zip(data, xs_offsets, legend_labels):
        xs = xs_common + offset
        bar = ax.bar(xs, subdata, width=bar_width, label=label, align='edge')
        bars.append(bar)

    ax.set_xticks(xs_common, labels=xlabels)
    if not all((x is None) for x in legend_labels):
        ax.legend()

    return fig, ax, bars



