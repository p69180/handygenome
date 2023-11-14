import collections
import functools

import numpy as np
import pandas as pd
import scipy.stats
import scipy.signal
import scipy.linalg
import matplotlib as mpl
import matplotlib.pyplot as plt

#import handygenome.deco as deco
#import handygenome.logutils as logutils


class NoPeakError(Exception):
    pass


class DensityGenerationFailure(Exception):
    pass


class PeaksBase(collections.UserDict):
    @property
    def npeaks(self):
        return len(self['peak_xs'])

    def get_height_cutoff(self):
        if 'height' in self['find_peaks_kwargs'].keys():
            return self['find_peaks_kwargs']['height'][0]
        else:
            return None

    def subset(self, selector):
        """Values of 'peak_indexes' are not changed."""
        srcdict = dict()

        # first level value
        for key in (
            'peak_indexes',
            'peak_xs',
            'peak_ys',
            'ndata',
        ):
            if key in self:
                srcdict[key] = self[key][selector]

        # second level
        for key in (
            'peak_properties',
            'xcoord_values',
        ):
            if key in self:
                srcdict[key] = dict()
                for subkey, subval in self[key].items():
                    srcdict[key][subkey] = subval[selector]

        # do not subset
        for key in (
            'data',
            'find_peaks_kwargs',
            'xs',
            'ys',
            'density',
            'bin_edges',
            'bin_midpoints',
            'hist',
        ):
            if key in self:
                srcdict[key] = self[key]

        result = self.__class__(srcdict)
        return result

    @classmethod
    def merge_two(cls, left, right):
        srcdict = dict()

        # first level value
        for key in (
            #'peak_indexes',
            'peak_xs',
            'peak_ys',
            'ndata',
        ):
            srcdict[key] = list()
            for peakobj in (left, right):
                if key in peakobj:
                    srcdict[key].extend(peakobj[key])

        # second level
        for key in (
            'peak_properties',
            'xcoord_values',
        ):
            srcdict[key] = dict()
            for peakobj in (left, right):
                if key in peakobj:
                    for subkey, subval in peakobj[key].items():
                        srcdict[key].setdefault(subkey, list())
                        srcdict[key][subkey].extend(subval)

        result = cls(srcdict)
        return result

    @classmethod
    def merge(cls, *args):
        peakobj_list = list(args)
        srcdict = dict()

        # first level value
        for key in (
            #'peak_indexes',
            'peak_xs',
            'peak_ys',
            'ndata',
        ):
            srcdict[key] = list()
            for peakobj in peakobj_list:
                if key in peakobj:
                    srcdict[key].extend(peakobj[key])
            srcdict[key] = np.array(srcdict[key])

        # second level
        for key in (
            'peak_properties',
            'xcoord_values',
        ):
            srcdict[key] = dict()
            for peakobj in peakobj_list:
                if key in peakobj:
                    for subkey, subval in peakobj[key].items():
                        srcdict[key].setdefault(subkey, list())
                        srcdict[key][subkey].extend(subval)
            for subkey in tuple(srcdict[key].keys()):
                srcdict[key][subkey] = np.array(srcdict[key][subkey])

        result = cls(srcdict)
        return result

    #def get_highest_peak_index(self):
    #    max_peak_idx = np.argmax(self['peak_properties']['peak_heights'])
    #    return self['peak_indexes'][max_peak_idx]

    #############
    # filtering #
    #############

    def remove_contained_peaks(self):
        peak_properties = self['peak_properties']
        margins = [
            (x, y) for (x, y) in
            zip(
                peak_properties['left_ips'],
                peak_properties['right_ips'],
            )
        ]

        idxs_to_retain = list()
        for idx, query_margin in enumerate(margins):
            target_iterator = (
                x[1] for x in enumerate(margins) if x[0] != idx
            )
            contained = any(
                (
                    (query_margin[0] >= target_margin[0])
                    and (query_margin[1] <= target_margin[1])
                )
                for target_margin in target_iterator
            )
            if not contained:
                idxs_to_retain.append(idx)

        return self.subset(idxs_to_retain)

    def filter_ndata(self, n):
        self.subset(self['ndata'] >= n)

    def plot_peakinfo(
        self, ax, peakline_kwargs, draw_ndata=False, 
        draw_peak_positions=True,
        annotate_kwargs=dict(),
    ):
        height_cutoff = self.get_height_cutoff()
        if height_cutoff is not None: 
            ax.axhline(height_cutoff, color='green', lw=2, alpha=0.4)

        for x, height in zip(
            self['peak_xs'], 
            self['peak_properties']['peak_heights'],
        ):
            ax.axvline(x, **peakline_kwargs)

        if draw_ndata:
            ylim = ax.get_ylim()
            y = ylim[1] + 0.05 * (ylim[1] - ylim[0])
            for (x, ndata) in zip(self['peak_xs'], self['ndata']):
                ax.text(x, y, str(ndata), va='center', ha='left', rotation=90)

        if draw_peak_positions:
            annotate_kwargs = (
                dict(
                    va='bottom', ha='center', size=8,
                    arrowprops=dict(arrowstyle='-', color='black', linewidth=0.5),
                )
                | annotate_kwargs
            )
            xtext_list = np.linspace(*ax.get_xlim(), self.npeaks, endpoint=True)
            for x, xtext in zip(self['peak_xs'], xtext_list):
                #ax.axvline(x, color='black', linewidth=0.5)
                ax.annotate(
                    round(x, 3), 
                    #(ax.get_xlim()[1], y),
                    (x, ax.get_ylim()[1]),
                    xytext=(xtext, ax.get_ylim()[1] * 1.1),
                    annotation_clip=False,
                    **annotate_kwargs,
                )

    def plot_peakinfo_lying(
        self, ax, peakline_kwargs, draw_ndata=False,
        draw_peak_positions=True,
        annotate_kwargs=dict(),
    ):
        height_cutoff = self.get_height_cutoff()
        if height_cutoff is not None: 
            ax.axvline(height_cutoff, color='green', lw=2, alpha=0.4)

        for y in self['peak_xs']:
            ax.axhline(y, **peakline_kwargs)

        if draw_ndata:
            xlim = ax.get_xlim()
            x = xlim[1] + 0.05 * (xlim[1] - xlim[0])
            for (y, ndata) in zip(self['peak_xs'], self['ndata']):
                ax.text(x, y, str(ndata), va='center', ha='left', rotation=0)

        if draw_peak_positions:
            annotate_kwargs = (
                dict(
                    va='center', ha='left', size=8,
                    arrowprops=dict(arrowstyle='-', color='black', linewidth=0.5),
                )
                | annotate_kwargs
            )
            ytext_list = np.linspace(*ax.get_ylim(), self.npeaks, endpoint=True)
            xlim = ax.get_xlim()
            for y, ytext in zip(self['peak_xs'], ytext_list):
                #ax.axhline(y, color='black', linewidth=0.5)
                ax.annotate(
                    round(y, 3), 
                    (xlim[1], y),
                    xytext=(xlim[1] * 1.1, ytext),
                    annotation_clip=False,
                    **annotate_kwargs,
                )

    def postprocess(self, reltoheight=None):
        if reltoheight is not None:
            self.postprocess_width_reltoheight(reltoheight)
        self.postprocess_to_xcoords()
        self.postprocess_ndata()

    def postprocess_ndata(self):
        left_borders = self['xcoord_values']['left_ips']
        right_borders = self['xcoord_values']['right_ips']

        gt_left = self['data'][:, np.newaxis] >= left_borders
        lt_right = self['data'][:, np.newaxis] <= right_borders

        self['ndata'] = np.logical_and(gt_left, lt_right).sum(axis=0)

    def postprocess_to_xcoords(self):
        def convert_point(peak_properties_key):
            xs_indexes = self['peak_properties'][peak_properties_key]
            xs = self['xs']
            portions = xs_indexes / len(xs)
            return xs[0] + (xs[-1] - xs[0]) * portions

        def convert_length(peak_properties_key):
            lengths = self['peak_properties'][peak_properties_key]
            xs = self['xs']
            portions = lengths / len(xs)
            return (xs[-1] - xs[0]) * portions

        xcoord_values = dict()
        for key in (
            'left_bases', 
            'right_bases',

            'left_ips', 
            'right_ips',

            'left_ips_reltoheight', 
            'right_ips_reltoheight',
        ):
            if key in self['peak_properties']:
                xcoord_values[key] = convert_point(key)

        for key in (
            'widths',
            'widths_reltoheight',
        ):
            if key in self['peak_properties']:
                xcoord_values[key] = convert_length(key)

        self['xcoord_values'] = xcoord_values

    def postprocess_width_reltoheight(self, reltoheight):
        peak_properties = self['peak_properties']
        prominence_data = (
            peak_properties['prominences'],
            peak_properties['left_bases'],
            peak_properties['right_bases'],
        )
        rel_heights = (
            (peak_properties['peak_heights'] * reltoheight)
            / peak_properties['prominences']
        )
        width_results = list()
        for idx, (x_index, rel) in enumerate(
            zip(self['peak_indexes'], rel_heights)
        ):
            partial_result = scipy.signal.peak_widths(
                self['ys'], 
                [x_index], 
                rel_height=rel,
                prominence_data=(
                    np.atleast_1d(peak_properties['prominences'][idx]),
                    np.atleast_1d(peak_properties['left_bases'][idx]),
                    np.atleast_1d(peak_properties['right_bases'][idx]),
                ),
            )
            width_results.append(partial_result)

        peak_properties['widths_reltoheight'] = np.concatenate(
            [x[0] for x in width_results]
        )
        peak_properties['width_heights_reltoheight'] = np.concatenate(
            [x[1] for x in width_results]
        )
        peak_properties['left_ips_reltoheight'] = np.concatenate(
            [x[2] for x in width_results]
        )
        peak_properties['right_ips_reltoheight'] = np.concatenate(
            [x[3] for x in width_results]
        )

    def undecal(self, x):
        '''peak_indexes are not modified'''
        def undecal_value(arr):
            return (2 * x) - arr

        def modify_pos(arr, gt_selector):
            arr[gt_selector] = undecal_value(arr[gt_selector])

        def modify_borders(lt_arr, rt_arr, gt_selector):
            undecaled_rt_values = undecal_value(rt_arr[gt_selector])
            undecaled_lt_values = undecal_value(lt_arr[gt_selector])
            lt_arr[gt_selector] = undecaled_rt_values
            rt_arr[gt_selector] = undecaled_lt_values

        gt_selector = self['peak_xs'] > x

        # peak_xs
        modify_pos(self['peak_xs'], gt_selector)

        # peak_properties, xcoord_values
        peak_properties = self['peak_properties']
        xcoord_values = self['xcoord_values']
        for lt_key, rt_key in (
            ('left_bases', 'right_bases'),
            ('left_ips', 'right_ips'),
            ('left_ips_reltoheight', 'right_ips_reltoheight'),
        ):
            if lt_key in peak_properties:
                modify_borders(
                    peak_properties[lt_key], 
                    peak_properties[rt_key], 
                    gt_selector,
                )
                modify_borders(
                    xcoord_values[lt_key], 
                    xcoord_values[rt_key], 
                    gt_selector,
                )


class HistPeaks(PeaksBase):
    def plot(
        self, 
        figsize=None, 
        ax=None, 
        peakline_kwargs=dict(),
        bar_kwargs=dict(),
        omit_peakline=False,
    ):
        peakline_kwargs = (
            {'lw': 0.5, 'color': 'blue'}
            | peakline_kwargs
        )
        bar_kwargs = (
            {'fill': 'tab:blue'}
            | bar_kwargs
        )

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = None

        x = self['bin_edges'][:-1]
        height = self['hist']
        width = np.diff(self['bin_edges'])
        ax.bar(x, height, width, align='edge', **bar_kwargs)

        if not omit_peakline:
            self.plot_peakinfo(ax, peakline_kwargs)

        return fig, ax

    def plot_lying(
        self,
        figsize=None, 
        ax=None, 
        peakline_kwargs=dict(),
        bar_kwargs=dict(),
        omit_peakline=False,
    ):
        peakline_kwargs = (
            {'lw': 0.5, 'color': 'blue'}
            | peakline_kwargs
        )
        bar_kwargs = (
            {'fill': 'tab:blue'}
            | bar_kwargs
        )

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = None

        y = self['bin_edges'][:-1]
        width = self['hist']
        height = np.diff(self['bin_edges'])
        ax.barh(y, width, height=height, align='edge', **bar_kwargs)

        if not omit_peakline:
            self.plot_peakinfo_lying(ax, peakline_kwargs)

        return fig, ax


class DensityPeaks(PeaksBase):
    def plot(
        self, 
        figsize=None, 
        ax=None, 
        peakline_kwargs=dict(),
        plot_kwargs=dict(),
        omit_peakline=False,
    ):
        peakline_kwargs = (
            {'lw': 0.5, 'color': 'red'}
            | peakline_kwargs
        )
        plot_kwargs = (
            {'color': 'tab:red'}
            | plot_kwargs
        )

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = None

        ax.plot(self['xs'], self['ys'], **plot_kwargs)

        if not omit_peakline:
            self.plot_peakinfo(ax, peakline_kwargs)

        return fig, ax

    def plot_lying(
        self,
        figsize=None, 
        ax=None, 
        peakline_kwargs=dict(),
        plot_kwargs=dict(),
        omit_peakline=False,
    ):
        peakline_kwargs = (
            {'lw': 0.5, 'color': 'red'}
            | peakline_kwargs
        )
        plot_kwargs = (
            {'color': 'tab:red'}
            | plot_kwargs
        )

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = None

        ax.plot(self['ys'], self['xs'], **plot_kwargs)
        if not omit_peakline:
            self.plot_peakinfo_lying(ax, peakline_kwargs)

        return fig, ax


def find_hist_peaks(
    data, 

    plotonly=False,

    weights=None, 
    limit=None, 
    reltoheight=None,
    inverse=False,

    as_density=True,

    bins=10, 
    **find_peaks_kwargs,
):
    data, weights, find_peaks_kwargs = findpeaks_arghandler(
        data, weights, limit, find_peaks_kwargs,
    )

    hist, bin_edges = np.histogram(
        data, bins=bins, weights=weights, density=as_density,
    )
    bin_midpoints = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    if not plotonly:
        if inverse:
            findpeaks_input = np.max(hist) - hist
        else:
            findpeaks_input = hist
        peaks, peak_properties = scipy.signal.find_peaks(
            findpeaks_input, **find_peaks_kwargs
        )
        if len(peaks) == 0:
            raise NoPeakError('Peaks could not be found from histogram')

        peakresult = HistPeaks(
            {
                'data': data,
                'find_peaks_kwargs': find_peaks_kwargs,
                'peak_indexes': peaks,
                'peak_xs': bin_midpoints[peaks],
                'peak_ys': peak_properties['peak_heights'],
                'peak_properties': peak_properties,
                'bin_edges': bin_edges,
                'bin_midpoints': bin_midpoints,
                'hist': hist,
                'xs': bin_midpoints,
                'ys': hist,
            }
        )
        peakresult.postprocess(reltoheight=reltoheight)
    else:
        peakresult = HistPeaks(
            {
                'data': data,
                'find_peaks_kwargs': find_peaks_kwargs,
                'bin_edges': bin_edges,
                'bin_midpoints': bin_midpoints,
                'hist': hist,
                'xs': bin_midpoints,
                'ys': hist,
            }
        )

    return peakresult


def find_density_peaks(
    data, 

    plotonly=False,

    weights=None, 
    limit=None, 
    reltoheight=None,
    inverse=False,

    xs=100, 
    bw_method=None, 
    **find_peaks_kwargs,
):
    data, weights, find_peaks_kwargs = findpeaks_arghandler(
        data, weights, limit, find_peaks_kwargs,
    )

    try:
        density = scipy.stats.gaussian_kde(
            data, weights=weights, bw_method=bw_method,
        )
    except Exception as exc:
        msg = f'data: {data}; weights: {weights}'
        raise DensityGenerationFailure(msg) from exc

    if isinstance(xs, int):
        xs = np.linspace(data.min(), data.max(), xs)
        #xs = np.linspace(np.quantile(data, 0.01), np.quantile(data, 0.99), xs)
    ys = density(xs)

    if not plotonly:
        if inverse:
            findpeaks_input = np.max(ys) - ys
        else:
            findpeaks_input = ys
        peaks, peak_properties = scipy.signal.find_peaks(
            findpeaks_input, **find_peaks_kwargs
        )

        if len(peaks) == 0:
            raise NoPeakError('Peaks could not be found from density')

        peakresult = DensityPeaks(
            {
                'data': data,
                'find_peaks_kwargs': find_peaks_kwargs,
                'peak_indexes': peaks,
                'peak_xs': xs[peaks],
                'peak_ys': peak_properties['peak_heights'],
                'peak_properties': peak_properties,
                'density': density,
                'xs': xs,
                'ys': ys,
            }
        )
        peakresult.postprocess(reltoheight=reltoheight)
    else:
        peakresult = DensityPeaks(
            {
                'data': data,
                'find_peaks_kwargs': find_peaks_kwargs,
                'density': density,
                'xs': xs,
                'ys': ys,
            }
        )

    return peakresult


def findpeaks_arghandler(data, weights, limit, find_peaks_kwargs):
    data = np.asarray(data)
    if data.ndim != 1:
        raise Exception(f'Input data must be 1-d.')
    if len(data) == 0:
        raise Exception(f'Input data length is 0')

    if weights is None:
        weights = np.repeat(1, len(data))
    else:
        weights = np.asarray(weights)

    if limit is not None:
        selector = np.logical_and(data > limit[0], data < limit[1])
        data = data[selector]
        weights = weights[selector]

    # find_peaks_kwargs
    find_peaks_kwargs = (
        {
            'height': (None, None), 
            'width': (None, None),
            'rel_height': 0.5,
        }
        | find_peaks_kwargs
    )

    return data, weights, find_peaks_kwargs


def draw_peaks(
    data=None, 

    ax=None,

    weights=None, 
    limit=None, 
    reltoheight=None,
    inverse=False,

    xs=100,
    bins=30,

    hist=False,

    omit_hist=None,
    omit_hist_peakline=None,
    omit_density=None,
    omit_density_peakline=None,

    histpeaks=None,
    densitypeaks=None,

    hist_kwargs=dict(),
    density_kwargs=dict(), 

    lying=False,
):
    # argument handling
    hist_kwargs = (
        dict()
        | hist_kwargs
        | dict(bins=bins)
    )
    density_kwargs = (
        dict()
        | density_kwargs
        | dict(xs=xs)
    )

    for dic in (hist_kwargs, density_kwargs):
        dic['weights'] = weights
        dic['limit'] = limit
        dic['reltoheight'] = reltoheight
        dic['inverse'] = inverse

    def helper(old, new):
        if old is None:
            return new
        else:
            return old

    if hist:
        omit_hist = helper(omit_hist, False)
        omit_hist_peakline = helper(omit_hist_peakline, False)
        omit_density = helper(omit_density, True)
        omit_density_peakline = helper(omit_density_peakline, True)
    else:
        omit_hist = helper(omit_hist, False)
        omit_hist_peakline = helper(omit_hist_peakline, True)
        omit_density = helper(omit_density, False)
        omit_density_peakline = helper(omit_density_peakline, False)

    # main
    if (histpeaks is None) and (not omit_hist):
        histpeaks = find_hist_peaks(
            data, 
            plotonly=omit_hist_peakline,
            **hist_kwargs,
        )
    if (densitypeaks is None) and (not omit_density):
        try:
            densitypeaks = find_density_peaks(
                data, 
                plotonly=omit_density_peakline,
                **density_kwargs,
            )
        except DensityGenerationFailure:
            omit_density = True

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    if not omit_hist:
        if lying:
            histpeaks.plot_lying(
                ax=ax, 
                omit_peakline=omit_hist_peakline,
            )
        else:
            histpeaks.plot(
                ax=ax, 
                omit_peakline=omit_hist_peakline,
            )
    if not omit_density:
        if lying:
            densitypeaks.plot_lying(
                ax=ax, 
                omit_peakline=omit_density_peakline,
            )
        else:
            densitypeaks.plot(
                ax=ax, 
                omit_peakline=omit_density_peakline,
            )

    return fig, ax, densitypeaks, histpeaks


