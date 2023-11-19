import itertools

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


class IntervalSelector:
    @classmethod
    def from_ax(cls, ax):
        return cls(
            fig=ax.figure, 
            ax=ax, 
            omit_ax_labels=list(),
        )

    @classmethod
    def from_fig(cls, fig, omit_ax_labels=list()):
        return cls(
            fig=fig, 
            ax=None, 
            omit_ax_labels=omit_ax_labels,
        )

    def __init__(self, fig, ax, omit_ax_labels):
        self.fig = fig
        self.ax = ax
        self.omit_ax_labels = omit_ax_labels
        self.init_params()

    def get_axes(self):
        if self.ax is None:
            return [
                ax for ax in self.fig.get_axes() 
                if ax.get_label() not in self.omit_ax_labels
            ]
        else:
            return [self.ax]

    def init_params(self):
        '''axes is hashable'''
        self.cids = {ax: {'press': None, 'release': None} for ax in self.get_axes()}
        self.lines = {ax: list() for ax in self.get_axes()}
        self.boxes = {ax: list() for ax in self.get_axes()}

        self.events = {'press': list(), 'release': list()}
        self.coords = {'press': list(), 'release': list()}

        self.pressed = False
        self.on_press_ax = None

    def check_has_axes(self, ax):
        return ax in self.get_axes()

    def connect(self):
        for ax in self.get_axes():
            self.cids[ax]['press'] = self.fig.canvas.mpl_connect('button_press_event', self.on_press)
            self.cids[ax]['release'] = self.fig.canvas.mpl_connect('button_release_event', self.on_release)

    def disconnect(self):
        for subdic in self.cids.values():
            for cid in subdic.values():
                self.fig.canvas.mpl_disconnect(cid)
            for k in tuple(subdic.keys()):
                subdic[k] = None

    def on_press(self, event):
        if event.button == mpl.backend_bases.MouseButton.LEFT:
            self.on_press_left(event)
        elif event.button == mpl.backend_bases.MouseButton.RIGHT:
            self.on_press_right(event)
        else:
            return
        
    def on_press_left(self, event):
        if event.inaxes not in self.get_axes():
            return

        self.pressed = True
        self.on_press_ax = event.inaxes
        
        self.events['press'].append(event)
        self.coords['press'].append(event.xdata)

        for ax in self.get_axes():
            line = ax.axvline(event.xdata, color='tab:red', alpha=0.5, linewidth=2)
            linelist = list()
            linelist.append(line)
            self.lines[ax].append(linelist)
        
    def on_release(self, event):
        if not self.pressed:
            return
        if event.inaxes != self.on_press_ax:
            return

        self.pressed = False
        self.on_press_ax = None

        self.events['release'].append(event)
        self.coords['release'].append(event.xdata)

        for ax in self.get_axes():
            line = ax.axvline(event.xdata, color='tab:red', alpha=0.5, linewidth=2)
            linelist = self.lines[ax][-1]
            linelist.append(line)
            
            box = ax.fill_betweenx(
                ax.get_ylim(), 
                self.coords['press'][-1], 
                self.coords['release'][-1], 
                color='tab:red', 
                alpha=0.1,
            )
            self.boxes[ax].append(box)

    def on_press_right(self, event):
        if event.inaxes not in self.get_axes():
            return

        for ax in self.get_axes():
            for l in self.lines[ax][-1]:
                l.remove()
            del self.lines[ax][-1]
            self.boxes[ax][-1].remove()
            del self.boxes[ax][-1]

        for sublist in self.events.values():
            del sublist[-1]
        for sublist in self.coords.values():
            del sublist[-1]
        
    def get_intervals(self):
        arr = np.stack([self.coords['press'], self.coords['release']], axis=1)
        arr.sort(axis=1)
        return arr
    
    def reset(self):
        for linelist in itertools.chain.from_iterable(self.lines.values()):
            for x in linelist:
                x.remove()
        for b in itertools.chain.from_iterable(self.boxes.values()):
            b.remove()
        self.disconnect()
        self.init_params()

        
