import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_axes_aligner import align
from wesanderson import wes_palettes
from matplotlib.ticker import FuncFormatter
import numpy as np

sns.set_context('paper')
rcParams['savefig.dpi'] = 900

def p_convert(x):
    if x < .001: return '***'
    elif x < .01: return '**'
    elif x < .05: return '*'
    elif x < .06: return '~'
    else: return ''

def mm2inch(*tupl):
    inch = 25.4
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
def paired_barplot_annotate_brackets(txt, x_tick, height, y_lim, dh=.05, barh=.05, fs=10, xtick_spread=.05, maxasterix=None, ax=None):
    """ 
    Annotate barplot with p-values.

    :param txt: string to write or number for generating asterixes
    :param x_tick: center of pair of bars
    :param height: heights of the errors in question
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    if type(txt) is str:
        text = txt
    else:
        # * is p < 0.05
        # ** is p < 0.005
        # *** is p < 0.0005
        # etc.
        text = ''
        p = .05

        while txt < p:
            text += '*'
            p /= 10.

            if maxasterix and len(text) == maxasterix:
                break

        if len(text) == 0:
            text = 'n. s.'

    ly, ry = height[0], height[1]

    if type(x_tick) == int:
        lx, rx = x_tick-xtick_spread, x_tick+xtick_spread
    else:
        lx, rx = x_tick[0], x_tick[1]
    # lx, ly = x_tick-.05, height[0]
    # rx, ry = x_tick+.05, height[1]

    ax_y0, ax_y1 = y_lim
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)

    y = max(ly, ry) + dh

    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y)

    ax.plot(barx, bary, c='black',linewidth=rcParams['lines.linewidth']*.75,zorder=100)

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    ax.text(*mid, text, **kwargs,zorder=100)

def align_yaxis(ax1, v1, ax2, v2, y2min, y2max):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1."""

    """where y2max is the maximum value in your secondary plot. I haven't
     had a problem with minimum values being cut, so haven't set this. This
     approach doesn't necessarily make for axis limits at nice near units,
     but does optimist plot space"""

    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    scale = 1
    while scale*(maxy+dy) < y2max:
        scale += 0.05
    ax2.set_ylim(scale*(miny+dy), scale*(maxy+dy))

def no_leading_zeros(x,pos):
    """Format 1 as 1, 0 as 0, and all values whose absolute values is between
    0 and 1 without the leading "0." (e.g., 0.7 is formatted as .7 and -0.4 is
    formatted as -.4)."""
    val_str = '{:g}'.format(x)
    if x > 0 and x < 1:
        # return val_str.replace("0", "", 1)
        val_str = val_str[1:]
        if len(val_str) == 2:
            val_str += '0'
        return val_str
    if x < 0 and x > -1:
        val_str = val_str[0] + val_str[2:]
        if len(val_str) == 3:
            val_str += '0'
        return val_str
    else:
        return val_str