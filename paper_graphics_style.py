import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_axes_aligner import align

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

    ax.plot(barx, bary, c='black',linewidth=rcParams['lines.linewidth']*.75)

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    ax.text(*mid, text, **kwargs)