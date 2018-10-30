# -*- coding: utf-8 -*-
""" Plotting tools for concentration measurement assessment. """
import numpy as np
import matplotlib.pyplot as plt
import windtunnel as wt


__all__ = [
    'plot_boxplots',
]

def plot_boxplots(data_dict, ylabel=None, **kwargs):
    """ Plot statistics of concentration measurements in boxplots. Expects
    input from PointConcentration class.
    @parameters: data, type = dict
    @parameters: ylabel, type = string
    @parameter ax: axis passed to function
    @parameter kwargs : additional keyword arguments passed to plt.boxplot()
    """
    # Set standard ylabel if none is specified
    if ylabel is None:
        ylabel = 'Concentration'

    # Generate namelist from dict keys
    namelist = list(data_dict.keys())
    data = [np.nan for i in range(len(namelist))]
    maxes = []
    for i, key in enumerate(namelist):
        data[i] = data_dict[key]['concentration']
        maxes.append(np.max(data[i]))

    numDists = len(data)
    fig, ax1 = plt.subplots(figsize=(10, 6))
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

    bp = ax1.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')

    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)

    # Hide these grid behind plot objects
    ax1.set_axisbelow(True)
    ax1.set_title('Your concentration measurements')
    ax1.set_xlabel('Measurement')
    ax1.set_ylabel(ylabel)

    # Now fill the boxes with desired colors
    boxColors = ['darkkhaki', 'royalblue']
    numBoxes = numDists
    medians = list(range(numBoxes))
    for i in range(numBoxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = np.column_stack([boxX, boxY])
        # Alternate between Dark Khaki and Royal Blue
        k = i % 2
        boxPolygon = patches.Polygon(boxCoords, facecolor=boxColors[k])
        ax1.add_patch(boxPolygon)
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            ax1.plot(medianX, medianY, 'k')
            medians[i] = medianY[0]
        # Finally, overplot the sample averages, with horizontal alignment
        # in the center of each box
        ax1.plot([np.average(med.get_xdata())], [np.average(data[i])],
                 color='w', marker='*', markeredgecolor='k')

    # Set the axes ranges and axes labels
    ax1.set_xlim(0.5, numBoxes + 0.5)
    top = np.max(maxes) + 0.1 * np.max(maxes)
    # if np.max(maxes) > 400:
    bottom = -10
    ax1.set_ylim(bottom, top)
    ax1.set_xticklabels(namelist, rotation=45, fontsize=8)

    # Due to the Y-axis scale being different across samples, it can be
    # hard to compare differences in medians across the samples. Add upper
    # X-axis tick labels with the sample medians to aid in comparison
    # (just use two decimal places of precision)
    pos = np.arange(numBoxes) + 1
    upperLabels = [str(np.round(s, 2)) for s in medians]
    weights = ['bold', 'semibold']
    for tick, label in zip(range(numBoxes), ax1.get_xticklabels()):
        k = tick % 2
        ax1.text(pos[tick], top - (top * 0.05), upperLabels[tick],
                 horizontalalignment='center', size='medium', weight=weights[k],
                 color=boxColors[k])