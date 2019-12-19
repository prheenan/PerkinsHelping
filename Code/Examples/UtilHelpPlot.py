# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from Lib.UtilGeneral import PlotUtilities

def make_FEATHER_plot(split_fec,pred_info,out_name,F_limit_pN=[-30,100]):
    """
    :param split_fec: first output from Detector._predict_full
    :param pred_info: second output from Detector._predict_full
    :param out_name: what to save the plot as
    :param F_limit_pN: y limits (useful if just looking for high or low events)
    :return:  nothing, makes plot
    """
    plt.close()
    fig = PlotUtilities.figure((3,7))
    ax_F = plt.subplot(3, 1, 1)
    X, F = split_fec.retract.Separation, split_fec.retract.Force
    T = split_fec.retract.Time
    F_limit = F_limit_pN
    plt.plot(T, F * 1e12,label="Raw data",color='k',alpha=0.3)
    plt.plot(T, pred_info.interp(T) * 1e12,color='r',
             label="Spline-smoothed")
    plt.ylim(F_limit)
    PlotUtilities.lazyLabel("", "$F$ (pN)", "",
                            legend_kwargs=dict(loc='upper left'))
    ax_P = plt.subplot(3, 1, 2)
    plt.semilogy(split_fec.retract.Time, pred_info.probabilities[-1])
    PlotUtilities.lazyLabel("Time (s)", "FEATHER Probability (au)", "")
    plt.subplot(3, 1, 3)
    plt.plot(X * 1e9, F * 1e12, color='k', alpha=0.3)
    PlotUtilities.no_x_label(ax=ax_F)
    for i in pred_info.event_idx:
        plt.axvline(X[i] * 1e9)
    plt.ylim(F_limit)
    PlotUtilities.lazyLabel("Separation (nm)", "$F$ (pN)", "")
    PlotUtilities.savefig(fig,out_name)
