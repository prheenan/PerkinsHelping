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

sys.path.append("../../")
from Lib.UtilForce.FEC import FEC_Util
# since we are importing the command line configs, we have to add the path to
# the main AppFEATHER location..
sys.path.append("../../Lib/AppFEATHER/")
from Lib.AppFEATHER.Code import _command_line_config
from Lib.UtilGeneral import PlotUtilities, CheckpointUtilities
from Lib.UtilForce.FEC import FEC_Plot
from Lib.AppFEATHER.Code import Detector
from Examples import UtilHelpPlot

def read_pxps(input_dir):
    """
    :param input_dir: directory with pxps
    :return: all the FECs in all the pxps in that directory
    """
    pxp_file_names, fecs = \
        FEC_Util.concatenate_fec_from_single_directory(input_dir)
    for f in fecs:
        yield f

def read_pxps_and_save_as_pkl(input_dir,cache_dir,force):
    """
    :param input_dir: input directory (full of pxps of beautiful AFM data)
    :param cache_dir: where to put the python-style 'pkl' data files
    :param force: if true, then force re-loading from the pxps
    :return: list of all the FECs
    """
    load_function = lambda :  read_pxps(input_dir)
    name_func = FEC_Util.fec_name_func
    return CheckpointUtilities.multi_load(load_func=load_function,
                                          cache_dir=cache_dir,
                                          force=force,
                                          name_func=name_func)

def demo_FEATHER_index_prediction(fecs,**kw_feather):
    """
    :param fecs: list of force-extension curves
    :return:  nothing, makes plot showing where the indices for each FEC are
    """
    n = len(fecs)
    # predict the indices where we have events
    predictions = [_command_line_config.predict_indices(fec=ex,**kw_feather)
                   for ex in fecs]
    # show all the data with the events overlayed
    plt.close()
    fig = PlotUtilities.figure((3, n * 2))
    axs = []
    for i, ex in enumerate(fecs):
        ax = plt.subplot(n, 1, (i + 1))
        FEC_Plot._fec_base_plot(ex.Time, ex.Force * -1e12)
        event_indices = predictions[i]
        for event_i in event_indices:
            plt.axvline(ex.Time[event_i])
        PlotUtilities.x_label_on_top(ax=ax)
        xlabel = "Time (s)" if i == 0 else ""
        PlotUtilities.lazyLabel(xlabel, "$F$ (pN)", "", useLegend=False)
        axs.append(ax)
        if i != 0:
            PlotUtilities.no_x_label(ax=ax)
    ys = np.concatenate([a.get_ylim() for a in axs])
    xs = np.concatenate([a.get_xlim() for a in axs])
    for a in axs:
        a.set_xlim([min(xs), max(xs)])
        a.set_ylim([min(ys), max(ys)])
    PlotUtilities.savefig(fig, "./FEATHER_index_prediction.png")

def demo_FEATHER_nuts_and_bolts(fecs,**kw_feather):
    """
    :param fecs: list of fecs
    :param kw_feather: input to Detector._predict_full (e.g. threshold, tau_f)
    :return:  nothing, makes plots showing how FEATHER works on each data point
    """
    for i,f in enumerate(fecs):
        name = FEC_Util.fec_name_func(i,f)
        # predict where the events are
        split_fec, pred_info = Detector._predict_full(f,**kw_feather)
        # make a helpful plot of the data
        save_name = "FEATHER_{:s}.png".format(name)
        F_limit_pN = [-20,550]
        UtilHelpPlot.make_FEATHER_plot(split_fec, pred_info,out_name=save_name,
                                       F_limit_pN = F_limit_pN)


def run():
    """

    """
    input_directory = "../../../Data/2019_MA/"
    reload_from_pxp = False
    cache_directory = "./cache/"
    fecs = read_pxps_and_save_as_pkl(input_dir=input_directory,
                                     cache_dir=cache_directory,
                                     force=reload_from_pxp)
    # set up the probability threshold (between 0 and 1)
    # and the smoothing fraction (between 0 and 1) for the data
    kw_feather = dict(threshold=1e-5,tau_fraction=0.003)
    # show the predictions FEATHER makes
    demo_FEATHER_index_prediction(fecs,**kw_feather)
    # show *how* FEATHER makes those predictions
    # (useful for debugging the parameters FEATHER uses)
    demo_FEATHER_nuts_and_bolts(fecs,**kw_feather)

if __name__ == "__main__":
    run()
