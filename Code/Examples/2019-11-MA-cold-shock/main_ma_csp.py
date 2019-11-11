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

def read_pxps(input_dir):
    pxp_file_names, fecs = \
        FEC_Util.concatenate_fec_from_single_directory(input_dir)
    for f in fecs:
        yield f

def read_pxps_and_save_as_pkl(input_dir,cache_dir,force):
    load_function = lambda :  read_pxps(input_dir)
    name_func = FEC_Util.fec_name_func
    return CheckpointUtilities.multi_load(load_func=load_function,
                                          cache_dir=cache_dir,
                                          force=force,
                                          name_func=name_func)

def run():
    """

    """
    input_directory = "../../../Data/2019_MA/"
    reload_from_pxp = False
    cache_directory = "./cache/"
    fecs = read_pxps_and_save_as_pkl(input_dir=input_directory,
                                     cache_dir=cache_directory,
                                     force=reload_from_pxp)
    n = len(fecs)
    predictions = [_command_line_config.predict_indices(fec=ex,
                                                        threshold=5e-4,
                                                        tau_fraction=0.001)
                   for ex in fecs]
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
        PlotUtilities.lazyLabel(xlabel, "$F$ (pN)", "",useLegend=False)
        axs.append(ax)
        if i != 0:
            PlotUtilities.no_x_label(ax=ax)
    ys = np.concatenate([a.get_ylim() for a in axs])
    xs = np.concatenate([a.get_xlim() for a in axs])
    for a in axs:
        a.set_xlim([min(xs), max(xs)])
        a.set_ylim([min(ys), max(ys)])
    PlotUtilities.savefig(fig, "./out.png")


if __name__ == "__main__":
    run()
