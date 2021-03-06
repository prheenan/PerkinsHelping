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
from Lib.UtilGeneral import PlotUtilities
from Lib.UtilForce.FEC import FEC_Plot

from Lib.AppWLC.Code import WLC
from Lib.AppFEATHER.Code import Detector
import re
from Lib.UtilForce.UtilIgor import PxpLoader

from Lib.UtilForce.UtilIgor import TimeSepForceObj
from Lib.UtilGeneral import CheckpointUtilities

name_pattern = re.compile(r"""
                          (\D+) # non-digits (like 'Image')
                          (\d+) # digits (like '1010' in 'Image1010"
                          (\D+) # ends in force_ext or whatever 
                          """, re.VERBOSE)

def to_fec(x, F):
    meta = F.Note
    time = np.arange(F.DataY.size) * F.DeltaX()
    fec_appr = TimeSepForceObj._cols_to_TimeSepForceObj(meta_dict=meta,
                                                        force=F.DataY,
                                                        sep=x.DataY,
                                                        time=time)
    return fec_appr

def convert_data_to_fec(ex):
    # get the fecs for the components
    fec_appr = to_fec(ex["sep_ext'"], ex["force_ext'"])
    fec_dwell = to_fec(ex["sep_towd'"], ex["force_towd'"])
    fec_retr = to_fec(ex["sep_ret'"], ex["force_ret'"])
    # make sure the  times makes sense
    dt = fec_appr.Time[1] - fec_appr.Time[0]
    fec_dwell.Time += fec_appr.Time[-1] + dt
    fec_retr.Time += fec_dwell.Time[-1] + dt
    # combine all the data
    fec_final = fec_appr._slice(slice(0, None, 1))
    f_all = [fec_appr, fec_dwell, fec_retr]
    fec_final.Time = np.concatenate([f.Time for f in f_all])
    fec_final.ZSnsr = np.concatenate([f.ZSnsr for f in f_all])
    fec_final.Force = np.concatenate([f.Force for f in f_all])
    fec_final.Separation = np.concatenate([f.Separation for f in f_all])
    return fec_final

def fit_based_on_FEATHER(pred_info,split_fec):
    feather_surface = pred_info.slice_fit.start
    surface = max(0, feather_surface)
    starts = [surface] + [e.stop for e in pred_info.event_slices[:-1]]
    # figure out where each event ends
    ends = [e for e in pred_info.event_idx]
    slices = [slice(i, f) for i, f in zip(starts, ends)]
    slices = [s for s in slices if s.start < s.stop - 2]
    t_fec, x_fec, F_fec = split_fec.retract.Time, split_fec.retract.Separation, split_fec.retract.Force
    x_fec -= x_fec[surface]
    x_arr = [x_fec[s] for s in slices]
    F_arr = [F_fec[s] for s in slices]
    fit_common = dict(kbT=4.1e-21,Lp=0.4e-9, K0=1000e12)
    min_x = [ max(1e-9,min(x)/2) for x in x_arr ]
    max_x = [ max(x)*2 for x in x_arr ]
    N_max = 100
    L0_ranges = [ slice(x_min_tmp,x_max_tmp,(x_max_tmp-x_min_tmp)/N_max)
                  for x_min_tmp, x_max_tmp in zip(min_x,max_x)]
    fits = [WLC.fit(separation=x_tmp, force=F_tmp,
                    brute_dict=dict(Ns=N_max, ranges=[range_tmp]),
                    **fit_common)
            for x_tmp, F_tmp,range_tmp in zip(x_arr, F_arr,L0_ranges)]
    # returning: for each event: (contour length, rupture force, separation/Force)
    to_ret = [ [L0,F[-1],x,F] for (L0,F),x in zip(fits,x_arr)]
    return to_ret

def read_and_convert(pxp_file,use_split_names=False):
    if use_split_names:
        data = PxpLoader.LoadPxp(pxp_file, name_pattern=name_pattern,
                                 valid_ext_func=lambda _: True)
        # convert to my preferred format
        fecs = [ convert_data_to_fec(ex) for _,ex in data.items()]
    else:
        fecs =  FEC_Util.ReadInData(FullName=pxp_file)
    return fecs

def plot_split_fec(split_fec,ruptures,pred_info,fits,all_contour_changes):
    # make a plot of FEC (just the last one) + contour length
    # the last rupture isn't asociated with
    ruptures_per_dL = [r[:-1] for r in ruptures]
    X, F = split_fec.retract.Separation, split_fec.retract.Force
    plt.close()
    fig = PlotUtilities.figure(figsize=(3, 7))
    plt.subplot(4, 1, 1)
    plt.plot(X * 1e9, F * 1e12, color='k', alpha=0.3)
    for i in pred_info.event_idx:
        plt.axvline(X[i] * 1e9)
    for L0, _, x, F in fits:
        plt.plot(x * 1e9, F * 1e12, linestyle='--')
    PlotUtilities.lazyLabel("Extension (nm)", "$F$ (pN)", "")
    plt.subplot(4, 1, 2)
    plt.hist(np.array(all_contour_changes) * 1e9)
    PlotUtilities.lazyLabel("$\mathbf{\Delta}L_\mathbf{0}$ (nm)", "$N$ ", "")
    plt.subplot(4, 1, 3)
    plt.hist(np.concatenate(ruptures) * 1e12)
    PlotUtilities.lazyLabel("$F_R$ (pN)", "$N$ ", "")
    plt.subplot(4, 1, 4)
    plt.semilogy(np.concatenate(ruptures_per_dL) * 1e12,
                 np.array(all_contour_changes) * 1e9,
                 'r.')
    PlotUtilities.lazyLabel("$F_R$ (pN)",
                            "$\mathbf{\Delta}L_\mathbf{0}$ ", "")
    plt.xlim([0, 75])
    PlotUtilities.savefig(fig, "./nick.png")

def run():
    """

    """
    # read in the data and convert
    pxp_file = r"""./Newsies.pxp"""
    force_read_pxp = True
    fecs = CheckpointUtilities.getCheckpoint("./nick.pkl",
                                             read_and_convert,
                                             force_read_pxp, pxp_file,
                                             use_split_names=True)
    # understand after this line
    all_contour_changes = []
    ruptures = []
    fits_all = []
    if len(fecs) == 0:
        assert False , "Didn't find any FECs to analyze. Check your names. "
    # POST: should have found some FECs
    for i,fec in enumerate(fecs):
        print("Analyzing curve {:d} / {:d}".format(i,len(fecs)))
        split_fec, pred_info = \
            Detector._predict_full(fec, threshold=2e-2, tau_fraction=2e-3,
                                   f_refs=[Detector.delta_mask_function])
        # fit to the eventsfeather sees
        # returninng: L0, F_rupture, x_wlc, F_wlc
        fits = fit_based_on_FEATHER(pred_info, split_fec)
        # record thecontour length changes
        contour_lengths = np.concatenate([f[0] for f in fits])
        contour_changes = np.diff(contour_lengths)
        all_contour_changes.extend(contour_changes)
        ruptures.append([f[1] for f in fits])
        fits_all.append(fits)
        """
        # uncomment to plot 
        plt.close()
        ax_F = plt.subplot(3, 1, 1)
        X, F = split_fec.retract.Separation, split_fec.retract.Force
        T = split_fec.retract.Time
        F_limit = [-30, 100]
        plt.plot(T, F * 1e12)
        plt.plot(T, pred_info.interp(T) * 1e12)
        plt.ylim(F_limit)
        PlotUtilities.lazyLabel("", "$F$ (pN)", "")
        ax_P = plt.subplot(3, 1, 2)
        plt.semilogy(split_fec.retract.Time, pred_info.probabilities[-1])
        PlotUtilities.lazyLabel("Time (s)", "FEATHER Probability (au)", "")
        plt.subplot(3, 1, 3)
        plt.plot(X * 1e9, F * 1e12, color='k', alpha=0.3)
        PlotUtilities.no_x_label(ax=ax_F)
        for i in pred_info.event_idx:
            plt.axvline(X[i] * 1e9)
        for L0, F_rupture, x, F in fits:
            plt.plot(x * 1e9, F * 1e12, linestyle='--')
        plt.ylim(F_limit)
        PlotUtilities.lazyLabel("Separation (nm)", "$F$ (pN)", "")
        plt.show()
        """
    plot_split_fec(split_fec, ruptures, pred_info, fits,
                   all_contour_changes)



if __name__ == "__main__":
    run()
