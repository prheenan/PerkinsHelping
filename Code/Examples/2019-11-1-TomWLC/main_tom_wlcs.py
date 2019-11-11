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
from Lib.AppWLC.Code import TracerWLC


from Lib.UtilGeneral import PlotUtilities
from Lib.UtilGeneral.Plot import Annotations
from Lib.UtilGeneral.Plot import Scalebar


def run():
    """

    """
    p_nm = 50
    dL_nm = 0.1
    L_arr_nm = [500,1000, 1500, 2000]
    n_ex = 4
    kw_wlc = dict(Lp=p_nm, dL=dL_nm)
    np.random.seed(42)
    wlcs = [[TracerWLC.wlc(L0=L0_nm_tmp, **kw_wlc)
             for _ in range(n_ex)]
            for L0_nm_tmp in L_arr_nm]
    n_rows = len(L_arr_nm)
    n_cols = n_ex
    plt.close()
    axs = []
    fig = PlotUtilities.figure(figsize=(7, 7))
    colors = ['crimson', 'rebeccapurple', 'forestgreen', 'navy']
    example_i = 1
    size_max = 700
    for i, wlc_list in enumerate(wlcs):
        for j, wlc_tmp in enumerate(wlc_list):
            color = colors[i]
            ax_tmp = plt.subplot(n_rows, n_cols, example_i)
            x, y = wlc_tmp.x, wlc_tmp.y
            x0, y0 = x - np.mean(x), y - np.mean(y)
            plt.plot(x0, y0, color=color)
            example_i += 1
            axs.append(ax_tmp)
            PlotUtilities.color_frame(ax=ax_tmp, color=color)
            str_L0 = "$L_\mathbf{0}$ = " + \
                     "{:.1f} ".format(L_arr_nm[i] / 1000) + \
                     PlotUtilities.upright_mu() + "m"
            ax_tmp.set_xlim([-size_max, size_max])
            ax_tmp.set_ylim([-size_max, size_max])
            if j == 0:
                Annotations.relative_annotate(s=str_L0, color=color,
                                              xy=(0.4, 0.85), ax=ax_tmp,
                                              fontsize=12)
            font_x, _ = Scalebar.font_kwargs_modified(
                x_kwargs=dict(color=color))
            Scalebar.x_scale_bar_and_ticks_relative(unit="nm", width=500,
                                                    offset_x=0.75,
                                                    offset_y=0.05,
                                                    ax=ax_tmp,
                                                    line_kwargs=dict(
                                                        color=color,
                                                        linewidth=1.5),
                                                    font_kwargs=font_x)

    for a in axs:
        PlotUtilities.format_image_axis(ax=a, remove_frame=False)
    PlotUtilities.savefig(fig, "./WormLikeChainSamples.png")
    pass
    # read in the data and convert


if __name__ == "__main__":
    run()
