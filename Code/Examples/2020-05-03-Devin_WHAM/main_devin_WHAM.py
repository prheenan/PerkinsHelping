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
sys.path.append("../../../")
from Code.Lib.AppWHAM.Code import WeightedHistogram, UtilWHAM
from Code.Lib.AppWHAM.Lib.SimulationFEC import Test
import matplotlib.pyplot as plt


def output_fec_object_to_csv(f,fname):
    header = "Time\t,\tSeparation\t,\tZ\t,\tForce"
    Z_fwd = f.Offset + (f.Time - f.Time[0]) * f.Velocity
    X_fwd = np.array([f.Time, f.Separation, Z_fwd, f.Force])
    np.savetxt(fname,X=X_fwd.T, header=header)

def convert_csv_to_object(fname,k,kT):
    example = np.loadtxt(fname)
    time,extension,Z,force = \
        example[:,0], example[:,1], example[:,2], example[:,3]
    v = (Z[1]-Z[0])/abs(time[1]-time[0])
    offset = Z[0]
    to_ret = Test.SimpleFEC(Time=time,Extension=extension,Force=force,
                            kT=kT,SpringConstant=k,Offset=offset,
                            Velocity=v)
    return to_ret

def run():
    """
    """
    # generate some test data
    n_fecs = 10
    k = 0.1e-3
    kT = 4.1e-21
    fecs_fwd = list(Test.get_simulated_ensemble(n=n_fecs,k=k,
                                                reverse=True))
    fecs_rev = list(Test.get_simulated_ensemble(n=n_fecs,k=k))
    # save it out to simple CSVs
    for i, (f, r) in enumerate(zip(fecs_fwd, fecs_rev)):
        output_fec_object_to_csv(f,"fwd_{:0d}.csv".format(i))
        output_fec_object_to_csv(r,"rev_{:0d}.csv".format(i))
    # read the FECS back in
    fecs_fwd_reread = []
    fecs_rev_reread = []
    for i in range(n_fecs):
        fwd = convert_csv_to_object("fwd_{:d}.csv".format(i),k, kT)
        rev = convert_csv_to_object("rev_{:d}.csv".format(i),k,kT)
        fecs_fwd_reread.append(fwd)
        fecs_rev_reread.append(rev)
    fwd_wham = UtilWHAM.to_wham_input(fecs_fwd_reread)
    rev_wham = UtilWHAM.to_wham_input(fecs_rev_reread)
    landscape_WHAM = WeightedHistogram.wham(fwd_input=fwd_wham,
                                            rev_input=rev_wham)
    plt.close()
    fig = plt.figure()
    plt.plot(landscape_WHAM.q_nm, landscape_WHAM.G0_kT)
    plt.xlabel("Extension (nm)")
    plt.ylabel("$\Delta G_0$ ($k_\mathrm{B}T)$")
    fig.savefig("out_WHAM.png")
    pass

if __name__ == "__main__":
    run()
