from pymatgen.util.testing import PymatgenTest
from atomate2.abinit.flows.gw import G0W0ConvergenceMaker
from atomate2.abinit.jobs.core import StaticMaker, NonSCFMaker
from atomate2.abinit.jobs.gw import ScreeningMaker, SigmaMaker
from jobflow import run_locally
from atomate2.abinit.files import load_abinit_input
from abipy.electrons.gw import SigresFile
import os

struct = PymatgenTest.get_structure("Si")


scf = StaticMaker.from_params(ecut=12.0, kppa=10, spin_mode="unpolarized", nband=4)
nscf = NonSCFMaker.from_params(nband=100, kppa=10, extra_abivars={"nbdbuf": 20}, mode="uniform", shifts="Gamma")
scr1 = ScreeningMaker.from_params(ecuteps=1.2, nband=10)
# scr2 = ScreeningMaker.from_params(ecuteps=1.2, nband=20)
# scr3 = ScreeningMaker.from_params(ecuteps=1.2, nband=30)
# scr1 = ScreeningMaker.from_params(ecuteps=4.0, nband=180)
# sigmas = []
# for ecuteps in [0.8, 1.6, 2.4, 3.2, 4.0]:
#     for nband in [20, 40, 60, 80, 100, 120, 140, 160, 180]:
#         sigmas.append(SigmaMaker.from_params(ecuteps=ecuteps, nband=nband, ecutsigx=12.0))
sigma1 = SigmaMaker.from_params(ecuteps=0.4, ecutsigx=8, nband=10)
sigma2 = SigmaMaker.from_params(ecuteps=0.4, ecutsigx=8, nband=20)
sigma3 = SigmaMaker.from_params(ecuteps=0.4, ecutsigx=8, nband=30)
sigma4 = SigmaMaker.from_params(ecuteps=0.4, ecutsigx=8, nband=40)
sigma5 = SigmaMaker.from_params(ecuteps=0.8, ecutsigx=8, nband=10)
sigma6 = SigmaMaker.from_params(ecuteps=0.8, ecutsigx=8, nband=20)
sigma7 = SigmaMaker.from_params(ecuteps=0.8, ecutsigx=8, nband=30)
sigma8 = SigmaMaker.from_params(ecuteps=0.8, ecutsigx=8, nband=40)
sigma9 = SigmaMaker.from_params(ecuteps=1.2, ecutsigx=8, nband=10)
sigma10 = SigmaMaker.from_params(ecuteps=1.2, ecutsigx=8, nband=20)
sigma11 = SigmaMaker.from_params(ecuteps=1.2, ecutsigx=8, nband=30)
sigma12 = SigmaMaker.from_params(ecuteps=1.2, ecutsigx=8, nband=40)
maker = G0W0ConvergenceMaker(scf_maker=scf,
                             nscf_maker=nscf,
                             # scr_makers=[scr1, scr2, scr3],
                             sigma_makers=[sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7, sigma8, sigma9, sigma10, sigma11, sigma12],
                             # sigma_makers=sigmas,
                             scr_makers=[scr1],
                             # sigma_makers=[sigma1]
                             )

import json
with open('gwconv_maker.json', 'w') as f:
    json.dump(maker.as_dict(), f)

flow = maker.make(structure=struct)

flow.draw_graph().show()

out = run_locally(flow, create_folders=True)


# print(out.keys())
sigma_dirs = []
for k in out.keys():
    # print(k, out[k].keys())
    for idx in out[k].keys():
        # print(out[k][idx].output.dict().keys())
        # print(out[k][idx].output.task_label)
        if "Sigma" in out[k][idx].output.task_label:
            sigma_dirs.append(out[k][idx].output.dir_name)
# print(sigma_dirs)

sigma_gamma_gap = {}
for sigma_dir in sigma_dirs:
    ai = load_abinit_input(sigma_dir)
    sigres = SigresFile(os.path.join(sigma_dir, "outdata", "out_SIGRES.nc"))
    qpgap_gamma = sigres.qpgaps.flatten()[0]
    ecuteps = ai["ecuteps"]
    nband = ai["nband"]
    if ecuteps not in sigma_gamma_gap:
        sigma_gamma_gap[ecuteps] = []
    sigma_gamma_gap[ecuteps].append((nband, qpgap_gamma, sigma_dir))


import matplotlib.pyplot as plt

fig, subplot = plt.subplots()
for ecuteps, nbands_gammagaps in sigma_gamma_gap.items():
    nbg = sorted(nbands_gammagaps, key=lambda x: x[0])
    nb = [item[0] for item in nbg]
    gp = [item[1] for item in nbg]
    subplot.plot(nb, gp, "o-", label=f"ecuteps {ecuteps}")

# print(sigma_gamma_gap)
subplot.legend()
plt.show()
