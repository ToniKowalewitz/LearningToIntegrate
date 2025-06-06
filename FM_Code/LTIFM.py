from torchcfm.models import MLP
import os

import torch.nn as nn

import matplotlib.pyplot as plt
import torch
from torchdyn.core import NeuralODE
import pandas as pd
from tqdm import tqdm
from torchcfm.conditional_flow_matching import *
import json
from configargparse import Namespace
from torchcfm.models.models import *
from torchcfm.utils import *
from torch.optim.lr_scheduler import StepLR, ReduceLROnPlateau
from ltifm_toolkit import *
import seaborn as sns
device="cuda"

new_experiment=False


if new_experiment:
    opt=get_opts("A2_bigamma_OTCFM_ds100000_b50_e600_l0.001_w500_d4_sched_0.1_200.0")
    opt.batch_size=1000
    
    opt.step_size=opt.n_epochs/3
    opt.dataset_loc="Resources/modes_"+opt.modes+"_L_"+opt.distribution+"_N_"+str(opt.count)+".csv"
    num_modes=9 if "A2" in opt.dataset_loc else 25
    opt.parameters=["z"+str(i+1) for i in range(num_modes)]
    opt.name=opt.modes+"_"+opt.distribution+"_"+opt.fm_type+"_ds"+str(opt.datasize)+"_b"+str(opt.batch_size)+"_e"+str(opt.n_epochs)+"_l"+str(opt.lr)+"_w"+str(opt.width)+"_d"+str(opt.depth)
    if opt.scheduler==True:
        opt.name+="_sched_"+str(opt.gamma)+"_"+str(opt.step_size)



    os.makedirs(os.path.join("Experiments",opt.name),exist_ok=True)

    with open(os.path.join("Experiments",opt.name,"config.json"),"w") as f:
        json.dump(vars(opt),f)


    model,opt=fm_train_pipeline(opt)

else:
    opt=get_opts("A3_bigamma_CFM_ds100000_b1000_e600_l0.001_w500_d4_sched_0.1_200.0")

    model,opt=fm_train_pipeline(opt)
    fm_inference(opt,plots=True)
    validate_rp_all(opt)
