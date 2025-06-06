import pandas as pd
from ltifm_toolkit import validate_rp_all
import numpy as np
import math as m
import os
import json
import multiprocessing
from configargparse import Namespace

df=pd.DataFrame(index=["A2","A3"],columns=["gamma","bigamma","norm","pois"])
df2=pd.DataFrame(index=["A2","A3"],columns=["gamma","bigamma","norm","pois"])
path="/homes/math/pkrueger/work/FlowMatching/Learning_to_Integrate/RP_Bootstrap_Res_Arrays"

for dist in ["gamma","bigamma","norm","pois"]:
    for i,modes in enumerate(["A2","A3"]):
        with open (os.path.join(path,modes+"_"+dist+".json")) as f:
            print(modes,dist)
            list=json.load(f)
            print(len(list))
            arr=np.array(list)
            df.iloc[i][dist]=np.percentile(arr,90)
            df2.iloc[i][dist]=np.percentile(arr,95)
print(df)
print(df2)

df.to_csv("CI_90_monomials.csv")
df2.to_csv("CI_95_monomials.csv")