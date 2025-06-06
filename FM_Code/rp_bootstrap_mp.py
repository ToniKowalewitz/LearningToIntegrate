import pandas as pd
from ltifm_toolkit import validate_rp_all
import numpy as np
import math as m
import os
import json
import multiprocessing
from configargparse import Namespace




def worker(result,start,end,modes,dist):
    
    path="Resources/modes_"+modes+"_L_"+dist+"_N_100000.csv"

    if "A2" in path:
        parameters=["z1", "z2", "z3", "z4", "z5", "z6", "z7", "z8", "z9"]
    else:
        parameters= ["z1", "z2", "z3", "z4", "z5", "z6", "z7", "z8", "z9", "z10", "z11", "z12", "z13", "z14", "z15", "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23", "z24", "z25"]



    df_mcmc=pd.read_csv(path)
    mcmc_arr=df_mcmc[parameters].to_numpy()

    
    for k in range(start,end):
 
        res=0
        df_bootstrap=df_mcmc.sample(100000,replace=True,random_state=k+opt.start)
        bootstrap_arr=df_bootstrap[parameters].to_numpy()



        for i in range(len(parameters)):
            counts = np.zeros(len(parameters))
            counts[i]=1
            print(counts)
            def func(x,counts=counts):
                return np.prod(x**counts)

            bootstrap_points=np.apply_along_axis(func,arr=bootstrap_arr,axis=1)
            mcmc_points=np.apply_along_axis(func,axis=1,arr=mcmc_arr)
            res+=np.abs(bootstrap_points.mean()-mcmc_points.mean())

        for i in range(len(parameters)):
            for j in range(i,len(parameters)):
                counts = np.zeros(len(parameters))
                counts[i]+=1
                counts[j]+=1
                print(counts)
                def func(x,counts=counts):
                    return np.prod(x**counts)
                bootstrap_points=np.apply_along_axis(func,arr=bootstrap_arr,axis=1)
                mcmc_points=np.apply_along_axis(func,axis=1,arr=mcmc_arr)
                res+=np.abs(bootstrap_points.mean()-mcmc_points.mean())

        for n in range(len(parameters)):   
            for i in range(n,len(parameters)):
                for j in range(i,len(parameters)):
                    counts = np.zeros(len(parameters))
                    counts[i]+=1
                    counts[j]+=1
                    counts[n]+=1
                    print(counts)
                    def func(x,counts=counts):
                        return np.prod(x**counts)
                    bootstrap_points=np.apply_along_axis(func,arr=bootstrap_arr,axis=1)
                    mcmc_points=np.apply_along_axis(func,axis=1,arr=mcmc_arr)
                    res+=np.abs(bootstrap_points.mean()-mcmc_points.mean())
        result[k]=(res/(m.factorial(len(parameters)+2)/(m.factorial(len(parameters)-1)*m.factorial(3)) + m.factorial(len(parameters)+1)/(m.factorial(len(parameters)-1)*m.factorial(2))+m.factorial(len(parameters))/(m.factorial(len(parameters)-1)*m.factorial(1))))
        print(result[k])
    
def main(opt):


    result=multiprocessing.Array("d",opt.num_tests)
    segment=opt.num_tests//opt.core_count
    processes=[]
    
    for q in range(opt.core_count):
        start = q * segment
        if q == opt.core_count - 1:
            end = len(result)  # Ensure the last segment goes up to the end
        else:
            end = start + segment
        # Creating a process for each segment
        p = multiprocessing.Process(target=worker, args=(result, start, end, opt.modes,opt.dist))
        processes.append(p)
        p.start()

    for p in processes:
        p.join()
    with open(os.path.join("/homes/math/pkrueger/work/FlowMatching/Learning_to_Integrate/RP_Bootstrap_Res_Arrays",opt.modes+"_"+opt.dist+"_1.json"),"w") as f:
        json.dump(list(result),f)
    return result
   
if __name__ == '__main__':
    opt=Namespace()
    opt.num_tests=50
    opt.start=50

    opt.dist="bigamma"
    opt.modes="A3"

    opt.core_count=multiprocessing.cpu_count()
    print(opt.modes)
    print(opt.dist)
    print("start ",opt.start)
    print(opt.core_count)
    result = main(opt)
   
    print(list(result))
     
