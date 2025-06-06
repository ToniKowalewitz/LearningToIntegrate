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
import seaborn as sns
import math as m
device="cuda"


def get_opts(experiment_name):
    with open(os.path.join("Experiments",experiment_name,"config.json")) as f:
        dict=json.load(f)
    opt=Namespace(**dict)
    return opt



def seed(opt):
    torch.manual_seed(opt.seed)
    torch.cuda.manual_seed(opt.seed)
    np.random.seed(opt.seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    os.environ["PYTHONHASHSEED"] = str(opt.seed)

def plot_groundtruth(opt):
    
    df=pd.read_csv(opt.dataset_loc)
    fig=sns.pairplot(df[opt.parameters],kind="scatter",corner=True)
    plt.savefig(os.path.join("Groundtruths",opt.dataset_loc.split("/")[1]+".png"))
    plt.close()

def make_dataloader(opt):
    seed(opt)
    df=pd.read_csv(opt.dataset_loc)[:opt.datasize]
  
    pos=df[opt.parameters].to_numpy()
    
    pos=torch.Tensor(pos)


    train_loader = torch.utils.data.DataLoader(pos,
        batch_size=opt.batch_size, shuffle=True, drop_last=True)

    
   
    return train_loader



class MLP(torch.nn.Module):

    #Define the neural network that is trained to approximate the target vector field.
    
    
    def __init__(self, dim, out_dim=None, w=64,d=5, activation="SELU",time_varying=False):
        super().__init__()
        self.time_varying = time_varying
        if out_dim is None:
            out_dim = dim
        self.dict_act={"SELU":torch.nn.SELU(),"Softplus":torch.nn.Softplus(),"ELU":torch.nn.ELU()}
        self.layers=[torch.nn.Linear(dim + (1 if time_varying else 0), w),self.dict_act[activation]]
        for i in range(d):
            self.layers.append(torch.nn.Linear(w,w))
            self.layers.append(self.dict_act[activation])
        self.layers.append(torch.nn.Linear(w,out_dim))

        self.net = torch.nn.Sequential(*self.layers)

    def forward(self, x):
        return self.net(x)


def fm_train_pipeline(opt):
    
    #Set global seed for repeatability of experiments
    seed(opt)
    
    #Initialize the Dataloader

    train_loader=make_dataloader(opt)

    #Initialiue Flow Matching loss and sample functions
    if opt.fm_type=="CFM":
        FM=ConditionalFlowMatcher(sigma=opt.sigma)
    if opt.fm_type=="OTCFM":
        FM=ExactOptimalTransportConditionalFlowMatcher(sigma=opt.sigma)
    if opt.fm_type=="SBCFM":
        FM=SchrodingerBridgeConditionalFlowMatcher(sigma=opt.sigma)
    #Initialize model and set to device
    model=MLP(dim=len(opt.parameters),out_dim=len(opt.parameters),activation=opt.activation,time_varying=True,w=opt.width,d=opt.depth)
    model=model.to(device)
        
    
    #Initialize optimizer
    optimizer=torch.optim.AdamW(model.parameters(),lr=opt.lr)
    #scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer=optimizer,factor=0.25,patience=100,verbose=True)
    scheduler=torch.optim.lr_scheduler.StepLR(optimizer=optimizer,gamma=opt.gamma,step_size=opt.step_size,verbose=True)
    #Training loop
    for epoch in tqdm(range(opt.n_epochs)):
        losses=np.array([0],dtype=np.float64)
        batch_idx=0
        for x in train_loader:
            batch_idx+=1
            optimizer.zero_grad()
            x1 = x.to(device)
            
            x0 = torch.randn_like(x1)
            t, xt, ut = FM.sample_location_and_conditional_flow(x0, x1)
            vt = model(torch.cat([xt, t[:, None]], dim=-1))
            loss = torch.mean((vt - ut) ** 2)
            losses[0]+=loss
            loss.backward()
            optimizer.step()
        #scheduler.step(loss)
        scheduler.step()
        with open(os.path.join("Experiments",opt.name,"log.txt"),"a") as fh:
            fh.write(f"Loss= {losses/batch_idx}\n")
    #Save model weights
    torch.save(model.state_dict(),os.path.join("Experiments",opt.name,"model.pth"))
    

    return model,opt

def fm_inference(opt,plots=True):

    seed(opt)
    FM=ConditionalFlowMatcher(sigma=opt.sigma)
    model=MLP(dim=len(opt.parameters),out_dim=len(opt.parameters),time_varying=True,w=opt.width,d=opt.depth,activation=opt.activation)
    model=model.to(device)
    model.load_state_dict(torch.load(os.path.join("Experiments",opt.name,"model.pth")))
    
    fun=lambda t, x,args: model(torch.cat([x, t.repeat(x.shape[0])[:, None],], dim=-1))
    node = NeuralODE(fun, solver="dopri5", sensitivity="adjoint", atol=1e-4, rtol=1e-4)


    size=10000

    traj = node.trajectory(
            torch.randn(size,len(opt.parameters), device=device),
            t_span=torch.linspace(0, 1, 2, device=device),
            
        )

    x_gen=traj[1,:,:][:]
    x_gen=x_gen.cpu().detach().numpy()
    df=pd.DataFrame()
    df[opt.parameters]=x_gen
    
    if "Inf" in plots or plots=="All":
        fig=sns.pairplot(df[opt.parameters],kind="scatter",corner=True)
        plt.savefig(os.path.join("Experiments",opt.name,"Inf.png"))
        plt.close()
    
    if "A3" in opt.name:
        df=pd.read_csv("/homes/math/pkrueger/work/FlowMatching/Learning_to_Integrate/Resources/sparse_grids_d_25.csv")
    else:
        df=pd.read_csv("/homes/math/pkrueger/work/FlowMatching/Learning_to_Integrate/Resources/sparse_grids_d_9.csv")

    grid=torch.Tensor(df[opt.parameters].to_numpy())
    grid=grid.to(device)

    traj_sparse = node.trajectory(
            grid,
            t_span=torch.linspace(0, 1, 2, device=device),
            
        )

    x_gen_sparse=traj_sparse[1,:,:][:]
    x_gen_sparse=x_gen_sparse.cpu().detach().numpy()
    #df=pd.DataFrame()
    df[opt.parameters]=x_gen_sparse
    df.to_csv(os.path.join("Experiments",opt.name,"Sparse_Inference_"+opt.name+".csv"))
    
    if "Sparse" in plots or plots=="All":
        fig=sns.pairplot(df[opt.parameters],kind="scatter",corner=True)
        plt.savefig(os.path.join("Experiments",opt.name,"Sparse.png"))
        plt.close()
    
 
    


    if "A2" in opt.name:
        df=pd.read_csv("/homes/math/pkrueger/work/FlowMatching/Learning_to_Integrate/Resources/linegrid.csv")
        grid_line=torch.Tensor(df[opt.parameters].to_numpy())
        grid_line=grid_line.to(device)

        traj_line = node.trajectory(
                grid_line,
                t_span=torch.linspace(0, 1, 2, device=device),
                
            )

        x_gen_line=traj_line[1,:,:][:]
        x_gen_line=x_gen_line.cpu().detach().numpy()
        #df=pd.DataFrame()
        df[opt.parameters]=x_gen_line
        df.to_csv(os.path.join("Experiments",opt.name,"Line_Inference_"+opt.name+".csv"))
        if "Line" in plots or plots=="All":
            fig=sns.pairplot(df[opt.parameters],kind="scatter",corner=True)
            plt.savefig(os.path.join("Experiments",opt.name,"Line.png"))
            plt.close()    
    
    return x_gen



def make_rp(dim,max_deg):
    if max_deg==0:
        q=max_deg
    else:
        q=np.random.choice(np.arange(1,max_deg+1,1),1)
        print(q)
    vec=np.arange(1,dim+1,1) #samplen hieraus

    samples=np.random.choice(vec,q)
    counts = [np.sum(samples == element) for element in vec]
    print(counts)
    def func(x,counts=counts):
        return np.prod(x**counts)

    return func

def validate_rp_all(opt,name="val_rp"):

    df_sparse=pd.read_csv(os.path.join("Experiments",opt.name,"Sparse_Inference_"+opt.name+".csv"))
    df_mcmc=pd.read_csv(os.path.join("Resources","modes_"+opt.modes+"_L_"+opt.distribution+"_N_100000.csv"))
    res_array=np.array([0.,0.,0.,0.])
    
    for k in range(4):
        print("K=",k,"\n","-----------------------------------")
        if k==0:
            res_array[k]=0
        
        elif k==1:
            for i in range(len(opt.parameters)):
                counts = np.zeros(len(opt.parameters))
                counts[i]=1
                print(counts)
                def func(x,counts=counts):
                    return np.prod(x**counts)
                sparse_points=df_sparse[df_sparse["k"]==k+1][opt.parameters].to_numpy()
                weights=df_sparse[df_sparse["k"]==k+1]["weight"].to_numpy()
                sparse_points=np.apply_along_axis(func,axis=1,arr=sparse_points)
                sparse_points*=weights
                mcmc_points=np.apply_along_axis(func,axis=1,arr=df_mcmc[opt.parameters].to_numpy())
                res_array[k]+=np.abs(sparse_points.sum()-mcmc_points.mean())
            #print(m.factorial(len(opt.parameters)+k-1)/m.factorial(len(opt.parameters)-1)*m.factorial(k))
            res_array[k]/=m.factorial(len(opt.parameters)+k-1)/m.factorial(len(opt.parameters)-1)*m.factorial(k)
        
        
        
        
        elif k==2:

            for i in range(len(opt.parameters)):
                counts = np.zeros(len(opt.parameters))
                counts[i]=1
                print(counts)
                def func(x,counts=counts):
                    return np.prod(x**counts)
                sparse_points=df_sparse[df_sparse["k"]==k+1][opt.parameters].to_numpy()
                weights=df_sparse[df_sparse["k"]==k+1]["weight"].to_numpy()
                sparse_points=np.apply_along_axis(func,axis=1,arr=sparse_points)
                sparse_points*=weights
                mcmc_points=np.apply_along_axis(func,axis=1,arr=df_mcmc[opt.parameters].to_numpy())
                res_array[k]+=np.abs(sparse_points.sum()-mcmc_points.mean())

            for i in range(len(opt.parameters)):
                for j in range(i,len(opt.parameters)):
                    counts = np.zeros(len(opt.parameters))
                    counts[i]+=1
                    counts[j]+=1
                    print(counts)
                    def func(x,counts=counts):
                        return np.prod(x**counts)
                    sparse_points=df_sparse[df_sparse["k"]==k+1][opt.parameters].to_numpy()
                    weights=df_sparse[df_sparse["k"]==k+1]["weight"].to_numpy()
                    sparse_points=np.apply_along_axis(func,axis=1,arr=sparse_points)
                    sparse_points*=weights
                    mcmc_points=np.apply_along_axis(func,axis=1,arr=df_mcmc[opt.parameters].to_numpy())
                    res_array[k]+=np.abs(sparse_points.sum()-mcmc_points.mean())
            #print((m.factorial(len(opt.parameters)+k-1)/(m.factorial(len(opt.parameters)-1)*m.factorial(k)) + m.factorial(len(opt.parameters)+k-2)/(m.factorial(len(opt.parameters)-1)*m.factorial(k-1))))
            res_array[k]/=(m.factorial(len(opt.parameters)+k-1)/(m.factorial(len(opt.parameters)-1)*m.factorial(k)) + m.factorial(len(opt.parameters)+k-2)/(m.factorial(len(opt.parameters)-1)*m.factorial(k-1)))
        
        
        elif k==3:
            for i in range(len(opt.parameters)):
                counts = np.zeros(len(opt.parameters))
                counts[i]=1
                print(counts)
                def func(x,counts=counts):
                    return np.prod(x**counts)
                sparse_points=df_sparse[df_sparse["k"]==k+1][opt.parameters].to_numpy()
                weights=df_sparse[df_sparse["k"]==k+1]["weight"].to_numpy()
                sparse_points=np.apply_along_axis(func,axis=1,arr=sparse_points)
                sparse_points*=weights
                mcmc_points=np.apply_along_axis(func,axis=1,arr=df_mcmc[opt.parameters].to_numpy())
                res_array[k]+=np.abs(sparse_points.sum()-mcmc_points.mean())
            
            for i in range(len(opt.parameters)):
                for j in range(i,len(opt.parameters)):
                    counts = np.zeros(len(opt.parameters))
                    counts[i]+=1
                    counts[j]+=1
                    print(counts)
                    def func(x,counts=counts):
                        return np.prod(x**counts)
                    sparse_points=df_sparse[df_sparse["k"]==k+1][opt.parameters].to_numpy()
                    weights=df_sparse[df_sparse["k"]==k+1]["weight"].to_numpy()
                    sparse_points=np.apply_along_axis(func,axis=1,arr=sparse_points)
                    sparse_points*=weights
                    mcmc_points=np.apply_along_axis(func,axis=1,arr=df_mcmc[opt.parameters].to_numpy())
                    res_array[k]+=np.abs(sparse_points.sum()-mcmc_points.mean())
            
            for n in range(len(opt.parameters)):   
                for i in range(n,len(opt.parameters)):
                    for j in range(i,len(opt.parameters)):
                        counts = np.zeros(len(opt.parameters))
                        counts[i]+=1
                        counts[j]+=1
                        counts[n]+=1
                        print(counts)
                        def func(x,counts=counts):
                            return np.prod(x**counts)
                        sparse_points=df_sparse[df_sparse["k"]==k+1][opt.parameters].to_numpy()
                        weights=df_sparse[df_sparse["k"]==k+1]["weight"].to_numpy()
                        sparse_points=np.apply_along_axis(func,axis=1,arr=sparse_points)
                        sparse_points*=weights
                        mcmc_points=np.apply_along_axis(func,axis=1,arr=df_mcmc[opt.parameters].to_numpy())
                        res_array[k]+=np.abs(sparse_points.sum()-mcmc_points.mean())
            #print(m.factorial(len(opt.parameters)+k-1)/(m.factorial(len(opt.parameters)-1)*m.factorial(k)) + m.factorial(len(opt.parameters)+k-2)/(m.factorial(len(opt.parameters)-1)*m.factorial(k-1))+m.factorial(len(opt.parameters)+k-3)/(m.factorial(len(opt.parameters)-1)*m.factorial(k-2)))
            res_array[k]/=(m.factorial(len(opt.parameters)+k-1)/(m.factorial(len(opt.parameters)-1)*m.factorial(k)) + m.factorial(len(opt.parameters)+k-2)/(m.factorial(len(opt.parameters)-1)*m.factorial(k-1))+m.factorial(len(opt.parameters)+k-3)/(m.factorial(len(opt.parameters)-1)*m.factorial(k-2)))
   
            
    with open(os.path.join("Experiments",opt.name,name+"_"+str(res_array.mean())+".json"),"w") as f:
        json.dump(res_array.tolist(),f)


    return res_array
