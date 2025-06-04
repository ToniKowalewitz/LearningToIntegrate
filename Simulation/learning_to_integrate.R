# script to simulate random fields
# and evaluate them on the nodes 
# and quadrature points of an FE-mesh

# written by hanno gottschalk gottschalk@math.tu-berlin.de
# june 2025
# for the joint project 'learning to integrate'
# with patrick krüger, toni kowalewitz and oliver ernst.

# required r packages
# rgl (for visualization)
# akima for interpolation
# SparseGrid for sparse Grid generation

# load required packages
require(rgl)
require(akima)
require(SparseGrid)
require(expm)
# define settings for mathern smoothing function and strength of noise
settings=list(
  # set grid dimensions
  LD= 100,
  # set shape parameter aka 
  # strength of noise aka variance depends linearly on lambda
  lambda=.05,
  # set Matern class params
  m=.1,# inverse correlation length scale in lattice step units
  nu=3 # smoothness of Mathern function (higher is smoother)
)

######################
#                    #
#  helper functions  #
#                    #
######################

# function to generate a mask for the modes of the random field. 
# by this mask, all nodes not centred at the zero mode are clipped out.
# note that centering in a torus geometry involves all four 
# angles of the array that stores the mode data
#
# parameters:
# LD:     edge size of the quadratic lattice
# approx: size of the non-clipped region is (2*(approx-1)+1)^2
#         i.e. approx = 1 one mode, approx=2 9 modes, approx 3 = 25 modes...
#         note that every complexmode contains a real and a imaginary part, i.e.
#         two real modes
#         approx = 0 returns a matrix with all entries = 1, i.e. no clipping
#
# value: a matrix of entries 1 and 0, where the nearest neighbor modes to 0 are
#        not clipped (entry 1) and al the other modes are clipped (entry 0)

mask.rf=function(LD,approx){
  mask=matrix(1,ncol=LD,nrow=LD)
  if(approx){
    mask[0:(LD-1)>=approx&0:(LD-1)<=LD-approx,]=0
    mask[,0:(LD-1)>=approx&0:(LD-1)<=LD-approx]=0
  }
  return(mask)
}

# function to generate grid values of a lattice
# green function of Matern class in Fourier space
#
# arguments
# setting: list with entries m and nu (see above)
# approx: level of masking modes in fourier space, 
#         defaults to 0 (no masking)

matern.rf=function(settings,approx=0){
  #extract settings
  LD=settings$LD;lambda=settings$lambda;m=settings$m;nu=settings$nu
  G_hat=function(k) m^(2*nu)/(2*(2-cos(k[1])-cos(k[2]))+m^2)^nu
  # grid in fourier space
  k_Grid=expand.grid(k1=0:(LD-1),k2=0:(LD-1))*2*pi/LD
  G_hat_Mat=apply(k_Grid,1,G_hat)
  G_hat_Mat=matrix(G_hat_Mat,ncol=LD)
  if(approx){
    mask=mask.rf(LD,approx)
    G_hat_Mat=G_hat_Mat*mask
  }
  return(G_hat_Mat)
}

###########################
#                         #
#  functions to simulate  #
#  random fields          #
#                         #
###########################

# function to generate a random field and its modes
#
# arguments
# law:   a code word that contains the law of the noise
#        "norm","pois","gamma","bigamma" distributions are supported
# settg: list of settings for the matern function and noise strength
# approx:integer that specifes the level of approximation - see above 
# seed: an optional setting for the random seed. defaults to NULL (no seed)
#
# value: a list containing 
#       $field: an array of real numbers containing the field values
#       $modes: the associated array of modes obtained by fft.

sim.rf=function(law="norm",settg=settings,approx=0,seed=NULL){
  if(!law %in% c("norm","pois","gamma","bigamma")) return("law of noise unknown - no output")
  #extract settings
  LD=settg$LD
  # simulate lattice noise
  if(!is.null(seed)) set.seed(seed)
  if(law=="norm") noise=rnorm(LD*LD,sd=sqrt(settg$lambda/LD)) # Gaussian WN
  if(law=="pois") noise=rpois(LD*LD,lambda = settg$lambda/LD) # Poisson
  if(law=="gamma") noise=rgamma(LD*LD,shape=settg$lambda/LD,scale = 1) # Gamma
  if(law=="bigamma") noise=rgamma(LD*LD,shape=settg$lambda/LD,scale = 1/sqrt(2)) -rgamma(LD*LD,shape=settg$lambda/LD,scale = 1/sqrt(2))# Variance Gamma
  noise=matrix(noise,ncol=LD,nrow=LD) # fill into spatial pattern
  # generate 
  G_hat_mat=matern.rf(settg,app=approx)
  # simulate random field
  noise_hat=fft(noise,inverse = FALSE)
  phi_hat=noise_hat*G_hat_mat
  phi=fft(phi_hat,inverse=TRUE)/LD^2
  phi=Re(phi)
  return(list(field=phi,modes=phi_hat,settings=settg))
}

# an utility function to approximate a fixed realization of a
# random field via its modes. the function in addition quantifies
# a number of distance measures 
# 
# arguments
# phi: a list given as the output of the sim.rf function
# approx: selects the level of approximation, see above
#
# value: a list with the input phi copied and in addition
#        $approx: the level of approimation used
#        $field_approx: the approximated field
#        $metrics: a list containing the L1, L2 and Linf 
#                  distances between $field and $field_approx

approx.rf=function(phi,approx=0){
  LD=ncol(phi$field)
  mask=mask.rf(LD=LD,approx=approx)
  phi_approx=fft(phi$modes*mask,inverse=TRUE)/LD^2
  phi_approx=Re(phi_approx)
  metrics=list(L1=sum(abs(phi$field-phi_approx)),
              L2=sqrt(sum(abs(phi$field-phi_approx)^2)),
              Linf=max(abs(phi$field-phi_approx))
  )
  phi$approx=approx
  phi$field_approx=phi_approx
  phi$metrics=metrics
  return(phi)
}
#########################
#                       #
#  visualizition tools  #
#                       #
#########################

# plot function for a random field
#
# arguments 
# phi: a list as the output of sim.rf 
# cont: logical, should also a 2D contour plot be drawn?
#
# value: NULL - but opens a rgl device and displays the field

plot.rf=function(phi,cont=FALSE){
  LD=phi$settings$LD
  try(close3d(),silent = TRUE)
  open3d()
  surface3d(x=seq(-1,1,length.out=LD),y=seq(-1,1,length.out=LD),z=phi$field,col="lightblue")
  bg3d(color = "white") #"#263F58"
  aspect3d(1,1,1)
  decorate3d()
  if(cont) contour(z=phi$field)
return(NULL)
}

# plot functions that compares a realization
# of a random field with its modal approximation
#
# arguments
# phi: a list as the output of approx.rf
# value: NULL - but opens a rgl device and displays the field
#        and its approximation

plot.approx.rf=function(phi){
  LD=phi$settings$LD
  plot.rf(phi,cont=FALSE)
  surface3d(x=seq(-1,1,length.out=LD),y=seq(-1,1,length.out=LD),z=phi$field_approx,col="lightgrey",alpha=0.3)
  return(NULL)
}

###################################
#                                 #
#  utility functions to simulate  #
#   distributions of modes        #
#                                 #
###################################

# function to extract complex modes k=(0,0),...,k=(0,approx-1)
# and combine them to a real multi-mode distribution
# const, cos-modes (real parts), sin-modes (imaginary parts)
#
# arguments
# n:      number of samples to be drawn
# law:    distribution of noise - see above
# approx: integer >0: controls modes to be extracted 2*(approx-1)+1
#
# value: table of mode samples with n rows and 2*(approx-1)+1 columns

mode.dist.rf=function(n=300,law="norm",approx=0){
  if(approx==0) return("approx has to have a value >0")
  modes=matrix(0,ncol=2*(approx-1)+1,nrow=n)
  for(i in 1:n){
    phi=sim.rf(law=law, approx=approx)
    modes_i=Re(phi$modes[1,1])
    if(approx>1) modes_i=c(modes_i,Re(phi$modes[1,2:approx]),Im(phi$modes[1,2:approx]))
    modes[i,]=modes_i
  }
  return(modes)
}

# this function provides statistical comparison
# of single modes with the normal distribution
# using the shapiro-wilkinson test and the ks test
#
# argument:
# dist: an array containing samples from the real modes of the 
#       random field as produced by mode.dist.rf
# 
# value
# A list containing the ks and sw test results for each marginal

stat.eval.rf=function(dist=mode.dist.rf()){
  tests=list()
  for(i in 1:ncol(dist)){
    tests[[i]]=list(KS=ks.test(dist[,i],'pnorm'),WS=wilcox.test(dist[,i]))
  }
  pvalues=matrix(0,ncol=ncol(dist),nrow=2)
  ncos=(ncol(dist)-1)/2
  if(ncos==0) colnames(pvalues)="constant_mode"
  else colnames(pvalues)=c("constant_mode",paste0("cos_mode",1:ncos),paste0("sin_mode",1:ncos))
  for(i in 1:ncol(dist)){
    pvalues[1,i]=tests[[i]]$KS$p.value
    pvalues[2,i]=tests[[i]]$WS$p.value
  }
  rownames(pvalues)=c("KS-test","WS-test")
  tests$pvalues=pvalues
  return(tests)
}

###############################################
#                                             #
#  functions to read and handle the FE mesh   #
#     and interpolate a random field          #
#  to the quadrature points of the mesh       #
#                                             #
###############################################

# function to read fe mesh data from files (2D triangular mesh)
# and computes the quadrature points as COM of the elements
#
# arguments
# node_file:     txt-file that contains xy coordinates of a node
# element_file:  txt-file that contains the node numbers of the
#                element in each row
#
# value
# A list with three elements
# $elements:   a matrix containing the node numbers of the elements
# $nodes:      a matrix containg the xy-coordinates of the nodes
# $quadpoints: a mateix containing the xy coordinates of the quadrature points 

read.fe=function(file_nodes="mesh-48-nodes.txt",file_elements="mesh-48-elements.txt"){
  nodes=read.csv(file_nodes,header=FALSE)
  colnames(nodes)=c("x","y")
  elements=read.csv(file_elements,header=FALSE)
  colnames(elements)=paste0("N",1:3)
  quadpoints=matrix(0,nrow=nrow(elements),ncol=2)
  colnames(quadpoints)=c("x","y")
  for(i in 1:nrow(elements)) quadpoints[i,]=colMeans(nodes[as.integer(elements[i,]),])
  mesh=list(nodes=as.data.frame(nodes),elements=as.data.frame(elements),quadpoints=as.data.frame(quadpoints))
  return(mesh)
}

# function to visualize the fe mesh
#
#
# argument
#
# mesh: a list containing the fe mesh as described in the 
#       output of read.fe
#
# value
# NULL - but displays a plot of fe mash and quadrature points

plot.fe=function(mesh=read.fe()){
  plot(mesh$nodes,col="lightblue",pch=21,cex=1,bg="lightblue")
  for(i in 1:nrow(mesh$elements)){
    polyg=mesh$nodes[as.integer(mesh$elements[i,]),]
    polygon(x=polyg$x,y=polyg$y,col="lightblue")
  }
  points(mesh$quadpoints,col="lightgrey",bg="lightgrey",pch=21)
  return(NULL)
}

# function that interpolates the random field on the points given in newpoints
# (usually they would be the quadrature points of the mesh)
#
# arguments
#
# phi: a list as in the output of sim.rf - $field gives the field values
# newpoints: a data frame with $x and $y arguments that containes the xy 
#            coordinates to interpolate to
#
# value: a data frame that contains newpoints with an additional $z entry
#        that contains the interpolated field values

interpolate.rf=function(phi,newpoints,app=FALSE){
  LD=nrow(phi$field)
  x=seq(-1,1,length.out=LD)
  y=x
  grid=expand.grid(x=x,y=y)
  if(app) newpoints$z=interpp(x=grid$x,y=grid$y,z=as.numeric(phi$field_approx),xo=newpoints$x,yo=newpoints$y)$z 
  else newpoints$z=interpp(x=grid$x,y=grid$y,z=as.numeric(phi$field),xo=newpoints$x,yo=newpoints$y)$z
  return(newpoints)
}

# function for the visualization, how a rf is interpolated to (quadrature) points
#
# arguments
# phi: a list as provided by the output of sim.rf
#
# value
# NULL - but opens a rgl device with the visualization


plot.interp.rf=function(phi,mesh,up_lines=TRUE){
  range_field=range(phi$field)
  zo=range_field[1]-.2*diff(range_field)
  plot.rf(phi)
  interp_field=interpolate.rf(phi,mesh$quadpoints)
  spheres3d(x=interp_field$x,y=interp_field$y,z=interp_field$z,radius = .02,col="lightgrey")
  spheres3d(x=mesh$quadpoints$x,y=mesh$quadpoints$y,z=zo,col="lightgrey",radius=.02)
  spheres3d(x=mesh$nodes$x,y=mesh$nodes$y,z=zo,col="lightblue",radius=.02)
  for(i in 1:nrow(mesh$elements)){
    polyg=mesh$nodes[as.integer(mesh$elements[i,]),]
    polyg=rbind(polyg,polyg[1,])
    lines3d(x=polyg$x,y=polyg$y,z=zo,col="lightblue",alpha=1)
    if(up_lines) lines3d(x=rep(mesh$quadpoints$x[i],2),y=rep(mesh$quadpoints$y[i],2),z=c(zo,interp_field$z[i]),col="lightgrey",alpha=.3)
  }
  return(NULL)
}

#################################
#                               #
#  functions for extraction of  #
#   independent modes from the  #
#  fft of a real valued funct.  #
#                               #
#################################

# utility function to check symmetry of a complex, quadratic matrix
# that has been obtained by 2Dfft from a real quadratic matrix

# arguments
# mat: a quadratic matrix with complex numbers
#
# value
# a quadratic matrix containing the torus-reflected and complex conjugated of mat

reflect.mat=function(mat){
  L=ncol(mat)
  index=1:L
  ref_index=c(1,rev(tail(index,-1)))
  ref_mat=mat[ref_index,ref_index]
  return(Conj(ref_mat))
}

# utility function to extract L^2 independent real modes from a quadratic matrix 
# that contains compex numbert obtained a the fft of a real valued quadratic matrix
# this utility serves, to eliminate redundant modes determined by the symmetry
#
# argument
# mat: a quadratic LxL matrix with complex values
# approx: Level of approximation as explained in mask.rf 
#
# value
# a vector of L^2 real numbers representing the independent modes.
# Usually, L=(2(approx-1)+1), but if approx =0, L=ncol(mat).

extract.modes=function(mat,approx=0){
  LD=ncol(mat)
  mask=mask.rf(LD = LD,approx = approx)
  if(approx>0){
    ind=rowSums(mask)>0
    mat=mat[ind,ind]
    if( !is.matrix(mat) ) mat=matrix(mat,ncol=1) 
  }
  L=ncol(mat)
  if(L==1) return(as.numeric(Re(mat)))
  index_mat=matrix(1:((L-1)^2),ncol=L-1)
  nReal=ceiling((L-1)/2)
  nIm=floor((L-1)/2)
  nReMat=ceiling((L-1)^2/2)
  nImMat=floor((L-1)^2/2)
  modes_re_null_null=Re(mat[1,1])
  modes_re_null_x=Re(mat[1,2:(1+nReal)])
  modes_re_x_null=Re(mat[2:(1+nReal),1])
  modes_re_x_x=as.numeric(Re(mat[2:L,2:L]))[index_mat<=nReMat]
  modes_im_null_x=Im(mat[1,2:(1+nIm)])
  modes_im_x_null=Im(mat[2:(1+nIm),1])
  modes_im_x_x=as.numeric(Im(mat[2:L,2:L]))[index_mat<=nImMat]
  modes=c(modes_re_null_null,
          modes_re_null_x,
          modes_re_x_null,
          modes_re_x_x,
          modes_im_null_x,
          modes_im_x_null,
          modes_im_x_x)
  return(modes)
}

# a utility function to reverse the previous operation,
# i.e. to reconstruct the complex LxL matrix with the symmetry of
# fft of a real valued matrix from L^2 real modes
#
# arguments
# modes: a real valued vector of legth L^2
# LD: edge size of the array, into which the modes are to be inserted
#
# value
# a LxL complex matrix with the prescribed symmetry


insert.modes=function(modes,LD=settings$LD){
  L=sqrt(length(modes))
  if(L<LD) approx=(L-1)/2+1
  else approx=0
  if(L>1){
    index_mat=matrix(1:((L-1)^2),ncol=L-1)
    nReal=ceiling((L-1)/2)
    nIm=floor((L-1)/2)
    nReMat=ceiling((L-1)^2/2)
    nImMat=floor((L-1)^2/2)
    modes_re_null_null=head(modes,1);modes=tail(modes,-1)
    modes_re_null_x=head(modes,nReal);modes=tail(modes,-nReal)
    modes_re_x_null=head(modes,nReal);modes=tail(modes,-nReal)
    modes_re_x_x=head(modes,nReMat);modes=tail(modes,-nReMat)
    modes_im_null_x=head(modes,nIm);modes=tail(modes,-nIm)
    modes_im_x_null=head(modes,nIm);modes=tail(modes,-nIm)
    modes_im_x_x=modes
    mat_re=matrix(0,ncol=L,nrow=L)
    mat_im=mat_re 
    mat_re[1,1]=modes_re_null_null
    if((L-1)%%2==0){
      mat_re[1,2:(1+2*nReal)]=c(modes_re_null_x,rev(modes_re_null_x))
      mat_re[2:(1+2*nReal),1]=c(modes_re_x_null,rev(modes_re_x_null))
      mat_re[2:L,2:L]=matrix(c(modes_re_x_x,rev(modes_re_x_x)),ncol=L-1)
      mat_im[1,2:(1+2*nReal)]=c(modes_im_null_x,-rev(modes_im_null_x))
      mat_im[2:(1+2*nReal),1]=c(modes_im_x_null,-rev(modes_im_x_null))
      mat_im[2:L,2:L]=matrix(c(modes_im_x_x,-rev(modes_im_x_x)),ncol=L-1)
      }
    else{ 
      mat_re[1,2:(2*nReal)]=c(modes_re_null_x,tail(rev(modes_re_null_x),-1))
      mat_re[2:(2*nReal),1]=c(modes_re_x_null,tail(rev(modes_re_x_null),-1))
      mat_re[2:L,2:L]=matrix(c(modes_re_x_x,tail(rev(modes_re_x_x),-1)),ncol=L-1)
      mat_im[1,2:(2*nReal)]=c(modes_im_null_x,0,-rev(modes_im_null_x))
      mat_im[2:(2*nReal),1]=c(modes_im_x_null,0,-rev(modes_im_x_null))
      mat_im[2:L,2:L]=matrix(c(modes_im_x_x,0,-rev(modes_im_x_x)),ncol=L-1)
    }
  } 
  else{
    mat_re=matrix(modes,ncol=1)
    mat_im=matrix(0,ncol=1,nrow=1)
  }
  mat=matrix(complex(real=mat_re,imaginary=mat_im),ncol=L,nrow=L)
  mat_full=matrix(0,ncol=LD,nrow=LD)
  if(approx>0) mask=mask.rf(LD=LD,approx=(L-1)/2+1)
  else mask=mask.rf(LD=LD,approx=0)
  ind=rowSums(mask)>0
  mat_full[ind,ind]=mat
  phi=sim.rf()
  phi$modes=mat_full
  phi$field=fft(phi$modes*mask,inverse=TRUE)/LD^2
  return(phi)  
}


###############################
#                             #
# scripts for data generation #
#                             #
###############################

# function to generate a monte carlo distribution of modes and interpolated paths for various degrees of approximation.
#
# arguments:
# nsample: number of MC samples
# approx: level of approximation - as in sim.rf
# law: code of the the levy-law, one from from "norm", "pois", "gamma" "bigamma"
# setg: The settings of the mathern covariance and the strength of noise - defaults to settings as defined in the header
# path: where to write the output files
#
generate.distribution=function(nsample=1e4,approx=0,law="norm",setg=settings,path="./data/",write_file=TRUE,write_modes=FALSE,interp=TRUE,mesh_resolution="M"){
  if(mesh_resolution=="C") {
    elements="mesh-48-elements.txt"
    nodes="mesh-48-nodes.txt"
  }
  if(mesh_resolution=="M") {
    elements="mesh-798-elements.txt"
    nodes="mesh-798-nodes.txt"
  }
  if(mesh_resolution=="F") {
    elements="mesh-8510-elements.txt"
    nodes="mesh-8510-nodes.txt"
  }
  options(warn=-1)
  phi=sim.rf(law = law,settg = setg,approx = approx)
  modes=extract.modes(phi$modes,approx=approx)
  mode_mat=matrix(0,ncol=length(modes),nrow=nsample)
  colnames(mode_mat)=paste0("mode",1:ncol(mode_mat))
  if(interp){
    mesh=read.fe(file_nodes = nodes,file_elements = elements)
    interp_mat=matrix(0,ncol=nrow(mesh$quadpoints),nrow=nsample)
    colnames(interp_mat)=paste0("quad_pt",1:ncol(interp_mat))
  }
  pb = txtProgressBar(min = 0, max = nsample, initial = 0) 
  for ( i in 1:nsample){
    setTxtProgressBar(pb,i)
    phi=sim.rf(law = law,settg = setg,approx = approx)
    modes=extract.modes(phi$modes,approx=approx)
    mode_mat[i,]=modes
    if(interp){
      interpolation=interpolate.rf(phi,mesh$quadpoints)
      interp_mat[i,]=interpolation$z
    }
  }
  if(write_file){
    if(write_modes) write.table(mode_mat,file=paste0(path,"modes_A_",approx,"_L_",law,"_N_",nsample,".txt"),col.names = TRUE,row.names = FALSE)
    if(interp) write.table(interp_mat,file=paste0(path,"mesh_",mesh_resolution,"_interp_A",approx,"_L_",law,"_N_",nsample,".txt"),col.names = TRUE,row.names = FALSE)
  }
  options(warn=0)
  close(pb)
  
  if(interp) return(list(mode_dist=mode_mat, rf_on_quad_pt=interp_mat))
  else return(list(mode_dist=mode_mat))
}

# function to generate interpolations of ONE random field on different FE meshes of resolution F=fine
# M= medium and C= coarse
#
# attention: the three mesh data sets are hard coded. 
# 
# arguments:
#
# nSample= The number of fields to be produced
# approx: The level of mode approximations
# law: one of "norm", "pois", "gamma" or "bigamma" selecting the Levy distribution
# setg: the setings of noise strength and Matern kernel, defaults to the settings list defined in the header
# path: path to write the output files to
# write_file: Locical: should the output be written to files?
# write_modes: Should mode files be written?
#
# value:
# A list with entries F,M,C containing data frames with the values interpolated to the quadrature
# points of th elements of the fine (F), medium (M) or coarse (C) grids.


generate.comparison=function(nsample=1e2,approx=c(0,2,3,4),law="norm",setg=settings,path="./FE_comparison//",write_file=TRUE,write_modes=FALSE){
  options(warn=-1)
  phi=sim.rf(law = law,settg = setg,approx = 0)
  mesh_C=read.fe(file_nodes = "mesh-48-nodes.txt",file_elements = "mesh-48-elements.txt")
  mesh_M=read.fe(file_nodes = "mesh-798-nodes.txt",file_elements = "mesh-798-elements.txt")
  mesh_F=read.fe(file_nodes = "mesh-8510-nodes.txt",file_elements = "mesh-8510-elements.txt")
  
  interp_mat_C=matrix(0,ncol=nrow(mesh_C$quadpoints),nrow=length(approx)*nsample)
  colnames(interp_mat_C)=paste0("quad_pt",1:ncol(interp_mat_C))
  interp_mat_C=as.data.frame(interp_mat_C)
  interp_mat_C$approx=0
  
  interp_mat_M=matrix(0,ncol=nrow(mesh_M$quadpoints),nrow=length(approx)*nsample)
  colnames(interp_mat_M)=paste0("quad_pt",1:ncol(interp_mat_M))
  interp_mat_M=as.data.frame(interp_mat_M)
  interp_mat_M$approx=0
  
  
  interp_mat_F=matrix(0,ncol=nrow(mesh_F$quadpoints),nrow=length(approx)*nsample)
  colnames(interp_mat_F)=paste0("quad_pt",1:ncol(interp_mat_F))
  interp_mat_F=as.data.frame(interp_mat_F)
  interp_mat_F$approx=0
  
  
  pb = txtProgressBar(min = 0, max = nsample, initial = 0) 
  for ( i in 1:nsample){
    setTxtProgressBar(pb,i)
    phi=sim.rf(law = law,settg = setg,approx = 0)
    for(j in 1:length(approx)){
    phi_a=approx.rf(phi,approx=approx[j])  
    
    interpolation_C=interpolate.rf(phi_a,mesh_C$quadpoints,app=TRUE)
    interp_mat_C[(i-1)*length(approx)+j,1:(ncol(interp_mat_C)-1)]=interpolation_C$z
    interp_mat_C$approx[(i-1)*length(approx)+j]=approx[j]
    
    
    interpolation_M=interpolate.rf(phi_a,mesh_M$quadpoints,app=TRUE)
    interp_mat_M[(i-1)*length(approx)+j,1:(ncol(interp_mat_M)-1)]=interpolation_M$z
    interp_mat_M$approx[(i-1)*length(approx)+j]=approx[j]
  
    interpolation_F=interpolate.rf(phi_a,mesh_F$quadpoints,app=TRUE)
    interp_mat_F[(i-1)*length(approx)+j,1:(ncol(interp_mat_F)-1)]=interpolation_F$z
    interp_mat_F$approx[(i-1)*length(approx)+j]=approx[j]
    }
  }
  if(write_file){
    write.table(interp_mat_C,file=paste0(path,"comparison_mesh_C_L_",law,"_N_",nsample,".txt"),col.names = TRUE,row.names = FALSE)
    write.table(interp_mat_M,file=paste0(path,"comparison_mesh_M_L_",law,"_N_",nsample,".txt"),col.names = TRUE,row.names = FALSE)
    write.table(interp_mat_F,file=paste0(path,"comparison_mesh_F_L_",law,"_N_",nsample,".txt"),col.names = TRUE,row.names = FALSE)
  }
  options(warn=0)
  close(pb)
  return(list(C=interp_mat_C,M=interp_mat_M,F=interp_mat_F))
}

# A utility function to create Smolyak sparse grids for multivariate standard normal distribution 
# and write files containing quadrature points and weights to disk.
# The degree of exactness is determined by the maximum number of allowed quadrature points
#
# arguments
# dimension: The dimension of the multivariate normal distribution
# path: where to write the files
# nmax: the maximum number of quadrature points allowed.
# write_file: should  results files be written?
#
# value:
# a data frame containing quadrature points, weights and the degree of exactness k 
# (total degree of exactly integrated polynoms is (2k-1)) in a row


generate.grids=function(dimension=9,path="./grids/", nmax=20000,write_file=TRUE){
  N=0
  k=1
  all_grids=data.frame(t(rep(0,dimension+2)))[FALSE,]
  names(all_grids)=c(paste0("Z",1:dimension),"weight", "k")
    
  while(N<=nmax){
  grid=createSparseGrid(type="GQN",dimension=dimension,k=k)
  N=nrow(grid$nodes)
  tmp_grid=as.data.frame(grid$nodes)
  names(tmp_grid)=paste0("Z",1:dimension)
  tmp_grid$weight=grid$weights
  tmp_grid$k=k
  if(N<=nmax) all_grids=rbind(all_grids,tmp_grid)
  k=k+1
  }
  if(write_file) write.table(file=paste0(path,"sparse_grids_d_",dimension,".txt"),x = all_grids,row.names = FALSE)
  return(all_grids)
}

# function that reads the INN output files of transformed quadratures, 
# creates a random field from the INN-transformed sparse grid quadrature points
# and interpolates it to a finite element mesh

interpolate.sparse.grid.fields=function(read_path="./INN_outputs/",write_path="./sparse_grid_fields/",Noise=FALSE){
  files=dir(read_path)
  for(i in 1:length(files)){
    cat("processing file ",i," of",length(files)," filename ", files[i],"\n")
    approx=0
    if(length(grep("A1",files[i]))>0) approx=1
    if(length(grep("A2",files[i]))>0) approx=2
    if(length(grep("A3",files[i]))>0) approx=3
    if(length(grep("A4",files[i]))>0) approx=4
    if(length(grep("A5",files[i]))>0) approx=5
    if(length(grep("A6",files[i]))>0) approx=6
    
    mesh_C=read.fe(file_nodes = "mesh-48-nodes.txt",file_elements = "mesh-48-elements.txt")
    mesh_M=read.fe(file_nodes = "mesh-798-nodes.txt",file_elements = "mesh-798-elements.txt")
    mesh_F=read.fe(file_nodes = "mesh-8510-nodes.txt",file_elements = "mesh-8510-elements.txt")
    
    sparse_grid_transformed=read.csv(file = paste0(read_path,files[i]))
    sparse_grid_transformed=sparse_grid_transformed[,2:ncol(sparse_grid_transformed)]
    
    nsample=nrow(sparse_grid_transformed)
    
    interp_mat_C=matrix(0,ncol=nrow(mesh_C$quadpoints),nrow=length(approx)*nsample)
    colnames(interp_mat_C)=paste0("quad_pt",1:ncol(interp_mat_C))
    interp_mat_C=as.data.frame(interp_mat_C)
  
    interp_mat_M=matrix(0,ncol=nrow(mesh_M$quadpoints),nrow=length(approx)*nsample)
    colnames(interp_mat_M)=paste0("quad_pt",1:ncol(interp_mat_M))
    interp_mat_M=as.data.frame(interp_mat_M)
  
    interp_mat_F=matrix(0,ncol=nrow(mesh_F$quadpoints),nrow=length(approx)*nsample)
    colnames(interp_mat_F)=paste0("quad_pt",1:ncol(interp_mat_F))
    interp_mat_F=as.data.frame(interp_mat_F)

    pb = txtProgressBar(min = 0, max = nrow(sparse_grid_transformed), initial = 0) 
      
    for(j in 1:nrow(sparse_grid_transformed)){
        setTxtProgressBar(pb,j)
        if(Noise) phi=insert.modes(modes = as.numeric(sparse_grid_transformed[j,]))
        else phi=insert.modes(modes = as.numeric(sparse_grid_transformed[j,1:(ncol(sparse_grid_transformed)-2)]))
        interp_mat_C[j,]=interpolate.rf(phi,newpoints = mesh_C$quadpoints)$z
        interp_mat_M[j,]=interpolate.rf(phi,newpoints = mesh_M$quadpoints)$z
        interp_mat_F[j,]=interpolate.rf(phi,newpoints = mesh_F$quadpoints)$z
    }
    if(!Noise){
    interp_mat_C=cbind(interp_mat_C,sparse_grid_transformed[,(ncol(sparse_grid_transformed)-1):ncol(sparse_grid_transformed)])
    interp_mat_M=cbind(interp_mat_M,sparse_grid_transformed[,(ncol(sparse_grid_transformed)-1):ncol(sparse_grid_transformed)])
    interp_mat_F=cbind(interp_mat_F,sparse_grid_transformed[,(ncol(sparse_grid_transformed)-1):ncol(sparse_grid_transformed)])
    }
    
    write.csv(interp_mat_C,file=paste0(write_path,"sparse_grid_fields_mesh_C_",strsplit(files[i],c("modes"))[[1]][2]))
    write.csv(interp_mat_M,file=paste0(write_path,"sparse_grid_fields_mesh_M_",strsplit(files[i],c("modes"))[[1]][2]))
    write.csv(interp_mat_F,file=paste0(write_path,"sparse_grid_fields_mesh_F_",strsplit(files[i],c("modes"))[[1]][2]))
    close(pb)
  }
  return(NULL)
  
}


sparse.grid.normal.test=function(modes="./data/modes_A2_L_norm_N_10000.txt",grids="./grids/sparse_grids_d_9.txt",out_file="./sparse_grid_test_linear/LINEAR_OUTPUT_INN_LTI2_modes_A2_L_norm_N.csv"){
  normal_modes=read.table(file=modes,header = TRUE)
  grid_data=read.table(file=grids,header = TRUE)
  covariance_mat=cov(normal_modes)
  means=colMeans(normal_modes)
  sqrt_covariance_mat=sqrtm(covariance_mat)
  grid_data[1:(ncol(grid_data)-2)]=t(sqrt_covariance_mat%*%t(as.matrix(grid_data[1:9])))
  write.csv(grid_data,file=out_file)
  return(grid_data)
}


#################
#               #
# End of Code   #
#               #
#################


##########################
#                        #  
# Scripts for evaluation #
# uncomment to run       #
#                        #
##########################

# interpolate.sparse.grid.fields()

# cat("Processing block - approx 0 \n")
# generate.distribution(law="norm",approx = 0,nsample=1e4,write_modes=FALSE,mesh_resolution = "M")
# generate.distribution(law="pois",approx = 0,nsample=1e4,write_modes=FALSE,mesh_resolution = "M")
# generate.distribution(law="gamma",approx = 0,nsample=1e4,write_modes=FALSE,mesh_resolution = "M")
# generate.distribution(law="bigamma",approx = 0,nsample=1e4,write_modes=FALSE,mesh_resolution = "M")


# cat("Processing block - approx 2")
# generate.distribution(law="norm",approx = 2,nsample=1e5,write_modes=TRUE,interp=FALSE,mesh_resolution = "C")
#generate.distribution(law="pois",approx = 2,nsample=1e5,write_modes=TRUE,interp=FALSE,mesh_resolution = "C",path="./data/large/")
#generate.distribution(law="gamma",approx = 2,nsample=1e5,write_modes=TRUE,interp=FALSE,mesh_resolution = "C",path="./data/large/")
#generate.distribution(law="bigamma",approx = 2,nsample=1e5,write_modes=TRUE,interp=FALSE,mesh_resolution = "C",path="./data/large/")

# cat("Processing block - approx 3")
# generate.distribution(law="norm",approx = 3,nsample=1e4,write_modes=FALSE,mesh_resolution = "M")
#
#generate.distribution(law="norm",approx = 3,nsample=1e5,write_modes=TRUE,interp=FALSE,mesh_resolution = "C")
#generate.distribution(law="pois",approx = 3,nsample=1e5,write_modes=TRUE,interp=FALSE,mesh_resolution = "M",path="./data/large/")
#generate.distribution(law="gamma",approx = 3,nsample=1e5,write_modes=TRUE,interp=FALSE,mesh_resolution = "M",path="./data/large/")
#generate.distribution(law="bigamma",approx = 3,nsample=1e5,write_modes=TRUE,interp=FALSE,mesh_resolution = "M",path="./data/large/")

# Testing SG without INN
# sparse.grid.normal.test()
#interpolate.sparse.grid.fields(read_path = "./INN_outputs_norm/",write_path = "./sparse_grid_fields_norm/")