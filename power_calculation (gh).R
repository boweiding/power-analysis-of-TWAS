## parameters
## x: the matrix of snps, rows are snps, cols are patients.
## x_id: a column with all IDs of snps in x
## kinshipx: the square kinship matrix of snps. 
## zmap: the data set which records the genes and their surrounded snps. rows are genes, cols are IDs of snps(start from the third col).
## i: the number of genes we are going to fit in the regression models.
## ht: trait heritability. (squared), x to y directly
## h1: expression heritability, x to z directly
## h2: PVX, z to y
## num_snp: the number of snps
## num_gene: the number of genes
## nind: the number of patients
## ql_x: when x is the independent varible, the critical value under the central t distribution.
## ql_z: when z is the independent varible, the critical value under the central t distribution.
## min, max: the range of contributed SNPs
## q: number of simulations
## pa: saving path
args=commandArgs(T)##parameter passing
i<-args[1]
i<-as.numeric(i)
h1<-args[2]
h1<-as.numeric(h1)
h2<-args[3]
h2<-as.numeric(h2)
min<-args[5]
min<-as.numeric(min)
max<-args[5]
max<-as.numeric(max)
pa<-args[5]

x<-read.table('snp.csv',sep="\t",header=F) ###snp data
x<-as.matrix(x)
nind<-length(x[1,])-1## number of patients
print(nind)
num_snp<-length(x[,1])

zmap<-read.table('zmap.csv',sep="\t")
num_gene<-length(zmap[,1]) ##number of genes, have already cancelled the header.
newx<-x[,-1]
newx<-apply(newx, 2, as.numeric)
kinshipx<-read.table("kinshipx.csv",sep = "\t",header = T)
kx<-eigen(kinshipx)
kxvector<-kx$vectors
kxs<-diag(kx$values)
##k=usu (eigenvector decomposition), then calculate D (decorrelation paramenter) 

##loading packages in need
library(lme4)
library(glmnet)
library(foreach)
library(doParallel)
library(erer)
library(Deriv)
print("all packages have been loaded")
registerDoParallel(4)

### y~u+e ;  z~u+e
sigma_est<-function(n,u,s,y,h2){
  s<-diag(s)
  s<-as.vector(s)
  y<-y-mean(y)##centered y
  x<-t(u)%*%y
  x<-as.vector(x)
  f1<-function(de1) sum(n*x^2/((de1+s)^2*sum(x^2/(de1+s)))-1/(de1+s))
  f2<-Deriv(f1,"de1")
  inde<-1
  de1<-1/(1*h2)
  while(TRUE){
    de1<-1/(inde*h2)
    if(eval(f1(de1))>0)
      break
    else
      inde<-inde+1
  }
  de2<-0
  diff<-1
  while(TRUE){
    de2<-de1-eval(f1(de1))/eval(f2(de1)) ##newton method
    if (eval(f1(de2))<=0){
      de2<-de1
      break
    }
    diff<-de2-de1
    de1<-de2
    if (abs(diff)<=0.000001)
      break
  }
  sg<-sum((x^2)/(de2+s))/n
  se<-sg*de2
  result<-list(sigmag=sg,sigmae=se)
  return(result)
}

x_id<-x[,1]
x<-newx#remove the id from x

simulatey<-function(q,x,x_id,z,i,min,max,ht,h1,h2,kxvector,kxs,nind){##simulate y directly from x 
  n_gene<-length(z[,1])## number of gene 
  set.seed(q)
  xnum<-runif(i,min,max)##select the number of snps around each gene
  xnum<-round(xnum)
  xgroup<-NULL
  zgroup<-NULL
  z_num<-NULL
  I<-diag(1,nind,nind)
  xgroupindex<-NULL
  zindex<-NULL
  z_initial<-0
  vnum<-100*q
  knum<-100*q
  v<-1
  while (v <=i) {##within each gene
    set.seed(vnum+v)
    znum<-runif(1,1,n_gene)
    znum<-round(znum)##selecting the genes for regression
    vnum<-vnum+v
    if (znum==z_initial)
      next ## next "for" loop
    else
      z_initial<-znum
    id_temp<-NULL
    for (u in 2:length(z[znum,])){
      id_temp<-c(id_temp,z[znum,u])
    }
    initialxindex<-0
    selectedX<-NULL
    n_temp<-length(id_temp)
    k<-1
    indexgroup_j<-NULL
    while (k<=xnum[v]){
      set.seed(knum+k)
      xindex<-runif(1,1,n_temp)
      xindex<-round(xindex)
      knum<-knum+k
      if (initialxindex==xindex)
        next ## next "while" loop
      else
        initialxindex<-xindex 
      sx<-x[id_temp[xindex]==x_id,]
      selectedX<-c(selectedX,sx)
      indexgroup_j<-c(indexgroup_j,id_temp[xindex])
      k<-k+1
    }
    xmatrix<-matrix(selectedX,nrow=nind,ncol=xnum[v],byrow=F)
    options(digits=3)
    set.seed(10*q*i+v)
    beta<-runif(xnum[v],0,1)##no intercept
    betax<-xmatrix%*%beta
    samplemean<-mean(betax)
    sigmag2_temp<-t(betax-samplemean)%*%(betax-samplemean)/(nind-1)
    sigmag2_temp<-as.vector(sigmag2_temp)
    sigmae2<-sigmag2_temp*(1-h1)/h1
    sumge<-sigmag2_temp+sigmae2
    set.seed(10*q*i+v)
    eps<-rnorm(nind,0,sqrt(sigmae2))
    z_temp<-(betax+eps)/sqrt(sumge)##rescale according to sigmag2+sigmae2=1
    zgroup<-cbind(zgroup,z_temp)
    xgroup<-cbind(xgroup,xmatrix)
    xgroupindex<-c(xgroupindex,indexgroup_j)
    zindex<-c(zindex,z[znum,1])
    z_num<-c(z_num,znum)
    v<-v+1
  } 
  ##already selected z and x and calculated simulated z 
  ##y simulated directly by x
  options(digits=3)
  set.seed(4000*q)
  beta_xy<-runif(sum(xnum),0,1)
  betaxy<-xgroup%*%beta_xy
  samplemeanxy<-mean(betaxy)
  sigmag2xy_temp<-t(betaxy-samplemeanxy)%*%(betaxy-samplemeanxy)/(nind-1)
  sigmag2xy_temp<-as.vector(sigmag2xy_temp)
  sigmae2xy<-sigmag2xy_temp*(1-ht)/ht
  sumge_xy<-sigmag2xy_temp+sigmae2xy
  set.seed(400*q)
  epsxy<-rnorm(nind,0,sqrt(sigmae2xy))
  y_x<-(betaxy+epsxy)/sqrt(sumge_xy)
  ##y is a (n*1) vector
  
  ##y simulated from z
  set.seed(400*q+1)
  beta_zy<-runif(i,0,1)
  betazy<-zgroup%*%beta_zy
  samplemeanzy<-mean(betazy)
  sigmag2zy_temp<-t(betazy-samplemeanzy)%*%(betazy-samplemeanzy)/(nind-1)
  sigmag2zy_temp<-as.vector(sigmag2zy_temp)
  sigmae2zy<-sigmag2zy_temp*(1-h2)/h2
  sumge_yz<-sigmag2zy_temp+sigmae2zy
  set.seed(400*q+1)
  epszy<-rnorm(nind,0,sqrt(sigmae2zy))
  y_z<-(betazy+epszy)/sqrt(sumge_yz)
  
  ##elastic net:
  z_hat<-NULL
  ela_net<-NULL
  options(digits = 5)
  a<-seq(0.1,0.9,0.05) ##alpha
  for (j in 1:i){
    xfit<-NULL
    snpid<-c("Intercept")
    for (k in 2:length(z[z_num[j],]) ) {
      fx<-x[z[z_num[j],k]==x_id,]
      snpid<-c(snpid,z[z_num[j],k])
      xfit<-cbind(xfit,fx)######to be debug
    }
    searcha<-foreach(q=a, .combine = rbind) %dopar% {
      cv<-cv.glmnet(xfit,zgroup[,j],family= "gaussian",nfold=10, type.measure = "deviance",parallel = TRUE, alpha=q)
      data.frame(cvm=cv$cvm[cv$lambda==cv$lambda.1se],lambda.1se=cv$lambda.1se,alpha=q)
    }
    cv3<-searcha[searcha$cvm==min(searcha$cvm),]
    md3<-glmnet(xfit,zgroup[,j],family = "gaussian",lambda = cv3$lambda.1se, alpha = cv3$alpha)
    newz<-predict(md3,xfit)
    mz<-mean(zgroup[,j])
    tss<-t(zgroup[,j]-mz)%*%(zgroup[,j]-mz)
    tss<-as.vector(tss)
    rss<-t(zgroup[,j]-newz)%*%(zgroup[,j]-newz)
    rss<-as.vector(rss)
    r2<-1-(rss/tss)
    coe<-as.vector(coef(md3))
    ela_temp<-data.frame(coef=coe,id=snpid,Rsquared=r2)
    if(j==1)
      ela_net<-ela_temp
    else
      ela_net<-cbind(ela_net,ela_temp)
    z_hat<-cbind(z_hat,newz)
  }
  ##y simulated by x
  sig_yxx<-sigma_est(nind,kxvector,kxs,y_x,h2)
  sig_yxx_g2<-sig_yxx$sigmag
  sig_yxx_e2<-sig_yxx$sigmae
  D1_te<-diag(sig_yxx_g2*kxs+sig_yxx_e2*I)
  D1_te<-diag(D1_te^(-0.5))
  D_yx<-D1_te%*%t(kxvector)
  ##dimension: nind * nind
  
  ##y simulated by z
  sig_yzx<-sigma_est(nind,kxvector,kxs,y_z,h2)
  sig_yzx_g2<-sig_yzx$sigmag
  sig_yzx_e2<-sig_yzx$sigmae
  D3_te<-diag(sig_yzx_g2*kxs+sig_yzx_e2*I)
  D3_te<-diag(D3_te^(-0.5))
  D_yz<-D3_te%*%t(kxvector)
  
  simresult<-list(z=zgroup,x=xgroup,xnum=xnum,xgroupindex=xgroupindex,y_x=y_x,y_z=y_z,z_hat=z_hat,
                  zindex=zindex,z_num=z_num,D_yx=D_yx,D_yz=D_yz,ela_net_result=ela_net)
  return(simresult)
}

ncp<-function(x,y,d,nind){##here x is a vector not a matrix, y is a vector, d and is a matrix.
  x<-d%*%x
  y<-d%*%y
  ##decorrelation, vector
  numerator1<-c(0)
  sumnumeratorxy<-t(x)%*%y
  sumnumeratorxy<-as.vector(sumnumeratorxy)
  sumnumeratory1<-c(0)
  sumnumeratorx1<-c(0)
  for (j in 1:nind){
    Di<-sum(d[j,])
    numerator1<-numerator1+sumnumeratorxy*(Di^2)
    sumnumeratory1<-sumnumeratory1+y[j]*Di
    sumnumeratorx1<-sumnumeratorx1+x[j]*Di
  }
  numerator<-numerator1-sumnumeratory1*sumnumeratorx1
  sumdenomx1<-t(x)%*%x
  sumdenomx1<-as.vector(sumdenomx1)
  sumdenomx2<-c(0)
  sumDi2<-c(0)
  for (j in 1:nind){
    Di<-sum(d[j,])
    sumDi2<-sumDi2+Di^2
    sumdenomx2<-sumdenomx2+Di*x[j]
  }
  denominator<-sqrt(sumdenomx1*(sumDi2^2)-sumDi2*(sumdenomx2^2))
  if (is.nan(denominator)==TRUE || denominator==0)
    denominator<-NaN
  ncp<-numerator/denominator
  return(ncp)
}

ncp_zx<-function(x,z,nind){##z~x is a univariate simple linear regression, no decorrelation needed
  mu_x<-mean(x)
  mu_z<-mean(z)
  tx<-x-mu_x
  beta1<-(t(tx)%*%tx)^(-1)%*%t(tx)%*%z ## least square
  beta1<-as.vector(beta1)
  sig<-0
  for (j in 1:nind){
    numerator<-(z[j]-mu_z+beta1*(mu_x-x[j]))^2
    sig<-sig+numerator
  }
  sig<-sig/(nind-2)
  ncp<-(t(tx)%*%tx)^(-0.5)%*%t(tx)%*%z/sqrt(sig)
  ncp<-as.vector(ncp)
  return(ncp)
}
ql_x<-qt(p=(1-0.05/num_snp),df=nind-2) 
ql_z<-qt(p=(1-0.05/num_gene),df=nind-2)

gwasmodel<-function(sim_x,y_x,y_z,D_yx,D_yz,xnum,xgroupindex,ql_x,nind,i){
  ##i is the num of gene (z)
  ##model 1 (GWAS): Y~X
  model1a<-NULL##to store all values of the statistical power for each linear mixed regression
  model1b<-NULL
  ##...a means y simulate from x
  ##...b means y simulate from z
  numtotalx<-sum(xnum)##the num of x in total
  for (j in 1:numtotalx){
    temp_ncp<-ncp(sim_x[,j],y_x,D_yx,nind)
    temp_power<-1-pt(q=ql_x,df=nind-2,ncp=temp_ncp)
    model1a<-c(model1a,temp_power)
  }
  m1_pt<-1
  for (j in 1:numtotalx){
    m1_pt<-m1_pt*(1-model1a[j])  ##at least find one snp
  }
  m1_pt<-1-m1_pt
  for (k in 1:numtotalx){
    temp_ncp<-ncp(sim_x[,k],y_z,D_yz,nind)
    temp_power<-1-pt(q=ql_x,df=nind-2,ncp=temp_ncp)
    model1b<-c(model1b,temp_power)
  }
  m2_pt<-1
  for (j in 1:numtotalx){
    m2_pt<-m2_pt*(1-model1b[j])  ##at least find one snp
  }
  m2_pt<-1-m2_pt
  m1a<-data.frame(id_snp=xgroupindex,power_a=model1a)
  m1a_t<-data.frame(n_gene=i,Genetic_Model="Pleiotropy", Estimate_Model="GWAS",power=m1_pt)
  m1b<-data.frame(id_snp=xgroupindex,power_b=model1b)
  m1b_t<-data.frame(n_gene=i,Genetic_Model="Causality", Estimate_Model="GWAS",power=m2_pt)
  gwas<-list(GWAS_Pt=m1a_t,GWAS_Ct=m1b_t,GWAS_ple=m1a,GWAS_cau=m1b)
  return(gwas)
}

egwasmodel<-function(sim_x,sim_z,y_x,y_z,D_yx,D_yz,xnum,i,xgroupindex,zindex,zmap,z_num,ql_z,nind){
  ##model 2 (EGWAS)
  power_x1<-NULL ##at least find one snp
  power_zx<-NULL##to store all values of the statistical power for each linear mixed regression
  power_yxz<-NULL
  power_yzz<-NULL
  ##...a means y simulate from x
  ##...b means y simulate from z
  previousnum<-0
  ##z~x
  for (j in 1:i){
    numx<-xnum[j]
    numx<-numx+previousnum
    ininum<-previousnum+1##first loop inimum=1
    model2_zx<-NULL
    n_snp<-length(zmap[z_num[j],])-1
    ql_zx<-qt(p=(1-0.05/n_snp),df=nind-2)
    power_x1_temp<-1
    for (k in ininum:numx){
      temp_ncp<-ncp_zx(sim_x[,k],sim_z[,j],nind)
      temp_power<-1-pt(q=ql_zx,df=nind-2,ncp=temp_ncp)#######
      power_zx<-c(power_zx,temp_power)
      power_x1_temp<-power_x1_temp*(1-temp_power)
    }
    power_x1_temp<-1-power_x1_temp
    power_x1<-c(power_x1,power_x1_temp)
    previousnum<-numx
    temp_ncp_yxz<-ncp(sim_z[,j],y_x,D_yx,nind)
    temp_power_yxz<-1-pt(q=ql_z,df=nind-2,ncp=temp_ncp_yxz)
    power_yxz<-c(power_yxz,temp_power_yxz)
    ###y~z, y simulated by x
    temp_ncp_yzz<-ncp(sim_z[,j],y_z,D_yz,nind)
    temp_power_yzz<-1-pt(q=ql_z,df=nind-2,ncp=temp_ncp_yzz)
    power_yzz<-c(power_yzz,temp_power_yzz)
    ###y~z, y simulated by z
  }
  id_gene<-rep(zindex,xnum)
  p_t_pl<-1
  for (j in 1:i){
    p_t_pl<-p_t_pl*(1-(power_yxz[j]*power_x1[j]))
  }
  p_t_pl<-1-p_t_pl
  p_t_ca<-1
  for (j in 1:i){
    p_t_ca<-p_t_ca*(1-(power_yzz[j]*power_x1[j]))
  }
  p_t_ca<-1-p_t_ca
  zx<-data.frame(id_gene=id_gene,id_snp=xgroupindex,power_zx=power_zx)
  yxz<-data.frame(id_gene=zindex,n_snp=xnum,power_yxz=power_yxz)
  powert_pl<-data.frame(n_gene=i,Genetic_Model="Pleiotropy",Estimate_Model="emGWAS",power=p_t_pl)
  yzz<-data.frame(id_gene=zindex,n_snp=xnum,power_yzz=power_yzz)
  powert_ca<-data.frame(n_gene=i,Genetic_Model="Causality",Estimate_Model="emGWAS",power=p_t_ca)
  simresult<-list(emGWAS_Pt=powert_pl,emGWAS_Ct=powert_ca,emGWA_zx=zx,emGWAS_yxz=yxz,emGWAS_yzz=yzz)
  return(simresult)
}

twasmodel<-function(z_hat,y_x,y_z,D_yx,D_yz,i,zindex,ql_z,nind,ela_net_result,xnum){
  ##model 3 (twas)
  power_yx_zhat<-NULL
  power_yz_zhat<-NULL
  for (j in 1:i){
    temp_ncp_yx<-ncp(z_hat[,j],y_x,D_yx,nind)
    temp_power_yx<-1-pt(q=ql_z,df=nind-2,ncp=temp_ncp_yx)
    power_yx_zhat<-c(power_yx_zhat,temp_power_yx)
  }
  p_yxz<-1
  for (j in 1:i){
    if (is.nan(power_yx_zhat[j])==TRUE)
      p_yxz<-p_yxz*1
    else
      p_yxz<-p_yxz*(1-power_yx_zhat[j])
  }
  p_yxz<-1-p_yxz
  ##y~zhat, y simulated by x
  for (k in 1:i){
    temp_ncp_yz<-ncp(z_hat[,k],y_z,D_yz,nind)
    temp_power_yz<-1-pt(q=ql_z,df=nind-2,ncp=temp_ncp_yz)
    power_yz_zhat<-c(power_yz_zhat,temp_power_yz)
  }
  p_yzz<-1
  for (j in 1:i){
    if (is.nan(power_yz_zhat[j])==TRUE)
      p_yzz<-p_yzz*1
    else
      p_yzz<-p_yzz*(1-power_yz_zhat[j])
  }
  p_yzz<-1-p_yzz
  ##y~zhat, y simulated by z
  yxz<-data.frame(id_gene=zindex,n_snp=xnum,power_yzhat=power_yx_zhat)
  yzz<-data.frame(id_gene=zindex,n_snp=xnum,power_yzhat=power_yz_zhat)
  pl<-data.frame(n_gene=i,Genetic_Model="Pleiotropy",Estimate_Model="TWAS",power=p_yxz)
  ca<-data.frame(n_gene=i,Genetic_Model="Causality",Estimate_Model="TWAS",power=p_yzz)
  simresult<-list(TWAS_Pt=pl,TWAS_Ct=ca,TWAS_yx_zhat=yxz,TWAS_yz_zhat=yzz,ela_net_result=ela_net_result)
  return(simresult)
}

name<-c(1:100)
q<-100
pl<-NULL
ca<-NULL
R2<-NULL
ptm<-proc.time()
for (j in 1:q){
  setwd(paste0(pa,name[j])) ## pa is the parents path, e.g. "/gpfs/qlong/GROUP_DATA/Genomes/Human_GWAS/ncbi/dbGaP-11604-CVD/20190807_qing_download/70297/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU/4Bowei/1000G/power_result/"
  sim<-simulatey(j,x,x_id,zmap,i,min,max,ht,h1,h2,kxvector,kxs,kzvector,kzs,nind)
  y_x<-sim$y_x
  y_z<-sim$y_z
  sim_x<-sim$x## all snps that were selected regarding different genes
  sim_z<-sim$z## i cols
  zindex<-sim$zindex
  z_hat<-sim$z_hat
  z_num<-sim$z_num
  xnum<-sim$xnum
  xgroupindex<-sim$xgroupindex
  D_yxx<-sim$D_yxx
  D_yxz<-sim$D_yxz
  D_yzx<-sim$D_yzx
  D_yzz<-sim$D_yzz
  ela_net_result<-sim$ela_net_result
  options(scipen = 200,digits=4)
  g<-gwasmodel(sim_x,y_x,y_z,D_yxx,D_yzx,xnum,xgroupindex,ql_x,nind,i) 
  e<-egwasmodel(sim_x,sim_z,y_x,y_z,D_yxz,D_yzz,xnum,i,xgroupindex,zindex,zmap,z_num,ql_z,nind)
  t<-twasmodel(z_hat,y_x,y_z,D_yxz,D_yzz,i,zindex,ql_z,nind,ela_net_result,xnum)
  for (j in 1:i){
    R2<-c(R2,t$ela_net_result[1,3*j])
  }
  pl<-rbind(pl,g$GWAS_Pt,e$EGWAS_Pt,t$TWAS_Pt)
  ca<-rbind(ca,g$GWAS_Ct,e$EGWAS_Ct,t$TWAS_Ct)
}
write.table(pl,"pl.csv")
write.table(ca,"ca.csv")
write.table(R2,"R2.csv")
