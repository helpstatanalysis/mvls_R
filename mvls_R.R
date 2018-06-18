###########################################
###########################################
#NEWS

#VERSION 0.3 RELEASED
#-NEW FUNCTION:
# -mvls.Mclust
# -mvls.pam
# -mvls.pvclust

#VERSION 0.2 RELEASED
#-BUG RESOLVED (attached library)
#-FUNCTION VISUALDIAGMVLS PERMIT TO PLOT, NOT ONLY ALL CLUSTER TOGETHER, BUT ALSO SINGLE PLOT
#-NEW MANUALS
#-NEW EXAMPLES

#THERE IS AN OPEN PROJECT AT https://github.com/helpstatanalysis/mvls_code/projects/1

###########################################
###########################################
#INSTALL PACKAGE FOR GIT-HUB (v 0.3)

library(devtools)
install_github("helpstatanalysis/mvls")

###########################################
###########################################

### Library ###

library(randomForestSRC)
library(ggplot2)
library(mice)
library(zoo)
library(reshape2)

#Function included in the final functions 

### Exclude no long data ###
#This function permit to delete rows with any values

exclude<-function(data){
  db.null<-data
  db.exclude<-data
  for(i in 1:dim(data)[1]){
    for(l in 1:dim(data)[2]){
      if(is.na(data[i,l])==T){db.null[i,l]<-0
      }else(db.null[i,l]<-1)
    }
  }
  for(i in 1:dim(data)[1]){
    if(sum(db.null[i,])==0){db.exclude<-data[-i,]}
  }
  return(db.exclude)
}

### Normalize function ###
#This function permit to normalize data (each row for the higher value f the row

normalize <- function(data){
  index<-NULL
  for (i in 1:dim(data)[1]){
    index[i]<-max(data[i,], na.rm = T)
    data[i,]<-data[i,]/max(data[i,], na.rm = T)
  }
  result<-return(list(data=data,index=index))
}

### Variation matrix ###
#This function permit to create the matrix of variation that it's used to cluster data

var.matrix<-function(data,d=.1){
  matrix.1<-data.frame(NULL)
  for (z in 1:(dim(data)[2]-1)){
    for (i in 1:dim(data)[1]){
      if(is.na(data[i,z])==T){
        matrix.1[i,z]<-3
      } else if(is.na(data[i,z+1])==T){
        matrix.1[i,z]<-3
      } else if(data[i,z+1]>(data[i,z]+d)){
        matrix.1[i,z]<-1
      } else if(data[i,z+1]<(data[i,z]-d)){
        matrix.1[i,z]<-2
      } else{matrix.1[i,z]<-0
      }
    }
  }
  matrix.2<-data.frame(NULL)
  for (z in 1:(dim(data)[2]-2)){
    for (i in 1:dim(data)[1]){
      if(is.na(data[i,z])==T){
        matrix.2[i,z]<-3
      } else if(is.na(data[i,z+2])==T){
        matrix.2[i,z]<-3
      } else if(data[i,z+2]>(data[i,z]+d)){
        matrix.2[i,z]<-1
      } else if(data[i,z+2]<(data[i,z]-d)){
        matrix.2[i,z]<-2
      } else{matrix.2[i,z]<-0
      }
    }
  }
  matrix.3<-data.frame(NULL)
  for (i in 1:dim(data)[1]){
    if(is.na(data[i,1])==T){
      matrix.3[i,1]<-3
    } else if(is.na(data[i,dim(data)[2]])==T){
      matrix.3[i,1]<-3
    } else if(data[i,dim(data)[2]]>(data[i,1]+d)){
      matrix.3[i,1]<-1
    } else if(data[i,dim(data)[2]]<(data[i,1]-d)){
      matrix.3[i,1]<-2
    } else{matrix.3[i,1]<-0
    }
  }
  matrix<-cbind(matrix.1,matrix.2,matrix.3)
  return(matrix)
}

###  Pre-impute data ###
#This function permit to pre-impute data to create variable matrix when all missing value rows are grouped into the same cluster

preimputation<-function(data, imp.method='mean'){
  if(imp.method=='mean'){
    data<-mice(data, method = 'mean',printFlag = F)
    datapreimput<-complete(data)
  }
  if(imp.method=='locf'){
    data<-t(db.prov)
    data<-na.locf(data)
    data<-t(data)
    data<-mice(data, method = 'mean', printFlag = F)
    datapreimput<-complete(data)
  }
  return(datapreimput)
}

### Function toclusterfunc ###
#These function includes other function to simplify end funcion code

toclusterfunc.noimp<-function(data, d){
  data<-exclude(data)
  matrix<-data
  result<-normalize(data)
  data<-result$data
  index<-result$index
  db.var<-var.matrix(data,d)
  result<-return(list(data=data,index=index,vari.matrix=db.var,matrix=matrix))
}

toclusterfunc.imp<-function(data, d, imp.method='mean'){
  data<-exclude(data)
  matrix<-data
  data.imp<-preimputation(data, imp.method)
  result.imp<-normalize(data)
  db.var<-var.matrix(result.imp$data,d)
  result<-normalize(data)
  data<-result$data
  index<-result$index
  result<-return(list(data=data,index=index,vari.matrix=db.var,matrix=matrix))
}

### Imputation of missing data ###
#This function permit to impute missing value using cluster from var.matrix builted on pre-imputed or not data

imputation.kh<-function(data, vari.matrix, method='k', cluster=6, nstart=20){
  if(method=='k'){
    clusters <- kmeans(vari.matrix, cluster, nstart)
    clusterCut <- clusters$cluster
  }else if(method=='h'){
    clusters <- hclust(dist(vari.matrix))
    clusterCut <- cutree(clusters, cluster) 
  }
  clu.matrix<-matrix(rep(clusterCut,dim(data)[2]),ncol =(dim(data)[2]))
  sd.1.j<-matrix(c(rep(NA,(dim(data)[2]*dim(data)[1]))),ncol=dim(data)[2],nrow=dim(data)[1])
  sd.1.k<-matrix(c(rep(NA,(dim(data)[2]*dim(data)[1]))),ncol=dim(data)[2],nrow=dim(data)[1])
  for(i in 1:dim(data)[2]){
    for(l in 1:dim(data)[1]){
      if(is.na(data[l,i])==T){
        if(i<(dim(data)[2]-1)){ 
          if(is.na(data[l,(i+1)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i+1)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i+1)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i+1)]+(j-k))>1){h<-1
              }else if((data[l,(i+1)]+(j-k))<0){h<-0
              }else(h<-data[l,(i+1)]+(j-k))
              data[l,i]<-h
            }
          }
          if(is.na(data[l,(i+1)])==T && is.na(data[l,(i+2)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i+2)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i+2)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i+2)]+(j-k))>1){h<-1
              }else if((data[l,(i+2)]+(j-k))<0){h<-0
              }else(h<-data[l,(i+2)]+(j-k))
              data[l,i]<-h
            }
          }
        }
        if(i==(dim(data)[2]-1)){ 
          if(is.na(data[l,(i-1)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i-1)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i-1)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i-1)]-(j-k))>1){h<-1
              }else if((data[l,(i-1)]-(j-k))<0){h<-0
              }else(h<-data[l,(i-1)]-(j-k))
              data[l,i]<-h
            }
          }
          if(is.na(data[l,(i-1)])==T && is.na(data[l,(i-2)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i-2)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i-2)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i-2)]-(j-k))>1){h<-1
              }else if((data[l,(i-2)]-(j-k))<0){h<-0
              }else(h<-data[l,(i-2)]-(j-k))
              data[l,i]<-h
            }
          }
        }
        if(i==(dim(data)[2])){ 
          if(is.na(data[l,(i-1)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i-1)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i-1)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i-1)]-(j-k))>1){h<-1
              }else if((data[l,(i-1)]-(j-k))<0){h<-0
              }else(h<-data[l,(i-1)]-(j-k))
              data[l,i]<-h
            }
          }
          if(is.na(data[l,(i-1)])==T && is.na(data[l,(i-2)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i-2)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i-2)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i-2)]-(j-k))>1){h<-1
              }else if((data[l,(i-2)]-(j-k))<0){h<-0
              }else(h<-data[l,(i-2)]-(j-k))
              data[l,i]<-h
            }
          }
        }
      } 
    }
  }
  
  for(i in 1:dim(data)[2]){
    for(l in 1:dim(data)[1]){
      if(is.na(data[l,i])==T){
        if(i<(dim(data)[2]-1)){ 
          if(is.na(data[l,(i+1)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i+1)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i+1)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i+1)]+(j-k))>1){h<-1
              }else if((data[l,(i+1)]+(j-k))<0){h<-0
              }else(h<-data[l,(i+1)]+(j-k))
              data[l,i]<-h
            }
          }
          if(is.na(data[l,(i+1)])==T && is.na(data[l,(i+2)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i+2)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i+2)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i+2)]+(j-k))>1){h<-1
              }else if((data[l,(i+2)]+(j-k))<0){h<-0
              }else(h<-data[l,(i+2)]+(j-k))
              data[l,i]<-h
            }
          }
        }
        if(i==(dim(data)[2]-1)){ 
          if(is.na(data[l,(i-1)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i-1)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i-1)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i-1)]-(j-k))>1){h<-1
              }else if((data[l,(i-1)]-(j-k))<0){h<-0
              }else(h<-data[l,(i-1)]-(j-k))
              data[l,i]<-h
            }
          }
          if(is.na(data[l,(i-1)])==T && is.na(data[l,(i-2)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i-2)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i-2)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i-2)]-(j-k))>1){h<-1
              }else if((data[l,(i-2)]-(j-k))<0){h<-0
              }else(h<-data[l,(i-2)]-(j-k))
              data[l,i]<-h
            }
          }
        }
        if(i==(dim(data)[2])){ 
          if(is.na(data[l,(i-1)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i-1)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i-1)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i-1)]-(j-k))>1){h<-1
              }else if((data[l,(i-1)]-(j-k))<0){h<-0
              }else(h<-data[l,(i-1)]-(j-k))
              data[l,i]<-h
            }
          }
          if(is.na(data[l,(i-1)])==T && is.na(data[l,(i-2)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i-2)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i-2)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i-2)]-(j-k))>1){h<-1
              }else if((data[l,(i-2)]-(j-k))<0){h<-0
              }else(h<-data[l,(i-2)]-(j-k))
              data[l,i]<-h
            }
          }
        }
      } 
    }
  }
  
  result<-list(data=data, clu.matrix=clu.matrix, sd.1.j=sd.1.j, sd.1.k=sd.1.k)
}

    
#Final MVLS fuction

### MVLS.PRINT ###
# This function permit to identify best k values using some plot (different for k and h method) using original dataset or var.matrix
    
mvls.print<-function(data, d, method='h', varmatrix=F, kmax=20){
  if(method=='h' && varmatrix==F){
    data<-exclude(data)
    clusters <- hclust(dist(data))
    return(plot(clusters))
  } else if(method=='h' && varmatrix==T){
    data<-toclusterfunc.noimp(data,d)$vari.matrix
    clusters <- hclust(dist(data))
    return(plot(clusters))
  }else if (method=='k' && varmatrix==F){
    wss <- sapply(1:kmax, function(k){kmeans(na.omit(data), k, nstart=50,iter.max = 15 )$tot.withinss})
    plot(1:kmax, wss, type="b", pch = 19, frame = FALSE, 
         xlab="Number of clusters K",
         ylab="Total within-clusters sum of squares")
  } else if (method=='k' && varmatrix==T){
    data<-toclusterfunc.noimp(data,d)$vari.matrix
    wss <- sapply(1:kmax, function(k){kmeans(na.omit(data), k, nstart=50,iter.max = 15 )$tot.withinss})
    plot(1:kmax, wss, type="b", pch = 19, frame = FALSE, 
         xlab="Number of clusters K",
         ylab="Total within-clusters sum of squares")
  }else{cat('Error: only k method and h method permitted')}
}

### MVLS ###
#This function permit to impute a single dataset and to obtain sd used to this imputation (sd.1)
    
mvls<-function(data, d=0.1, method='k', cluster=6, nstart=10, pre.imp=F, imp.method='mean'){
  if(pre.imp==T){
    results<-toclusterfunc.imp(data,d,imp.method)
  }else if(pre.imp==F){
    results<-toclusterfunc.noimp(data,d)
  }
  data.f<-results$data
  index.f<-results$index
  vari.matrix.f<-results$vari.matrix
  result<-imputation.kh(data.f, vari.matrix.f, method, cluster, nstart)
  db<-result$data*index.f
  sd.1.j<-result$sd.1.j*index.f
  sd.1.k=result$sd.1.k*index.f
  return(list(data=db, data.norm=result$data, cluster=result$clu.matrix, matrix=results$matrix, sd.1.j=sd.1.j, sd.1.k=sd.1.k))
}

### MVLSBOOT ###
#This function permit to relize multiple imputation (n=5,10,15) and to obtain sd used for missing value imputation (sd.2)

mvlsboot<-function(data, d=0.1, method='k', cluster=12, nstart=20, imp=F, imp.method='mean', boot='medium'){
  sd.2<-matrix(c(rep(NA,(dim(data)[2]*dim(data)[1]))),ncol=dim(data)[2],nrow=dim(data)[1])
  result.boot<-data
  if(boot=='low'){
    db.boot.1<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.2<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.3<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.4<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.5<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    for (i in 1:dim(result.boot)[2]){
      for (l in 1:dim(result.boot)[1]) {
        sd.2[l,i]<-rnorm(1,mean=0,sd=sd(c(db.boot.1[l,i],db.boot.2[l,i],db.boot.3[l,i],db.boot.4[l,i],db.boot.5[l,i]),na.rm = T))
        result.boot[l,i]<-mean(c(db.boot.1[l,i],db.boot.2[l,i],db.boot.3[l,i],db.boot.4[l,i],db.boot.5[l,i]),na.rm = T)+sd.2[l,i]
      }
    }
  }
  if(boot=='medium'){
    db.boot.1<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.2<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.3<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.4<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.5<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.6<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.7<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.8<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.9<-mvls(data, d, method, cluster, nstart)$data
    db.boot.10<-mvls(data, d, method, cluster, nstart)$data
    for (i in 1:dim(result.boot)[2]){
      for (l in 1:dim(result.boot)[1]) {
        sd.2[l,i]<-rnorm(1,mean=0, sd=sd(c(db.boot.1[l,i],db.boot.2[l,i],db.boot.3[l,i],db.boot.4[l,i],db.boot.5[l,i],db.boot.6[l,i],db.boot.7[l,i],db.boot.8[l,i],db.boot.9[l,i],db.boot.10[l,i]),na.rm = T))
        result.boot[l,i]<-mean(c(db.boot.1[l,i],db.boot.2[l,i],db.boot.3[l,i],db.boot.4[l,i],db.boot.5[l,i],db.boot.6[l,i],db.boot.7[l,i],db.boot.8[l,i],db.boot.9[l,i],db.boot.10[l,i]),na.rm = T)+sd.2[l,i]
      }
    }
  }
  if(boot=='high'){
    db.boot.1<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.2<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.3<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.4<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.5<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.6<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.7<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.8<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.9<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.10<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.11<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.12<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.13<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.14<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    db.boot.15<-mvls(data, d, method, cluster, nstart,imp, imp.method)$data
    for (i in 1:dim(result.boot)[2]){
      for (l in 1:dim(result.boot)[1]) {
        sd.2[l,i]<-rnorm(1,mean=0, sd=sd(c(db.boot.1[l,i],db.boot.2[l,i],db.boot.3[l,i],db.boot.4[l,i],db.boot.5[l,i],db.boot.6[l,i],db.boot.7[l,i],db.boot.8[l,i],db.boot.9[l,i],db.boot.10[l,i],db.boot.11[l,i],db.boot.12[l,i],db.boot.13[l,i],db.boot.14[l,i],db.boot.15[l,i]),na.rm = T))
        result.boot[l,i]<-mean(c(db.boot.1[l,i],db.boot.2[l,i],db.boot.3[l,i],db.boot.4[l,i],db.boot.5[l,i],db.boot.6[l,i],db.boot.7[l,i],db.boot.8[l,i],db.boot.9[l,i],db.boot.10[l,i],db.boot.11[l,i],db.boot.12[l,i],db.boot.13[l,i],db.boot.14[l,i],db.boot.15[l,i]),na.rm = T)+sd.2[l,i]
      }
    }
  }
  return(list(data=result.boot,sd.2=sd.2))
}

### VISUALDIAGMLVS ###
#This function permit to realize a plot showing cluster pattern with different color and permit to identify adeguate cluserization
#New releasing in version 0.2

visualdiagmvls<-function(mvls, method="multiple", cluster, norm=F){
  if(method=="multiple"){
    if(norm=="TRUE"){matrix<-mvls$data.norm}else if(norm=="FALSE"){matrix<-mvls$matrix}
    Cluster<-as.character(as.vector(mvls$cluster[,1]))
    id<-seq(1,dim(matrix)[1], by=1)
    matrix<-data.frame(id,matrix,Cluster)
    matrix.ggplot<-reshape(na.omit(matrix),idvar ="id", varying=list(2:(dim(matrix)[2]-1)), direction = "long")
    ggplot(data=matrix.ggplot, aes(x=time, y=a, colour=Cluster, group=id))+geom_line(alpha=.5)+ggtitle("Distribuzione dei pattern")+labs (x="Time", y = "Values")+theme_classic()
  }else if(method=="single"){
    if(norm=="TRUE"){matrix<-mvls$data.norm}else if(norm=="FALSE"){matrix<-mvls$matrix}
    Cluster<-as.vector(mvls$cluster[,1])
    data<-cbind(matrix,Cluster)
    data.1<-subset(data, data$Cluster == cluster)
    matrix.1<-data.1[,-dim(data.1)[2]]
    id<-seq(1,dim(matrix.1)[1], by=1)
    matrix.ggplot<-reshape(na.omit(matrix.1),idvar ="id", varying=list(1:(dim(matrix.1)[2])), direction = "long")
    ggplot(data=matrix.ggplot, aes(x=time, y=a, group=id))+geom_line(alpha=.5)+ggtitle(paste0("Distribuzione del pattern",cluster))+labs (x="Time", y = "Values")+theme_classic()
  }
}

### MVLS.MCLUST ###
#This function permits to use Mclust on dataset.
#New releasing in version 0.3
    
mvls.Mclust<-function(data, imp.method='mean', G=1:9){
  data<-exclude(data)
  data.pi<-preimputation(data,imp.method)
  result<-Mclust(data.pi,G)
  clusterCut<-result$classification
  plot(result, what = "classification")
  clu.matrix<-matrix(rep(clusterCut,dim(data)[2]),ncol =(dim(data)[2]))
  sd.1.j<-matrix(c(rep(NA,(dim(data)[2]*dim(data)[1]))),ncol=dim(data)[2],nrow=dim(data)[1])
  sd.1.k<-matrix(c(rep(NA,(dim(data)[2]*dim(data)[1]))),ncol=dim(data)[2],nrow=dim(data)[1])
  for(i in 1:dim(data)[2]){
    for(l in 1:dim(data)[1]){
      if(is.na(data[l,i])==T){
        #[...]
  
### MVLS.PAM ###
#This function permits to use pamk and pam on dataset.
#New releasing in version 0.3
        
mvls.pam<-function(data, imp.method='mean', krange=2:12, scaling=T){
  data<-exclude(data)
  data.pi<-preimputation(data,imp.method)
  k<-pamk(data.pi,krange, scaling)$nc
  pam<-pam(data.pi, k)
  pam$clustering
  clusterCut<-pam$clustering
  clu.matrix<-matrix(rep(clusterCut,dim(data)[2]),ncol=(dim(data)[2]))
  sd.1.j<-matrix(c(rep(NA,(dim(data)[2]*dim(data)[1]))),ncol=dim(data)[2],nrow=dim(data)[1])
  sd.1.k<-matrix(c(rep(NA,(dim(data)[2]*dim(data)[1]))),ncol=dim(data)[2],nrow=dim(data)[1])
  for(i in 1:dim(data)[2]){
    for(l in 1:dim(data)[1]){
      if(is.na(data[l,i])==T){
        #[...]
        
### MVLS.PVCLUST ###
#This function permits to test h-clustering with pvclust function.
#New releasing in version 0.3
        
mvls.pvclust<-function(mvls, data="data",nboot=1000){
  if(data=="data"){
    result<-pvclust(mvls$data, method.hclust="average",method.dist="correlation",use.cor="pairwise.complete.obs",nboot=10)
  }else if(data=="vari.matrix"){
    result<-pvclust(mvls$vari.matrix, method.hclust="average",method.dist="correlation",use.cor="pairwise.complete.obs",nboot=10)
  }
  return(result)
}

############################################
###########      EXAMPLES       ############
    
db.prov<-data.frame(sample(1:100, 300, replace=T),sample(1:100, 300, replace=T),sample(1:100, 300, replace=T),sample(1:100, 300, replace=T))
names(db.prov)<-c('a','b','c','d')

db.prov[240,1]<-NA 
db.prov[241,2]<-NA
db.prov[242,3]<-NA
db.prov[243,4]<-NA
db.prov[244,1:2]<-NA
db.prov[245,3:4]<-NA
db.prov[246,2:3]<-NA
db.prov[247,1:2]<-NA
db.prov[248,2:4]<-NA
db.prov[249,1:4]<-NA

library(mvls)

mvls.print(db.prov, d=0.1, method = "k", varmatrix = F)
mvls.print(db.prov, d=0.1, method = "k", varmatrix = T)
mvls.print(db.prov, d=0.1, method = "h", varmatrix = F)
mvls.print(db.prov, d=0.1, method = "h", varmatrix = T)

mvls<-mvls(db.prov,d=0.1,cluster = 8, method = "k")
mvls<-mvls(db.prov,d=0.1,cluster = 8, method = "k", pre.imp = T, imp.method = "locf")
mvls<-mvls(db.prov,d=0.1,cluster = 8, method = "k", pre.imp = T)

mvls<-mvls(db.prov,d=0.1,cluster = 8, method = "h")
mvls<-mvls(db.prov,d=0.1,cluster = 8, method = "h", pre.imp = T, imp.method = "locf")
mvls<-mvls(db.prov,d=0.1,cluster = 12, method = "h", pre.imp = T)

visualdiagmvls(mvls, method = "multiple", norm = T)
visualdiagmvls(mvls, method = "multiple", norm = F)
for(i in 1:12){print(visualdiagmvls(mvls, method = "single", cluster=i, norm = T))}
visualdiagmvls(mvls, method = "single", cluster=4)

mvlsboot(db.prov, d=0.1, method = "k", cluster=8, nstart = 20, pre.imp = T, imp.method = "locf", boot="high")
        
#NEW!!!
        
result<-mvls.pvclust(mvls, data="vari.matrix",nboot=10000)
plot(result)
pvrect(result,alpha = 0.95)
seplot(result)
        
result<-mvls.pam(db.prov)
result$data
plot(result$pam)
result$pam$clusinfo

result<-mvls.Mclust(db.prov,imp.method="locf",G=1:12)
result$data
result$sd.1
