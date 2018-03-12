# Created on Mar 8th 2018

# @author: feardonn@ie.ibm.com

## The objective of this script is to create pseudo observations
## for testing of data assimilation module of a simple
## harbour configuration applied to EFDC

## The basic idea is that the EFDC model is run once WITHOUT data assimilatio
## to generate the "correct" solution
## Then it is run in data assimilation mode (with perturbed input or boundary conditions)
## to evaluate and demonstrate data assimilation capabilities


## Set the working directory to the same location as the script
script.dir <- dirname(sys.frame(1)$ofile)
setwd(script.dir)
require("ncdf4")

# Loop over the generated netCDF files and 
out_mat <- array(0, dim = c(550, 6))
for (ifil in seq(1,12)){
  fin_name = paste("efdcout_2017-01-01-",sprintf("%02d",ifil),"0000.nc",sep="")
  print(fin_name)
  ncin <- nc_open(fin_name)
  u <- ncvar_get(ncin,"u",verbose=F) # data, an array dim 20x225x455
  v <- ncvar_get(ncin,"v",verbose=F) # data, an array dim 20x225x455
  ii <- 0
  fout_name = paste("observations_2017-01-01-",sprintf("%02d",ifil),"00.csv",sep="")
  for (i in seq(3, 13)){
    for (j in seq(4, 53)){
      ii = ii + 1
      out_mat[ii, ] <- cbind(i,j, round(u[i,j,3],digits=5)*100 ,round(v[i,j,3],digits=5)*100, round(t(runif(2)),digits=3))       
    }
  }
  write.csv(out_mat, fout_name, row.names=FALSE, col.names=FALSE)
}
