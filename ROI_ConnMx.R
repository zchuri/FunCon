
# ROI_ConnMx(fMRI,ROIs,outdir)
# This function read filename fMRI and a corresponding mask with labeled ROIs
# to compute a ROI-to-ROI functional connectivity matrix. 
#
# -------------------------------------------------------------------------
# Inputs:
# fMRI - fMRI 4D NIfTI filename.
# ROIs - 3D NIfTI filename.
# outdir - output directory.
# -------------------------------------------------------------------------
#
#Zeus Chiripa
#07 Dic 2017
#zgtabuenca@gmail.com 

ROI_ConnMx <- function(fMRI,ROIs,outdir) {
  
  # Check if those inputs exists
  if(!file.exists(fMRI)) stop("STOP: input fMRI volume not found!")
  if(!file.exists(ROIs)) stop("STOP: input native atlas volume not found!")

  # First of all create a output folder
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  
  # 'oro.nifti' package is needed
  if(!require(oro.nifti)){
    install.packages("oro.nifti")
    library(oro.nifti)
  }
  
  cat("\n#####################################################")
  cat("\n1) Compute Connectivity Matrix")
  cat("\n#####################################################\n")
  
  atlas <- readNIfTI(ROIs,reorient=F)
  atlas <- atlas@.Data
  atlas_dim <- dim(atlas)
  # Convert atlas img into a vector in order to simplify the computation
  atlas <- c(atlas)
  # Get ROIs labels
  ROI_lab <- as.integer(names(table(atlas)))
  # Discard voxels with label zero
  ROI_lab <- ROI_lab[which(ROI_lab!=0)]
  ROI_n <- length(ROI_lab)
  cat(paste0("\nNumber of ROIs: ",ROI_n,"\n"))
  
  # Read fMRI NIfTI file
  fMRI_img <- readNIfTI(fMRI,reorient=F)
  fMRI_img <- fMRI_img@.Data
  fMRI_dim <- dim(fMRI_img)
  # Check if fMRI and ROIs files have the same dimensions
  if(!identical(atlas_dim,fMRI_dim[1:3])){
    stop("STOP: fMRI and ROIs volumes have different dimensions!!!")
  } 
  
  # Extract the temporal series mean of each ROI.
  cat("\nExtracting timeseries mean of each ROI\n")
  # Convert fMRI img into a matrix in order to simplify the computation
  fMRI_img <- array(fMRI_img, dim=c(prod(fMRI_dim[1:3]),fMRI_dim[4]))
  fMRI_ROI <- array(0,dim=c(ROI_n,fMRI_dim[4]))
  # Extract timeseries mean
  for(jj in 1:ROI_n){
    ROI_idx <- which(atlas==ROI[jj])
    if(length(ROI_idx)==1) {
      # The ROI only covers one voxel
      fMRI_ROI[jj,] <- fMRI_img[ROI_idx,]
    } else fMRI_ROI[jj,] <- apply(fMRI_img[ROI_idx,],2,mean)
  }
  
  # Calculate ROI's connectivity matrix (Pearson's correlation approach)  
  ConnMx <- cor(t(fMRI_ROI))
  # If NA's convert to zeros
  if(sum(is.na(ConnMx))>0) ConnMx[is.na(ConnMx)] <- 0
  
  # Write output
  cat("\nWrite connectivity matrix\n")
  rownames(ConnMx) <- as.character(ROI_lab)
  colnames(ConnMx) <- as.character(ROI_lab)
  saveRDS(ConnMx,paste0(outdir,"/ConnMx_",ROI_n,".rds"))
  write.table(ConnMx,paste0(outdir,"/ConnMx_",ROI_n,".txt"),quote=F,sep="\t")
}
