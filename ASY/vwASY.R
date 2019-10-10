# Function: vwASY

# VWASI(IN_NII,X_MIDLINE)
# This function computes the voxelwise asymmetry index between brain hemispheres.
# A symmetrical corregistration is assumed.
#
# -------------------------------------------------------------------------
# Inputs:
# IN_NII - Input variable of interest volume (3D NIfTI file).
# X_MIDLINE - Numeric vector indicating position of the midline in the
#             x-axis.
# ASI - Logical input. If TRUE, asymmetry index is computed too.
# -------------------------------------------------------------------------

vwASY <- function(in_nii,x_midline=23:24,ASI=T) {
  
  # 1) Read 3D Graph Theory NIfTI.
  if(require(oro.nifti)==0){
    install.packages("oro.nifti")
    library(oro.nifti)
  }
  
  in_img <- readNIfTI(in_nii)
  dim_nii <- dim(in_img@.Data)
  
  # 2) Take each hemisphere and flip one
  # Check x-axis length of each hemisphere from midline
  # and take the shortest as reference.
  mid_l <- length(x_midline)
  x_hem <- dim_nii[1]-x_midline[mid_l]
  if(x_hem > x_midline[1]-1) x_hem <- x_midline[1]-1
  
  # Take each hemisphere.
  hem_der <- in_img@.Data[(x_midline[1]-x_hem):(x_midline[1]-1),,]
  hem_izq_flip <- in_img@.Data[(x_midline[mid_l]+x_hem):(x_midline[mid_l]+1),,]
  dim_hem <- dim(hem_der)
  
  # Convert each hemisfere into a vector
  hem_der <- c(hem_der)
  hem_izq_flip <- c(hem_izq_flip)
  
  # 3) Compute asymmetry indices (AsI)
  nz_idx <- which((hem_der+hem_izq_flip)!=0)
  AsI_hem <- vector("numeric",prod(dim_hem))
  AsI_hem[nz_idx] <- (hem_der[nz_idx]-hem_izq_flip[nz_idx])
  
  # 4) Create and write 3D AsI NIfTI file
  AsI_hem <- array(AsI_hem,dim=dim_hem)
  
  # Create 3D AsI array
  in_img@.Data[x_midline,,] <- 0
  in_img@.Data[(x_midline[1]-x_hem):(x_midline[1]-1),,] <- AsI_hem
  in_img@.Data[(x_midline[mid_l]+x_hem):(x_midline[mid_l]+1),,] <- AsI_hem
  
  # Write 3D AsI NIfTI files
  out_file <- unlist(strsplit(in_nii,"[.]"))[1]
  writeNIfTI(in_img,paste0(out_file,"_ASY"))
  
  # Compute Asymmetry Index (ASI) on demand
  if(ASI){
    
    # Compute asymmetry index (standardtize asymmetry substraction)
    AsI_hem <- vector("numeric",prod(dim_hem))
    AsI_hem[nz_idx] <- (hem_der[nz_idx]-hem_izq_flip[nz_idx])/(hem_der[nz_idx]+hem_izq_flip[nz_idx])
    
    # 4) Create and write 3D AsI NIfTI file
    AsI_hem <- array(AsI_hem,dim=dim_hem)
    
    # Create 3D AsI array
    in_img@.Data[x_midline,,] <- 0
    in_img@.Data[(x_midline[1]-x_hem):(x_midline[1]-1),,] <- AsI_hem
    in_img@.Data[(x_midline[mid_l]+x_hem):(x_midline[mid_l]+1),,] <- AsI_hem
    
    # Write 3D AsI NIfTI files
    out_file <- unlist(strsplit(in_nii,"[.]"))[1]
    writeNIfTI(in_img,paste0(out_file,"_ASI"))
    
    # Apply Fisher's r-to-z transfrom
    fisher.r2z <- function(r) { 0.5 * (log(1+r) - log(1-r)) }
    # Apply r2z transformation
    in_img@.Data <- fisher.r2z(in_img@.Data)
    # Substitute infinite values for the maximum finite number
    if(sum(in_img@.Data==Inf)>0) in_img@.Data[which(in_img@.Data==Inf)] <- max(in_img@.Data[which(in_img@.Data!=Inf)])
    if(sum(in_img@.Data==-Inf)>0) in_img@.Data[which(in_img@.Data==-Inf)] <- min(in_img@.Data[which(in_img@.Data!=-Inf)])
    # Write results
    writeNIfTI(in_img,paste0(out_file,"_AsI_r2z"))
  }
}


# Creator: Zeus Chiripa
# Date: 10 Oct 2019
# Contact: zgtabuenga@comunidad.unam.mx