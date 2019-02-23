function AAL2nii(vector_AAL,nombre)
 
% Creator: Zeus Gracia
% Date: 14 November 2015
% Function: AAL2nii
 

% AAL2NII
% This function returns a brain 3D NIfTI with the vector values to each of
% the AAL brain regions.
%
% -------------------------------------------------------------------------
% Inputs:
% VECTOR_AAL - 116's length vector with a value corresponding to each of 
% the AAL areas.
% NOMBRE - String with the file prefix
% -------------------------------------------------------------------------

% Load AAL map.
AAL = load_nii('AAL_4mm.nii.gz');
[nx,ny,nz]=size(AAL.img);
auxAAL=reshape(AAL.img,[nx*ny*nz,1]);
% Load ROI number.
ROI=tdfread('ROI_MNI_V4.txt');
    
auxnii=zeros([nx*ny*nz,1]);

for ii=1:length(vector_AAL)
    auxnii(auxAAL==ROI.Num(ii))=vector_AAL(ii);
end
nii=AAL;
nii.img=reshape(auxnii,[nx ny nz]);
% Change data type
nii.hdr.dime.bitpix=64;nii.hdr.dime.datatype=64;
nii.fileprefix=nombre;
% Save NIfTI.
save_nii(nii,[nombre,'_AAL4mm.nii.gz']);

end