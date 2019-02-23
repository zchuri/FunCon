function read_ConnMx(dirdatos,reg_exp,mapa,roifile,Tp,outfile)
 
% Creator: Zeus Gracia
% Date: 17 November 2015
% Function: read_ConnMx
 

% READ_CONNMX(DIRDATOS,REGEXP,MAPA,OUTFILE)
% Produces a correlation matrix from a set of NIfTI files. 
%
% -------------------------------------------------------------------------
% Inputs:
% DIRD - pathway of the data folder.
% REG_EXP - regular expresion to access to relevant NifTi files.
% MAPA - select the type of brain atlas. The atlas has to be in NifTi
%        format .nii.gz
% ROIFILE - text file with the atlas' ROIs name and numbers.
% Tp - refers to the time points to include that subject in the
%      further analyses.
% OUTFILE - Name of data output.
% -------------------------------------------------------------------------

% By default
% mapa='AAL_4mm.nii.gz';

% Load AAL brain atlas.
AAL = load_nii(mapa, [], 1);
AAL_map=AAL.img;
[nx,ny,nz]=size(AAL_map);
auxAAL=reshape(AAL_map,[nx*ny*nz,1]); % aux vector with AAL ROIs

% Asigns a number to each AAL ROIs.
% roifile='ROI_MNI_V4.txt';
ROI=tdfread(roifile);
ROI_num=ROI.Num;

% Read the NifTi files of interest.
% dirdatos='/misc/sherrington/zgtabuenca/Adolescentes/FunImgARCWF';
% reg_exp='/*.gz';
D=regexpdir(dirdatos,reg_exp);
if isempty(D)
    error('Empty data folder')
end

% Initialitation of the output variables.
total_ConnMx=zeros(length(ROI_num),length(ROI_num),length(D));
total_names=cell(length(D),1);
%total_index=uint16(zeros(length(D),1));

% Aux variable with the final number of subjects.
aux_ii=0;

% Read each NifTi file.
for ii=1:length(D)
    nii = load_nii(char(D(ii)), [], 1);

    % Matrix dimensions.
    [nx,ny,nz,nt]=size(nii.img);

    % Less than Tp time points then the subject is rejected.
    if nt >= 90
        % Extract the temporal series mean of each ROI 
        aux_ii=aux_ii+1;
        auxnii=reshape(nii.img,[nx*ny*nz,nt]);
        auxnii=auxnii(:,1:Tp);

        ROIst=zeros(Tp,length(ROI_num));
        for ii=1:length(ROI_num)
            ROIst(:,ii)=mean(auxnii(auxAAL==ROI_num(ii),:))';
        end

        % Correlation matrix of all ROIs.
        total_ConnMx(:,:,aux_ii)=corr(ROIst);
        total_names(aux_ii,:)=cellstr(num2str(ii));
        %total_index(aux_ii,:)=ii;
    end
end

% Remove empty space in output data (due to rejected subjects).
total_ConnMx=total_ConnMx(:,:,1:aux_ii);
total_names=total_names(1:aux_ii,:);
% total_index=total_index(1:aux_ii);

% Data overview.
%figure; imagesc(mean(total_MxCorr,3))
%title(['AAL ROIs - Correlation Matrix (n=',num2str(aux_ii),')']);

% Save output variables in an struct array.
results_ConnMx.name=total_names;
% results_ConnMx.index=total_index;
results_ConnMx.ConnMx=total_ConnMx;

if exist('ID.csv','file')
    id=csvread('ID.csv');
    results_ConnMx.index=id(1:end,2);
    results_ConnMx.etapa=id(1:end,1);
    results_ConnMx.edad=id(1:end,3);
    results_ConnMx.sexo=id(1:end,4);
end
results_ConnMx.edad(results_ConnMx.edad==0)=NaN;

save(outfile,'results_ConnMx')

