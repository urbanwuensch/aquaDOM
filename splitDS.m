function [qy]=splitDS(DS,metadata)
% split drEEM dataset based on supplied metadata
%   
% USEAGE:
%           [qy]=splitDS(DS,metadata)
%
% INPUTS
%         DS:	drEEM-compartible data structure, containing .X, and .Abs
%   metadata:	varaible containing metadata produced by assembleMetadata.m
%
% OUTPUTS
%         qy:	split dataset ready for aqy.m
%             
% Examples:
%       [qy_R]=splitDS(EEMcor,metadataStandards); split of EEMcor with metadata supplied in metadataStandards
%
% Notice:
% This mfile is part of the aquaDOM toolbox.
%
%
% splitDS.m: Copyright (C) 2015 Urban J Wünsch
% Technical University of Demark
% National Institute of Aquatic Resources
% Section for Marine Ecology and Oceanography
% Kavalergården 6
% 2920 Charlottenlund, Denmark
% urbw@aqua.dtu.dk
%
% Version 2, December 2015 Second Release


%% Initial data format check
if strcmp(formatCheck(DS),'NOT PASSED')
    error('Your Dataset is incomplete. Cannot continue. Please include above mentioned variables')
end
%% Tagging of DS files with matching ID number for identification of dilution series / replicates in dataset
DS.id=nan(DS.nSample,1);
for i=1:DS.nSample
    for n=1:size(metadata.tag,2)
        if cell2mat(strfind(DS.filelist(i,1),char(metadata.tag(1,n))));
        DS.id(i,1)=n;
        hit=1;
        end
    end
end

%% Reduction of variables to identified samples (deletion of samples that were not mentioned in metadata)
% Condition is ~isnan(DS.id)
DS.filelist=DS.filelist(~isnan(DS.id),:,:);
DS.X=DS.X(~isnan(DS.id),:,:);
DS.Abs=DS.Abs(~isnan(DS.id),:,:);

if isfield(DS,'offset')
    qy.MiscMetadata.AbsOffset=DS.offset(~isnan(DS.id),:,:);
end

if isfield(DS,'RamanArea')
    DS.RamanArea=DS.RamanArea(~isnan(DS.id),:,:);
end



% Abs SNR calculation (if not done manually before aquaDOM)
if isfield(DS,'noise')
    qy.noise=DS.noise(~isnan(DS.id),:,:);
else 
    qy.noise=AbsSNR(DS);
end
DS.id=DS.id(~isnan(DS.id),:,:);

%% Transfer Abs and Flu data into new sample x concentration x data scheme
Abs=NaN(size(metadata.tag,2),size(metadata.conc,2),size(DS.Abs_wave,2));
X=NaN(size(metadata.tag,2),size(metadata.conc,2),DS.nEm,DS.nEx);
i=1;
for n=1:size(DS.id,1)
    Abs(DS.id(n),i,:)=DS.Abs(n,:);
    X(DS.id(n),i,:,:)=DS.X(n,:,:);
    RamanArea(DS.id(n),1)=DS.RamanArea(n,1);
    % If EEMs were normalized to integration times
    if isfield(DS,'intT')
        intT(DS.id(n),i,1)=DS.intT(n,1);
    end
    filelist(DS.id(n),1)=DS.filelist(n);
    id(DS.id(n),1)=DS.id(n);
    i=i+1;
    if i>metadata.nConc(DS.id(n))
       i=1;
    end

end


%% Creation of new qy structure
% Transfer of metadata into new variable
qy.name=metadata.name;
qy.tag=metadata.tag;
qy.conc=metadata.conc;
if isfield(metadata, 'nConc')
    qy.nConc=metadata.nConc;
end
qy.Em_range=metadata.Em_range;
qy.id=id;
qy.filelist=filelist;

% Transfer of vital DS data into new structure
qy.Abs_wave=DS.Abs_wave;
qy.nSample=size(metadata.tag,2);
qy.nEx=DS.nEx;
qy.nEm=DS.nEm;
qy.Ex=DS.Ex;
qy.Em=DS.Em;


% Transfer of nice to have data into nested structure
if isfield(DS,'intT')
    qy.MiscMetadata.intT=intT;
end


if ischar(DS.Smooth)
    qy.MiscMetadata.smootheem=DS.Smooth;
end


if isfield(DS,'IntensityUnit')
    qy.IntensityUnit=DS.IntensityUnit;
end

% This is only relevant if dataset contains qy standards
if isfield(metadata,'ref')
    qy.ref=metadata.ref;
    qy.ref_wave=metadata.ref_wave;
end



% Transfer Abs and Flu data and Raman Areas
qy.Abs=Abs;
qy.X=X;
qy.RamanArea=RamanArea;

% Last: Order alphabetically for easier navigation
qy=orderfields(qy);

end