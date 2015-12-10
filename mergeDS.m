function [ DS ] = mergeDS(varargin)
% merge and transform multiple (from different measurement days), or transform single datasets from the aquaDOM variable format to drEEM compartible format
%   
% USEAGE:
%           [DS]=mergeDS(varargin)
%
% INPUTS
%               varargin: any number of structures containing .qy and .qy_precision
%
% OUTPUTS
%                 DS: merged dataset
%             
% Examples:
%       DS=mergeDS(qy_X,qy_X1,qy_X3); merging of 3 datasets into a single structure
%
% Notice:
% This mfile is part of the aquaDOM toolbox.
%
%
% mergeDS.m: Copyright (C) 2015 Urban J Wünsch
% Technical University of Demark
% National Institute of Aquatic Resources
% Section for Marine Ecology and Oceanography
% Kavalergården 6
% 2920 Charlottenlund, Denmark
% urbw@aqua.dtu.dk
%
% Version 2, December 2015 Second Release

nDS=nargin;
% New nSample
nSample=0;
for i=1:nDS
    nSample=nSample+varargin{i}.nSample;
end

DS.nSample=nSample;
DS.Abs_wave=varargin{1}.Abs_wave;
DS.nEx=varargin{1}.nEx;
DS.nEm=varargin{1}.nEm;
DS.Em=varargin{1}.Em;
DS.Ex=varargin{1}.Ex;
DS.i(:,1)=1:DS.nSample;

i=1; % i is the counter for nSample in the new dataset, m counts the old dataset (nSample)
    for n=1:nDS
        for m=1:size(varargin{n}.filelist,1)
        DS.filelist(i,1)=varargin{n}.filelist(m,1);
        if isfield(varargin{n},'ref_name') % In case user wants to use function without qy variables to merge datasets
            DS.name(1,i)=varargin{n}.name(1,m);
            DS.qy(i,:)=varargin{n}.qy(m,:);
            DS.qy_precision(i,:)=varargin{n}.qy_precision(m,:);
        end
        DS.Abs(i,:)=squeeze(varargin{n}.Abs(m,1,:));
        DS.X(i,:,:)=squeeze(varargin{n}.X(m,1,:,:));
        if varargin{n}.IntensityUnit=='AU'
            DS.X(i,:,:)=DS.X(i,:,:)./varargin{n}.RamanArea(m,1);
            DS.IntensityUnit='RU';
        end
        DS.RamanArea(i,1)=varargin{n}.RamanArea(m,1);

        if isfield(varargin{n},'ref_name')
            DS.ref_name(i,1)=varargin{n}.ref_name(1,1);
            try DS.qy_method(i,1)={cellstr(varargin{n}.qy_method(1,:))};
            catch warning('Could not transer the aqy calculation method.')
                DS.qy_method(i,1)={'N.A.'};
            end
        else
            DS.ref_name(i,1)={'N.A.'};
            DS.qy_method(i,1)={'N.A.'};
        end
        
        if isfield(varargin{n},'epsilon')
            DS.epsilon(i,:)=varargin{n}.epsilon(m,:);
        end
        if isfield(varargin{n},'Xmf')
            DS.Xmf(i,:,:)=varargin{n}.Xmf(m,:,:);
        end
        
        i=i+1;
        end
    end

disp('Merging complete.');
disp(['Datasets merged: ', num2str(nDS)]);
disp(['New number of samples in dataset: ', num2str(nSample)]);
end
