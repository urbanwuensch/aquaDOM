cd('Data directory')
load(Xin.mat)
%% USER INPUT
% Standards:
% Info on Standard 1: QY: 0.xx at xxx nm
%                [conc1 conc2 conc3 conc4 conc5]
% Info on Stabdard 2: QY: 0.XX at xxx nm
%                [conc1 conc2 conc3 conc4 conc5]
tagsStandards={'Standard 1' 'Standard 2'};
metadataStandards=assembleMetadata(tagsStandards,[240 600],'STD');

% Provide details for the samples here
tagsSamples={'Sample 1' 'Sample 2'}
metadataSamples=assembleMetadata(tagsSamples,[240 600],'SAM');

% Assignment of sample IDs, and dataset split
[qyr]=splitDS(EEMcor,metadataStandards);
[qys]=splitDS(EEMcor,metadataSamples);

%% Quantum yield calculations
% Cross calibration of the standards
[qyr]=aqy(qyr,1,2);

% Sample
[qys]=aqy(qys,2,2,qyr,1);
[qys]=epsilon(qys);

Xout=mergeDS(qys);
aqyview(Xout,[],'SUP')

%%% END OF QUANTUM YIELD CALCULATIONS %%%