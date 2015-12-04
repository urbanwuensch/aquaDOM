%% aquaDOM tutorial Version 1
% To run this tutorial, execute the entire script.
% NOTE: Please load the file 'tutorial.mat', in order for the script to run properly
clc
disp('This is the tutorial for aquaDOM. Whenever the tutorial pauses (indicated by ''... ...''), press any key to continue... ...');pause;
disp(' ');
disp('Starting point for aquaDOM is a fully corrected EEM and absorbance dataset in the drEEM or DOMfluor format');
disp('For the tutorial, we are going to work with a dataset containing 5 measurements for Quinine sulfate, Salicylic acid, and Suwannee River NOM. ... ...');pause;
disp('The 5 measurements of each sample were replicates at the same concentration.')
disp(' ')
disp('This setup is chosen to perform quantum yield calculations via the zero-intercept approach.')
disp('This is a faster method, particularly suitable for CDOM quantum yields.')
disp('However, the calculation of molar fluorescence and absorbance is not possible with this data.')
disp('For the calculation of these parameters, aquaDOM requires a e.g. 5-point dilution series.')
disp('We would also choose the variable-intercept approach in this scenario.')
disp('Now back to the demo... ...'); pause();
disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp(EEMcor.filelist)
disp('Take a look at the filenames. They have been chosen to be unique between different samples, ');
disp('while allowing to distinguish between replicates of the same sampe with ''.1'', ''.2'', etc. ... ...');pause;
disp(' ');
disp('Since this dataset is fully corrected, we can proceed with the aquaDOM toolbox.')
disp(' ');
disp(' ');
disp('Next: To split the dataset, aquaDOM needs information on ''tags'', i.e. a sequence of letters in filenames, that separates samples from one another... ...');pause;
disp('This information is provided by the user in the variable tagsStandards, or tagsSamples... ...');pause;
disp(' ');
disp(' ');
disp('Next: aquaDOM needs additional data about the experiment, assembleMetadata.m guides the user step by step through this proccess,')
disp('the necessary information is provided in the comments below (lines 23-43, copy & paste the numbers)... ... ');pause

%% USER INPUT
% Standards:
% Info for Standard 1:      Qunine sulfate (NIST SRM 936a):
% Tag:                      QS 
% Reference quantum yield:  0.51
% Reference wavelength:     350
% Number of measurements:   5
% Concentrations:           [1.38E-06 1.38E-06 1.38E-06 1.38E-06 1.38E-06]
%
%
% Info for Stanard 2:       Salicylic acid (CAS: 69-72-7, AppliChem > 99.5% purity)
% Tag:                      SA
% Reference quantum yield:  0.36
% Reference wavelength:     295
% Number of measurements:   5
% %Concentrations:          [1.47E-05 1.47E-05 1.47E-05 1.47E-05 1.47E-05]
tagsStandards={'QS' 'SA'};
metadataStandards=assembleMetadata(tagsStandards,[325 600],'STD');

% Provide details for the samples here
% Info on Sample:           Suwannee River XAD-8 Fulvic acid dissolved in ultrapure water
% Tag:                      SwR
% Number of measurements:   5
% Concentration:            [8.3E-06 8.3E-06 8.3E-06 8.3E-06 8.3E-06]
tagsSamples={'SwR'}
metadataSamples=assembleMetadata(tagsSamples,[240 600],'SAM');

disp(' ');
disp(' ');
disp('Next: The EEM dataset has to be devided into quantum yield standards and samples. splitDS.m takes care of that. ... ...');pause()

% Assignment of sample IDs, and dataset split
[qy_R]=splitDS(EEMcor,metadataStandards);
[qy_X]=splitDS(EEMcor,metadataSamples);

disp('Next: The daily cross calibration. The script calculates qunatum yields based on both standards and tells the user if quality standards are met.');
disp('It also produces plots to provide a visual representation for the user to judge... ... ');pause

%% Quantum yield calculations
% Cross calibration of the standards
[qy_R]=aqy(qy_R,1,2);
disp(' ');
disp(' ');
disp(' ');
disp('The cross calibration was passed sucessfully. If not, the measurements should be repeated with freshly prepared solutions. ');
disp('Next: The calculation of apparent quantum yields of a Suwannee River XAD-8 sample... ... ');pause;close all

% Sample
[qy_X]=aqy(qy_X,2,2,qy_R,1);
disp(' ');
disp(' ');
disp(' ');
disp('All done. Now, the dataset should be reformated to be compatible with drEEM again using mergeDS.m... ... ');pause; close all
EEMcor=mergeDS(qy_X);

disp('Next, the user can inspect quantum yields using the aqyview.m script... ...');pause;
aqyview(EEMcor,[],'SUP');
disp(' ');
disp(' ');
disp(' ');
disp('Done. You successfully completed the tutorial for aquaDOM. If you need help with any of the functions, type ''help function.m'' ');