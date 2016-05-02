function [DS]=SPCalib(DS,mode,method,qy_Std,qy_std_choice)
% Calibration for quantum yield calculations bypassing the cross calibration
% This method should only be used when measuring certified NIST fluorescence
% standards such as quinine sulfate (NIST 936a) in sealed cuvettes.
% These are not expected to change their properties over time and therefore
% do not strictly need to be cross calibrated daily.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE AUTHORS DO NOT RECOMMEND THE USE OF THIS FUNCTION, AS IT BYBASSES QUALITY
% CONTROL MEASURES. YOUR DATA MIGHT CONTAIN SIGNIFICANT ERRORS WITHOUT DAILY
% CROSS CALIBRATION.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% USEAGE:
%           [DS]=aqy(DS)
%
% INPUTS
%               DS: Structure containing variables .X, .Abs, and processed by splitDS.m
%
% OUTPUTS
%               DS: original structure, supplemented with quantum yields in variable .qy, and quantum yield precision in .qy_precision
%             
% Examples:
%       [qyCout]=aqy(qyCin); Calibration, variable intercept approach
%
% Notice:
% This mfile is part of the aquaDOM toolbox.
%
%
% aqy.m: Copyright (C) 2016 Urban J Wünsch
% Technical University of Demark
% National Institute of Aquatic Resources
% Section for Marine Ecology and Oceanography
% Kavalergården 6
% 2920 Charlottenlund, Denmark
% urbw@aqua.dtu.dk
%
% Version 1, May 2016 First Release

% Welcome messages
disp(' ');close all;
disp('Single point calibration.')
warning('While this method may work for you, we strongly advice against skipping the cross calibration. ')
disp('For further information read the funtion help and doi: 10.3389/fmars.2015.00098')
disp('Press any key to proceed...')
pause()

%% Initialization of new variables / functions / thresholds
% These tresholds make sure that the calculated quantum yields are of satisfying quality
DS.Abs_s_n_treshold=25; 
Rcutoff=0.95;
% Plotting function for confidence interval




%% Allocation of variables
DS.F_a=NaN(DS.nSample,DS.nEx);
DS.qy_precision=NaN(DS.nSample,DS.nEx);
DS.Abs_sqE=NaN(DS.nSample,DS.nEx);
DS.AUC_sqE=NaN(DS.nSample,DS.nEx);


tic; 

%% Calculations for QY
% Calculation of the AUC
disp('CALCULATIONS:')
disp('F'' emission integration...')
for n=1:size(DS.id,1)
    % Range for the integration of emission integrals
    min_em=knnsearch(DS.Em,DS.Em_range(n,1)); max_em=knnsearch(DS.Em,DS.Em_range(n,2));
    for i=1:DS.nConc(n,1)
        for m=1:DS.nEx
            AUC(n,i,m) = nansum(DS.X(n,i,min_em:max_em,m));
        end
    end
end

disp('F''/Abs ratio calculation...')
for n=1:DS.nSample
    for m=1:DS.nEx

        % temporary allocation of AUC and Abs data to local variable for better overview
        tempAUC=(AUC(n,1:DS.nConc(n),m));
        tempAbs=(DS.Abs(n,1:DS.nConc(n),m));
        % Calculation of mean AUC and Abs without zeros
        tempAUC_mean=mean(tempAUC(tempAUC~=0));
        tempAbs_mean=mean(tempAbs(tempAbs~=0));
        % S/N filter
        if (tempAbs/DS.noise(n,1))>(DS.Abs_s_n_treshold*DS.noise(n,1))
            DS.F_a(n,m)= tempAUC_mean/tempAbs_mean;
        end
        if DS.nConc(n)>1
            DS.AUC_sqE(n,m)=(range(tempAUC,2)/tempAUC_mean)^2;
            DS.Abs_sqE(n,m)=(range(tempAbs,2)/tempAbs_mean)^2;
        end
        clearvars tempAUC tempAbs tempAUC_mean tempAbs_mean
    end
end

%% Quantum yield calculations
%% Cross-calibration using two qy standards
disp('Calibration in progress...')
disp(' ')
disp(' ')
disp('RESULTS:')
format shortEng
disp(['Slope at reference wavelength: ',char(num2str(DS.F_a(1,knnsearch(DS.Abs_wave',DS.ref_wave))))])

% Determination of the reference excitation position (or closest match)
for n=1:DS.nSample
    DS.posRef_wave(n)=knnsearch(DS.Ex,DS.ref_wave(1,n));
end


   
    
disp(' ')
toc
end