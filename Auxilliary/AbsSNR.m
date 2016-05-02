function [noise]=AbsSNR(DS)
% %Calculation of absorbance noise
% USEAGE:
%           [noise] = AbsSNR(DS)
%
% INPUTS
%               DS: Data structure containing the field Abs
%
% OUTPUTS
%            noise: Absorbance noise calculated using the range of 2nd derivative of
%                   the 10 highest absorbance measurements (CDOM will prodce a flat 2nd derivative at high wavelengths)
%             
%
% Notice:
% This mfile is part of the aquaDOM toolbox.
%
% (c) Urban Wuensch, DTU AQUA 06/2015
%
% AbsSNR.m: Copyright (C) 2015 Urban J Wünsch
% Technical University of Demark
% National Institute of Aquatic Resources
% Section for Marine Ecology and Oceanography
% Kavalergården 6
% 2920 Charlottenlund, Denmark
% urbw@aqua.dtu.dk
%
% Version 1, July 2015 First Release
AbsMax=max(DS.Abs_wave);
AbsMin=AbsMax-20;
AbsMax=knnsearch(DS.Abs_wave.',AbsMax);
AbsMin=knnsearch(DS.Abs_wave.',AbsMin);

noise=zeros(DS.nSample,1);
%% Absorbance noise calculation
for n=1:size(DS.Abs,1)
    % SNR as range of the 2nd derivative of the spectrum
    noise(n,1)=range(diff(DS.Abs(n,AbsMin:AbsMax),2),2);
end
end