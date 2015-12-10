function [result] = formatCheck(DS)
% Check data format to ensure aquaDOM compatibility
%   
% USEAGE:
%           [result] = formatCheck(DS)
%
% INPUTS
%               DS: Structure containing fluorescence and absorbance data
%
% OUTPUTS
%               result: 'PASSED' or 'NOT PASSED'
%               Additional warning with specifics will be shown if result='NOT PASSED'
%             
%
% Notice:
% This mfile is part of the aquaDOM toolbox.
%
%
% formatCheck.m: Copyright (C) 2015 Urban J Wünsch
% Technical University of Demark
% National Institute of Aquatic Resources
% Section for Marine Ecology and Oceanography
% Kavalergården 6
% 2920 Charlottenlund, Denmark
% urbw@aqua.dtu.dk
%
% Version 1, December 2015 First Release

% Definition of fields that will be checked
checkfield={'Abs' 'X' 'filelist' 'Ex' 'Em' 'nEm' 'nEx' 'IntensityUnit' 'nSample'};
nChecks=size(checkfield,2);
result=zeros(nChecks,1);
for n=1:nChecks
    result(n,1)=isfield(DS,checkfield{1,n});
end

if all(result)
    result='PASSED';
else
    missing=checkfield(1,find(not(result)));
    result='NOT PASSED';
        warning('Missing field(s):')
        disp(char(missing));
end
end