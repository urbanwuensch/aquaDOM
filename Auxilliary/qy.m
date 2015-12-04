function [aqy] = qy(a,b,c)
% %Calculation of quantum yields
% USEAGE:
%           [aqy] = qy(a,b,c)
%
% INPUTS
%               a:  reference quantum yield
%               b:  sample F/A or slope
%               c:  reference F/A or slope
%
% OUTPUTS
%             aqy: resulting quantum yield
%             
%
% Notice:
% This mfile is part of the aquaDOM toolbox.
%Calculation of quantum yield

% (c) Urban Wuensch, DTU AQUA 06/2015

% qy.m: Copyright (C) 2015 Urban J Wünsch
% Technical University of Demark
% National Institute of Aquatic Resources
% Section for Marine Ecology and Oceanography
% Kavalergården 6
% 2920 Charlottenlund, Denmark
% urbw@aqua.dtu.dk
%
% Version 1, June 2015 First Release

aqy=a*(b/c);
end