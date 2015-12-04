function qy_plot_errBar(x,y,qy_Precision,mode,n,DS)
% Function for plotting quantum yields with error bars 
% USEAGE:
%           [DS]=mergeDS(varargin)
%
% INPUTS
%               x: Excitiation wavelength
%               y: quantum yield estimate
%    qy_Precision: precision estimate
%            mode: Information needed for distinguishing between standard or sample mode for plotting
%               n: counter, function usually is called in nSample loop multiple times
%              DS: Original DS in qycalc, needed for some additional info
%
% OUTPUTS
%            plot
%             
% Examples:
%       DS=mergeDS(qy_X,qy_X1,qy_X3); merging of 3 datasets into a single structure
%
% Notice:
% This mfile is part of the aquaDOM toolbox.
%
%
% qy_plot_errBar.m: Copyright (C) 2015 Urban J Wünsch
% Technical University of Demark
% National Institute of Aquatic Resources
% Section for Marine Ecology and Oceanography
% Kavalergården 6
% 2920 Charlottenlund, Denmark
% urbw@aqua.dtu.dk

%% Data formatting
y_err=qy_Precision;
x=x';
%% Plotting
h=errorbar(x,y,y_err);h.LineStyle='none'; h.LineWidth=1; h.Color=[0.8 0.8 0.8];
hold on;
p=plot(x,y); p.Color='k'; p.LineWidth=0.5;

% Plot formatting
if mode==1
    scatter(x(1,knnsearch(DS.Ex,DS.ref_wave(1,n))),y(1,knnsearch(DS.Ex,DS.ref_wave(1,n))),'filled');
    tnames=fliplr(DS.name);
    disp(['Deviation when using ',char(tnames(1,n)),' as QY standard: ',num2str(DS.error(1,n)),'%']);
end

title(sprintf('%s', DS.name{1,n})); xlabel('Excitation wavelength'); ylabel('\Phi');
xlim([min(DS.Ex) max(DS.Ex)]);





end