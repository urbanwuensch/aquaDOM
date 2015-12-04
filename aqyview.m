function aqyview(DS,selection,view)
% View calculated quantum yields
%   
% USEAGE:
%           aqyview(DS,selection,view)
%
% INPUTS
%               DS: Structure containing quantum yields and quantum yield precision estimates
%        selection: selection of samples to be viewed (e.g. [1 2])
%             view: viewing preference, 'SUP' for superposition of AQYs or 'SEQ' for sequential display
%
% OUTPUTS
%               plot
%             
% Examples:
%       aqyview(qy_X,[2 17],'SUP'); superposition view, sampels 2 and 17
%       aqyview(qy_X,[],'SUP'); sequential view of all quantum yields in qy_X
%
% Notice:
% This mfile is part of the aquaDOM toolbox.
%
%
% aqyview.m: Copyright (C) 2015 Urban J Wünsch
% Technical University of Demark
% National Institute of Aquatic Resources
% Section for Marine Ecology and Oceanography
% Kavalergården 6
% 2920 Charlottenlund, Denmark
% urbw@aqua.dtu.dk
%
% Version 1, June 2015 First Release

% Welcome message
disp(' ')
disp('aqyview.m');close all;
disp(' ')


% Initial subsetting of data
if isempty(selection)
    qy=DS.qy();
    qyE=DS.qy_precision;
    names=DS.name();
else
    qy=DS.qy(selection,:);
    qyE=DS.qy_precision(selection,:);
    names=DS.name(:,selection);
end
if isfield(DS,'ref_name')
    disp(['AQY standard: ',char(DS.ref_name(1,1))]);
end
try disp(['AQY calculation method: ',char(DS.qy_method(1,:))]);
catch
end
if size(names,2)>48
    error('A superimposed plot with more than 48 samples is not supported. Try a smalller selection instead, or use the sequential plotting option')
end
% Plotting
if view=='SUP'
    %% Superimposed plotting with differnet line and marker styles
    linestyles = cellstr(char('-',':','-.','--',...
            '-',':','-.','--',...
            '-',':','-.','--',...
            '-',':','-.','--',...
            '-',':','-.','--',...
            '-',':','-.','--',...
            '-',':','-.','--',...
            '-',':','-.','--',...
            '-',':','-.','--',...
            '-',':','-.','--',...
            '-',':','-.','--',...
            '-',':','-.','--'));
    Markers=['o','x','+','*','s','d','v','^','<','>','p','h','.',...
            '+','*','o','x','^','<','h','.','>','p','s','d','v',...
            'o','x','+','*','s','d','v','^','<','>','p','h','.',...
            'o','x','+','*','s','d','v','^','<','>','p','h','.',...
            'o','x','+','*','s','d','v','^','<','>','p','h','.',...
            '+','*','o','x','+','*','s','d','v','^','<','>','p','h','.',...
            '+','*','+','*',];
    %MarkerEdgeColors=jet(DS.nSample);
    MarkerEdge=['r', 'g', 'b','c','m','k',...
            'k' ,'m', 'c' ,'b' ,'g','r',...
            'r', 'g', 'b','c','m','k',...
            'k' ,'m', 'c' ,'b' ,'g','r',...
            'r', 'g', 'b','c','m','k',...
            'k' ,'m', 'c' ,'b' ,'g','r',...
            'r', 'g', 'b','c','m','k',...
            'k' ,'m', 'c' ,'b' ,'g','r'];
    figure('InvertHardcopy','off','Color',[1 1 1],'units','normalized','outerposition',[0.25 0.4 0.5 0.5]); set(gcf,'numbertitle','off','name','aqyview: AQY spectra');
    hold on;grid on;box on;
    xlabel('Excitation wavelength')
    ylabel('AQY')
    set(gca,'fontsize',10,'FontName','times new roman');

    xlim([min(DS.Ex) max(DS.Ex)])
    for i=1:size(qy,1)
        h=errorbar(DS.Ex,(qy(i,:)),qyE(i,:));h.LineStyle='none'; h.LineWidth=1; h.Color=[0.7 0.7 0.7];
        p(i)=plot(DS.Ex,(qy(i,:)),[linestyles{i}...
            %,Markers(i)...
            ],'Color',MarkerEdge(i));%p(i).Color='k';
        p(i).LineWidth=1;
    end
    legend(p,names)
elseif view=='SEQ' % Sequential plotting of AQY spectra with jump-back function
    %% Sequential plotting
    figure('units','normalized','outerposition',[0.25 0.4 0.5 0.5]); set(gcf,'numbertitle','off','name','aqyview: AQY spectra');
    for i=1:size(qy,1)
            if i==1;
                plt=1;
                f=1;
            end
    
    try
        figure(1)
        grid on;
        xlabel('Excitation wavelength')
        ylabel('AQY')
        set(gca,'fontsize',10,'FontName','times new roman');
        h=errorbar(DS.Ex,qy(i,:),qyE(i,:));h.LineStyle='none'; h.LineWidth=1; h.Color=[0.7 0.7 0.7];
        hold on;
        p=plot(DS.Ex,qy(i,:)); p.Color='k'; p.LineWidth=0.5;
        legend(names(i))
        hold off;
        xlim([min(DS.Ex) max(DS.Ex)])
        ylim([0 nanmax(nanmax(qy))])
        set(gca,'fontsize',10,'FontName','times new roman');

    catch
        disp(' ')
        disp(['Could not plot the selected sample (',plt,') try again'])
        disp(' ')
        plt=plt+f;
    end
    
    
        try f=input('Press Enter to proceed to the next sample, or a number to go back n samples: ');
        catch disp('Please only use numbers or press enter. Try again ;) ');
        end
        plt=plt+1;
        if ~isempty(f);
            plt=(plt-1)-f;
        end
    end
end
end