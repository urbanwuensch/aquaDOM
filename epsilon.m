function [DS]=epsilon(DS)
% Calculate molar fluorescence and absorbance 
%   
% USEAGE:
%           [DS]=epsilon(DS)
%
% INPUTS
%               DS: Structure containing variables .X, .Abs in 4 and 3D data matrices (converted from drEEM format by splitDS.m)
%
% OUTPUTS
%               DS: Structure with varable Xmf and epsilon for molar fluorescence (in R.U. L mol-1)and absorbance (in "L mol-1 cm-1")
%             
% Examples:
%       [qys]=epsilon(qys);
%
% Notice:
% This mfile is part of the aquaDOM toolbox.
%
%
% epsilon.m: Copyright (C) 2015 Urban J Wünsch
% Technical University of Demark
% National Institute of Aquatic Resources
% Section for Marine Ecology and Oceanography
% Kavalergården 6
% 2920 Charlottenlund, Denmark
% urbw@aqua.dtu.dk
%
% Version 1, August 2015 First Release
DS.epsilon=zeros(DS.nSample,DS.nEx);

%% Molar absorbance
for i=1:DS.nSample
    for m=1:DS.nEx
    p= polyfit(sort(DS.conc(i,~isnan(DS.conc(i,:)))),sort(DS.Abs(i,~isnan(DS.Abs(i,:,m)),m)),1);
    DS.epsilon(i,m)=p(1);
    end
end

%% Molar fluorescence
% Counters for progress
total=DS.nSample*DS.nEx*DS.nEm;
current=1;
% EEMs contain large numbers of data points, hence linear fits can take a
% while. The user is informed with a waitbar
h = waitbar(0,['Please wait! Processing ',char(num2str(total)),' linear fits']);
for i=1:DS.nSample
    for m=1:DS.nEx
        for j=1:DS.nEm
            % Raman normalization for the molar fluorescence
            Xrn=squeeze(DS.X(i,~isnan(DS.conc(i,:)),j,m))/DS.RamanArea(i*DS.nConc(i),1);
            p=polyfit(DS.conc(i,~isnan(DS.conc(i,:))),Xrn,1);
            DS.Xmf(i,j,m)=p(1);
            % Update waitbar
            waitbar(current/total)
            current=current+1;
        end
    end
end
close(h)
figure('InvertHardcopy','off','Color',[1 1 1]);
set(gcf, 'units', 'centimeters', 'pos', [5 5 18 10])
subplot(1,2,1)
for n=1:DS.nSample
    plot(DS.Abs_wave,DS.epsilon(n,:),'LineWidth',1);hold on;
    title(DS.name(1,n));
    xlim([240 600]);
end
xlabel('Absorbance wavelength');
ylabel('\epsilon [L mol^-^1 cm^-^1]');
set(gca,'fontsize',8,'FontName','times new roman');
title('Molar absorbance')
ylim([0 inf])
legend(DS.name)


subplot(1,2,2)
for n=1:DS.nSample
    [~,idx]=nanmax(nanmax(DS.Xmf(n,:,:)))
    plot(DS.Em,DS.Xmf(n,:,idx),'LineWidth',1);hold on;
    title(DS.name(1,n));
    xlim([240 600]);
end
xlabel('Emission wavelength');
ylabel('Molar fluorescence (R.U. mol-1 L)');
set(gca,'fontsize',10,'FontName','times new roman');
title('Molar fluorescence')
ylim([0 inf])
legend(DS.name)
end

