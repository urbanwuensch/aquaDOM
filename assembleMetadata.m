function [ DS ] = assembleMetadata(tags,Em_range,DStype)
% assemble metadata to process EEMs by aquaDOM
%   
% USEAGE:
%           [DS]=assembleMetadata(tags,Em_range,DStype)
%
% INPUTS
%               tags: Cell array that contains data tags, separated by space
%           Em_range: matrix containing the range for emission signal integration
%             DStype: differentiator, 'STD' for AQY standards, 'SAM' for AQY samples
%             NOTE: For DStype 'STD', the reference QY, the reference QY
%                       wavelength, the number of measurement per standard and the
%                       concentration of the standard will be required during the execution
%                   For DStype 'SAM', the number of measurement per sample and the
%                       concentration of the sample will be required during the execution
%
% OUTPUTS
%                 DS: structure containing assembled metadata
%             
% Examples:
%       metadataStandards=assembleMetadata({QS SA},[325 600],'STD');  Metadata assembly for 2 AQY standards with tags QS and SA, emission integration from 325 to 600 nm
%       metadataSamples=assembleMetadata({'DANA15.2.15' 'DANA15.2.420' 'DANA15.4.2'},[240 600],'SAM'); Metadata assembly for 3 samples with the full Aqualog emission range
%
% Notice:
% This mfile is part of the aquaDOM toolbox.
%
%
% assembleMetadata.m: Copyright (C) 2015 Urban J Wünsch
% Technical University of Demark
% National Institute of Aquatic Resources
% Section for Marine Ecology and Oceanography
% Kavalergården 6
% 2920 Charlottenlund, Denmark
% urbw@aqua.dtu.dk
%
% Version 1, June 2015 First Release

DS.conc=nan(size(tags,2),10);

disp('  ');
disp('assembleMetadata.m');
% Just to let user know about the function input in DStype.
if DStype=='STD'
   disp('Metadata for AQY standards:')
end
errorflag=0;
% Input in tags is processed one by one
for n=1:size(tags,2)
    % Two extra metadata variables are created in case of AQY standards
    if DStype=='STD'
        try
        prompt = ['Enter the reference quantum yield for: ', char(tags(1,n)),' (e.g 0.51):  '];
        DS.ref(1,n) = input(prompt);
        catch
            error('The reference quantum yield has to be a number, e.g. 0.51')
        end
        try
        prompt = ['Enter the reference qy wavelength in nm for: ', char(tags(1,n)),':   '];
        DS.ref_wave(1,n) = input(prompt);
        catch
            error('The reference wavelength has to be a number, e.g. 250.')
        end
    end
    % assignment data tags and names
    DS.tag(1,n)=tags(1,n);
    DS.name(1,n)=DS.tag(1,n);
    DS.Em_range(n,:)=[Em_range(1);Em_range(2)];
    
    % Details for n of measurements for sample i
    try
    prompt = ['Number of measurements for: ', char(DS.name(1,n)),':   '];
    x = input(prompt);
    catch
        error('You have to provide a number here.')
    end
    
    % for case n measurement > 1: Details on concentration (if desired)
    if x>1
        try
        disp('Provide molar concentration if desired, otherwise, press enter (dummy values will be assigned)')
        prompt = '(e.g. [c1 c2 c3]):   ';
        x1 = input(prompt);
        catch
            error('You have to provide a number here. FORMAT: [conc1 conc2 conc3 ...], eg. [6.12E-06 5.10E-06 4.25E-06]')
        end
        if ~isempty(x1)
            DS.conc(n,1:size(x1,2))=transpose(x1);
        else 
            DS.conc(n,1:x)=1:x;
        end
    elseif x==1
        DS.conc(n,1)=1;
    end
    % # of conc for sample n
    DS.nConc(n,1)=size(find(~isnan(DS.conc(n,:))),2);
    if x~=DS.nConc(n,1)
        warning('Missmatch between number of concentrations and provided concentrations. Try again');
        n=n-1;
        errorflag=1;
    end
end
if errorflag~=1
    disp('Done. Metadata variable created. You can now run splitDS.m')
elseif errorflag==1
    warning('The script finished, but problems were encountered.')
end