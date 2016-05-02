function [DS]=aqy(DS,mode,method,qy_Std,qy_std_choice)
% Calculate quantum yields
%   
% USEAGE:
%           [DS]=aqy(DS,mode,method,qy_Std,qy_std_choice)
%
% INPUTS
%               DS: Structure containing variables .X, .Abs, and processed by splitDS.m
%             mode: usage mode, 1: cross-calibration, 2: AQY calculation of samples
%           method: calculation method, 1: variable intercept, 2: zero-intercept
%           qy_Std: Structure containing the information for AQY standards after cross calibtration (leave empty for method 1)
%    qy_std_choice: Choice of standard to be uesd for AQY calculations (leave empty for method 1)
%
% OUTPUTS
%               DS: original structure, supplemented with quantum yields in variable .qy, and quantum yield precision in .qy_precision
%             
% Examples:
%       [qy_R]=aqy(qy_R,1,1); Calibration, variable intercept approach
%       [qy_X]=aqy(qy_X,2,1,qy_R,1); Calculation, variable intercept approach, QY Std.
%
% Notice:
% This mfile is part of the aquaDOM toolbox.
%
%
% aqy.m: Copyright (C) 2015 Urban J Wünsch
% Technical University of Demark
% National Institute of Aquatic Resources
% Section for Marine Ecology and Oceanography
% Kavalergården 6
% 2920 Charlottenlund, Denmark
% urbw@aqua.dtu.dk
%
% Version 1, June 2015 First Release

% Welcome messages
disp(' ');close all;
disp('SETTINGS:')
if mode==1
    disp('Cross-calibration mode')
elseif mode==2
    disp('Sample mode')
end

if method==1
    disp('Variable intercept method selected')
    disp(' ')
    disp(' ')
    DS.qy_method='Variable intercept method';
elseif method==2
    disp('F''/A zero intercept method selected')
    disp(' ')
    disp(' ')
    DS.qy_method='zero intercept method';
end

%% Initialization of new variables / functions / thresholds
% These tresholds make sure that the calculated quantum yields are of satisfying quality
DS.Abs_s_n_treshold=25;
Rcutoff=0.95;
% Plotting function for confidence interval




%% Allocation of variables
DS.qy=NaN(DS.nSample,DS.nEx);
AUC=zeros(DS.nSample,size(DS.conc,2),DS.nEx);
if method==1
    slope=NaN(DS.nSample,DS.nEx);
    DS.qy_precision=NaN(DS.nSample,DS.nEx);
    slope_95=NaN(DS.nSample,DS.nEx,2);
end
if method==2
    DS.F_a=NaN(DS.nSample,DS.nEx);
    DS.qy_precision=NaN(DS.nSample,DS.nEx);
    DS.Abs_sqE=NaN(DS.nSample,DS.nEx);
    DS.AUC_sqE=NaN(DS.nSample,DS.nEx);
end

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

% F'/A slope calculation
% Robust fit with intercept ~=0
if method==1
    warning off;
    
    disp('F''/Abs slope calculation...')
    for n=1:DS.nSample
        i=1;
        for m=1:DS.nEx
            % Robust fit of F'/A
            if DS.nConc(n)>3
                mdlr = fitlm(DS.Abs(n,:,m).',AUC(n,:,m).','RobustOpts','on');
            else
                mdlr = fitlm(DS.Abs(n,:,m).',AUC(n,:,m).');
            end
%             scatter(DS.conc(n,:),AUC(n,:,m));pause;
%             %95% confidence intervall of slope parameter (alpha=0.05)
%             plot(mdlr); pause();close all;
%             disp(mdlr.Rsquared.Adjusted);disp((DS.Abs(n,DS.nConc(n),m)/DS.noise(n,1)));
             coef=coefCI(mdlr,0.05);
             slope_95(n,m,:)=coef(2,:);
            % Retrieval of slope parameter
            if mdlr.Rsquared.Adjusted>Rcutoff&&(DS.Abs(n,DS.nConc(n),m)/DS.noise(n,1))>(DS.Abs_s_n_treshold*DS.noise(n,1));
                slope(n,m)=mdlr.Coefficients{2,1};
                i=1+1;
            end
            clear mdlr coeff
        end
        
        if i==1
            warning on;
            warning('Quality standards during slope calculations were never met. Something might be wrong with your data. Check plots for diagnosis')
            subplot(1,2,1)
            scatter(DS.conc(n,:),squeeze(DS.X(n,:,1,1)),'filled');title('Fluorescence vs. conc');hold on;
            plot(DS.conc(n,:),squeeze(DS.X(n,:,1,1)))
            subplot(1,2,2)
            scatter(DS.conc(n,:),squeeze(DS.Abs(n,:,5)),'filled');title('Absorbance vs. conc');hold on;
            plot(DS.conc(n,:),squeeze(DS.Abs(n,:,5)))
            error(['No quantum yields calculated for ',DS.filelist{1,:}])
        end
    end
end


% Linear fit with intercept =0 and n=1
if method==2
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
end

%% Quantum yield calculations
%% Cross-calibration using two qy standards
if mode==1
    disp('Cross calibration in progress...')
    disp(' ')
    disp(' ')
    disp('RESULTS:')

    % Determination of the reference excitation position (or closest match)
    for n=1:DS.nSample
        DS.posRef_wave(n)=knnsearch(DS.Ex,DS.ref_wave(1,n));
    end
    
    % For the cross calibration, flips will flip some variables
    posRef_wave_inv = fliplr(DS.posRef_wave);
    ref_inv = fliplr(DS.ref);
    
    if method==1
        slope_inv = flipud(slope);
    end
    
    if method==2
        F_a_inv = flipud(DS.F_a);
        AUC_sqE_inv= flipud(DS.AUC_sqE);
        Abs_sqE_inv= flipud(DS.Abs_sqE);
    end
    
    % AQY standards
    for n=1:DS.nSample
        for m=1:DS.nEx
            if min(DS.Abs(n,:,m))>=DS.Abs_s_n_treshold*DS.noise(n)
                if method==1
                    % Creation of temporary variables for better overview
                    a=ref_inv(1,n);
                    b=slope(n,m);
                    c=slope_inv(n,posRef_wave_inv(1,n));
                    
                    a1=ref_inv(1,n);
                    b1=slope_95(n,m,1);
                    c1=slope_inv(n,posRef_wave_inv(1,n));
                                      
                    % Quantum yield estimate
                    DS.qy(n,m)=qy(a,b,c);
                    % Precision estimate based on 95% confidence interval +/- value
                    DS.qy_precision(n,m,1)=qy(a,b,c)-qy(a1,b1,c1);
                    clearvars a b c a1 b1 c1 a2 b2 c2 
                    DS.slope=slope;
                end
                if method==2
                    
                    a=ref_inv(1,n);
                    b=DS.F_a(n,m);
                    c=F_a_inv(n,posRef_wave_inv(1,n));
                    
                    
                    DS.qy(n,m)=qy(a,b,c);
                    clearvars a b c
                    
                    % Uncertainty of QY based on standard error propagation, fist: definition of terms (to avoid confusion)
                    const=ref_inv(1,n);
                    a=AUC_sqE_inv(n,m);
                    b=DS.AUC_sqE(n,m);
                    c=Abs_sqE_inv(n,m);
                    d=DS.Abs_sqE(n,m);
                    
                    % Calculation of uncertainty
                    relErr=const * sqrt( a + b + c + d );
                    DS.qy_precision(n,m)= relErr*DS.qy(n,m);
                end
            end
        end
        DS.error(n) = abs(100-100/DS.ref(1,n)*DS.qy(n,DS.posRef_wave(1,n)));
    end
    %% Evaluation of results and plotting
    figure; set(gcf,'numbertitle','off','name','Results: Quantum yield cross calibration');
    % Plot for F' & A vs. measurement 
    warning off;
    for n=1:DS.nSample
        subplot(3,DS.nSample,n) ;plotAdded(fitlm(DS.conc(n,:),squeeze(DS.Abs(n,:,DS.posRef_wave(1,n))))); 
        legend('off');title(' ');
        if method==1
            title(['Dilution series for ', char(num2str(DS.name{1,n}))]);
        elseif method==2
            title(['Data variablility for ', char(num2str(DS.name{1,n}))]);
        end
        ylabel(['Absorbance ',char(num2str(DS.ref_wave(1,n))),' nm']);
        xlabel('Concentration');
    end
    for n=1:DS.nSample
        subplot(3,DS.nSample,2+n); plotAdded(fitlm((DS.conc(n,:)),(squeeze(AUC(n,:,DS.posRef_wave(1,n))))));
        legend('off'); title(' ');
        ylabel(['F'' ',char(num2str(DS.ref_wave(1,n))),' nm']); xlabel('Concentration')
    end
    warning on;

    for n=1:DS.nSample
        subplot(3,DS.nSample,4+n)
        if method==1
           qy_plot_errBar(DS.Ex,DS.qy(n,:),DS.qy_precision(n,:,:),mode,n,DS)
        end
        if method==2
           qy_plot_errBar(DS.Ex,DS.qy(n,:),DS.qy_precision(n,:),mode,n,DS)
        end
    end
    
    if max(DS.error(1,:))<10
        disp('Cross-calibration PASSED. (Deviation from lit. values <10%)')
        disp('___________________')
    else
        warning('Cross-calibration FAILED. (Deviation from lit. values <10%)')
        disp('___________________')
    end
end


%% Calculation of sample quantum yields
if mode==2
    % AQY calculation for samples
    disp('Sample quantum yield calculation in progress...')
    for n=1:DS.nSample
        for m=1:DS.nEx
            if min(DS.Abs(n,:,m))>=DS.Abs_s_n_treshold*DS.noise(n)
                if method==1
                    a=qy_Std.ref(qy_std_choice);
                    b=slope(n,m);
                    c=qy_Std.slope(qy_std_choice,qy_Std.posRef_wave(1,qy_std_choice));
                    
                    a1=qy_Std.ref(qy_std_choice);
                    b1=slope_95(n,m,1);
                    c1=qy_Std.slope(qy_std_choice,qy_Std.posRef_wave(1,qy_std_choice));
                                       
                    DS.qy(n,m)=qy(a,b,c);
                    DS.qy_precision(n,m)=qy(a,b,c)-qy(a1,b1,c1);
                    clearvars a b c a1 b1 c1 a2 b2 c2
                end
                if method==2
                    
                    a=qy_Std.ref(qy_std_choice);
                    b=DS.F_a(n,m);
                    c=qy_Std.F_a(qy_std_choice,qy_Std.posRef_wave(1,qy_std_choice));
                    
                    DS.qy(n,m)=qy(a,b,c);
                    
                    clearvars a b c
                    
                    % Just to make it less confusing, the terms for SE will temp. be stored in variables
                    const=qy_Std.ref(qy_std_choice);
                    a=DS.AUC_sqE(n,m);
                    b=qy_Std.AUC_sqE(qy_std_choice,qy_Std.posRef_wave(1,qy_std_choice));
                    c=DS.Abs_sqE(n,m);
                    d=qy_Std.Abs_sqE(qy_std_choice,qy_Std.posRef_wave(1,qy_std_choice));
                    
                    % Actual SE calculation
                    relErr= const * sqrt( a + b + c + d );
                    DS.qy_precision(n,m)= relErr*DS.qy(n,m);
                    
                    %DS.errors(n,m,:)=[a b c d relErr];
                end
            end
        end
    end
    
    % plotting
    figure
    for n=1:DS.nSample
        subplot(round(DS.nSample/2),2,n)
        if method==1
            qy_plot_errBar(DS.Ex,DS.qy(n,:),DS.qy_precision(n,:),mode,n,DS)
            DS.ref_name=qy_Std.name(1,qy_std_choice);
        end
        if method==2

            qy_plot_errBar(DS.Ex,DS.qy(n,:),DS.qy_precision(n,:),mode,n,DS);
            DS.ref_name=qy_Std.name(1,qy_std_choice);
        end
    end
    disp(' ')
    disp(' ')
    disp('Quantum yield calculation complete.')
    disp('___________________________________') 
end
disp(' ')
toc
end