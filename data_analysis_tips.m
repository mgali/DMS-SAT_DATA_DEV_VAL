% SOME TWEAKS FOR THE DMS_P SATELLITE MATCH-UP DATABASE
% 
% Marti Gali Tapias, 2018-11-30
% marti.gali.tapias@gmail.com
% 
% Check these references (especially the supporting online materials) for
% descriptions of quality control procedures and filtering criteria
% 
% Galí, M., Devred, E., Levasseur, M., Royer, S. J., & Babin, M. (2015)
% A remote sensing algorithm for planktonic dimethylsulfoniopropionate
% (DMSP) and an analysis of global patterns. Remote Sensing of Environment,
% 171, 171-184.
% https://doi.org/10.1016/j.rse.2015.10.012
% 
% Galí, M., Levasseur, M., Devred, E., Simó, R., & Babin, M. (2018)
% Sea-surface dimethylsulfide (DMS) concentration from satellite data
% at global and regional scales. Biogeosciences, 15(11), 3497-3519.
% https://doi.org/10.5194/bg-15-3497-2018
% 
% NOTE1: I can add the code for quantile regression on request, as well as
% any other code used to produce the curated database. But this may take
% some time.
% 
% NOTE2: yes, I could have created a function to do the stepwise merging
% of satellite matchups obtained at different time-space resolutions to
% save many lines of code...

%% Load dataset
% Uncomment below to use *.mat data (faster!)
load ./standard_QCed_DMS_P_db_11-Dec-2018.mat
openvar('head')     % Display header variable

% % Uncomment to save as comma separated values
% csvwrite('./standard_QCed_DMS_P_db_11-Dec-2018.txt',DATA);

%% Tips

% % Uncomment to assign vector of values to each of the header variables
% for j=1:length(head), eval([head{j,1} ' = DATA(:,j);']); end

% % Example of how to select a column using its variable name (creating a
% % logical vector) to avoid the use of numeric indices
% myvar = 'dms';
% mydata = DATA(:,strcmp(head,myvar));

%% Define conversion factors, thresholds, constants, to be applied later on
% ----------------------------- EDIT TO PLAY ------------------------------
convfact.hplc2turner = 1.39; % [relative units], based on Sathyendranath et al. 2009 MEPS
threshold.sal.min = 30; % [PSU]
threshold.zbot.max = -200; % [m], note more negative is deeper
threshold.dms.min = 0.1; % [nM]
threshold.dms.max = 100; % [nM]
threshold.dmspt.min = 1; % [nM]
threshold.minfract = 0.1; % [relative units]
quota.micro = 0.01; % [mol DMSPcarbon / mol cell carbon]
quota.nano = 0.055; % [mol DMSPcarbon / mol cell carbon]
quota.pico = 0.055; % [mol DMSPcarbon / mol cell carbon]

%% Merge SeaWiFS and MODIS-A daily matchups
% NOTE: try to analyze 1 pixel, 3x3 and 5x5 pixels matchups separately, for each sensor
chl_day = nanmean(DATA(:,[87 88]),2);                               % start with 1 pixel matchups
fprintf('Daily 1 pix matchups N = %i\n',sum(~isnan(chl_day)))       % Check N: should give 3324
chl_day(isnan(chl_day)) = nanmean(DATA(isnan(chl_day),[89 90]),2);	% fill gaps with 3x3 pixel matchups
fprintf('Daily 3x3 pix matchups N = %i\n',sum(~isnan(chl_day)))     % Re-check N: should give 4369
chl_day(isnan(chl_day)) = nanmean(DATA(isnan(chl_day),[91 92]),2);	% fill gaps with 5x5 pixel matchups
fprintf('Daily 5x5 pix matchups N = %i\n',sum(~isnan(chl_day)))     % Re-check N: should give 5183

%% Merge SeaWiFS and MODIS-A 8-day matchups
% NOTE: as above
chl_8d = nanmean(DATA(:,[51 56]),2);                                % start with 1 pixel matchups
fprintf('8-day 1 pix matchups N = %i\n',sum(~isnan(chl_8d)))        % Check N: should give 14327
chl_8d(isnan(chl_8d)) = nanmean(DATA(isnan(chl_8d),[52 57]),2);     % fill gaps with 3x3 pixel matchups
fprintf('8-day 3x3 pix matchups N = %i\n',sum(~isnan(chl_8d)))      % Re-check N: should give 14327
chl_8d(isnan(chl_8d)) = nanmean(DATA(isnan(chl_8d),[53 58]),2);     % fill gaps with 5x5 pixel matchups
fprintf('8-day 5x5 pix matchups N = %i\n',sum(~isnan(chl_8d)))      % Re-check N: should give 14338

%% Check overall relationship between in situ and satellite chlorophyll

% Day matchups, display fit (fprintf)
y = real(log10(chl_day));
X=[ones(size(DATA(:,9))) real(log10(DATA(:,9)))];
[b,bint,r,rint,stats]=regress(y,X);
fprintf('\nIntercept_day = %0.2f (%0.2f__%0.2f),\nSlope_day = %0.2f (%0.2f__%0.2f),\nR^2 = %0.2f,\np = %0.1e\n',...
    b(1), bint(1,1), bint(1,2),...
    b(2), bint(2,1), bint(2,2),...
    stats(1), stats(3))
regr.oneD = struct('b',b,'bint',bint,'r',r,'rint',rint,'stats',stats);   % save in structure   

% 8-day matchups, display fit (fprintf)
y = real(log10(chl_8d));
X=[ones(size(DATA(:,9))) real(log10(DATA(:,9)))];
[b,bint,r,rint,stats]=regress(y,X);
fprintf('\nIntercept_8d = %0.2f (%0.2f__%0.2f),\nSlope_8d = %0.2f (%0.2f__%0.2f),\nR^2 = %0.2f,\np = %0.1e\n',...
    b(1), bint(1,1), bint(1,2),...
    b(2), bint(2,1), bint(2,2),...
    stats(1), stats(3))
regr.eightD = struct('b',b,'bint',bint,'r',r,'rint',rint,'stats',stats);   % save in structure

% Plot
figure()
loglog(DATA(:,9),chl_8d,'.g'), hold on
loglog(DATA(:,9),chl_day,'.b')
plot([0.01 30],[0.01 30],'-k')
axis([0.009 30 0.009 30])
xlabel('Chl in situ'); ylabel('Chl satellite')
legend('8-day matchups','day matchups','location','northwest')

%% Identify DMS/DMSPt/Chl in situ according to flags
% Note: Chl < 0.02 mg m-3 was removed beforehand (flagChl, column 14).
% This removed mainly chl = 0.01 mg m-3 from the Atlantic Meridional Transect program.
ilist.flagestuarine = DATA(:,13)==0;             % estuarine sample
ilist.flagAcidPhaeo = DATA(:,15)==0;             % acidified Phaeocystis-containing sample
ilist.flagSratio = DATA(:,16)==0;                % oustide 0.01-0.99 loglog linear quantile regression DMS vs. DMSPt
ilist.flagS2Chl = DATA(:,17)==0;                 % oustide 0.01-0.99 loglog linear quantile regression DMSPt vs. Chl
ilist.allflags = sum(DATA(:,[13 15 16 17]),2)<4; % all quality control flags at once

%% Identify DMS/DMSPt/Chl in situ according to salinity and bottom depth
% Using sal > 30 and zbot > 200 filters out several continental shelf data
% and Arctic waters with strong riverine and/or ice-melt influence
% NOTE: only ~80% of chl-dmspt data pairs have salinity data
% Apply to in situ dms, dmspt and chl only
ilist.sal.min = ~isnan(DATA(:,11)) & DATA(:,11) < threshold.sal.min;    % sal
ilist.zbot.max = ~isnan(DATA(:,37)) & DATA(:,37) > threshold.zbot.max;  % zbot

%% Identify DMSP measurements below minimal expected DMSP
% Calculate expected DMSP based on phytoplankton size classes, Chl to 
% carbon conversion (same for all size classes) and DMSP:C molar quotas
% Based on combination of Uitz et al. 2006 (JGR) size classes and
% Sathyendranath et al. 2009 (MEPS) for Chl to carbon conversion
% NOTE: Cmicro+Cnano+Cpico=CHL_in_situ (prior to HPLC or satellite
% adjustments made below. Tha is:
% nansum(DATA(:,43:45),2) is approximately equal to DATA(:,9)
microdmsp = 10.^(1.81 + 0.63*log10(DATA(:,43)))*1000/12*quota.micro*(1/5);
nanodmsp = 10.^(1.81 + 0.63*log10(DATA(:,44)))*1000/12*quota.nano*(1/5);
picodmsp = 10.^(1.81 + 0.63*log10(DATA(:,45)))*1000/12*quota.pico*(1/5);
% Multiply expected DMSP by 10% (or selected %) for conservative exclusion
minsumdmsp = threshold.minfract*(picodmsp + nanodmsp + microdmsp); 
ilist.minsumdmsp = DATA(:,8)<minsumdmsp;

%% Identify too low and too high DMS and too low DMSPt
ilist.dmsext = DATA(:,7)<threshold.dms.min | DATA(:,7)>threshold.dms.max;
ilist.dmspt.min = DATA(:,8)<threshold.dmspt.min;

%% Identify coccolithophore bloom datasets. Force into the "stratified" subset
% Identified by contribnum
% Constrain to chl concentrations between 0.35 and 4 (range based on
% inspection of the cocco-dominated datasets)
ilist.cocco = (DATA(:,1)==41 | DATA(:,1)==77 | DATA(:,1)==133 | DATA(:,1)==135 | DATA(:,1)==176 | DATA(:,1)==181)...
    & DATA(:,9)<4 & DATA(:,9)>0.35;
DATA(ilist.cocco,41) = 1;

%% Plot to check data flagged or deemed non-representative of pelagic ecosystem
figure(), set(gcf,'units','centimeters','position',[1 1 33 23])
loglog(DATA(:,8),DATA(:,7),'.k','markersize',10); hold on
loglog(DATA(ilist.flagestuarine,8),DATA(ilist.flagestuarine,7),'.r','markersize',10)
loglog(DATA(ilist.sal.min,8),DATA(ilist.sal.min,7),'sr','markersize',5)
loglog(DATA(ilist.zbot.max,8),DATA(ilist.zbot.max,7),'.y','markersize',5)
loglog(DATA(ilist.flagAcidPhaeo,8),DATA(ilist.flagAcidPhaeo,7),'+m','markersize',5)
loglog(DATA(ilist.flagSratio,8),DATA(ilist.flagSratio,7),'xb','markersize',8)
loglog(DATA(ilist.dmsext,8),DATA(ilist.dmsext,7),'sb','markersize',8)
loglog(DATA(ilist.dmspt.min,8),DATA(ilist.dmspt.min,7),'+c','markersize',5)
loglog(DATA(ilist.minsumdmsp,8),DATA(ilist.minsumdmsp,7),'og','markersize',8)
axis([0.1 3000 0.01 300])
legend('all data','estuarine','low salinity','shallow bottom',...
    'acidified Phaeo','outside 1-99 quantile regr',...
    'extreme DMS','too low DMSPt','too low DMSP based on Chl',...
    'location','southeast')
xlabel('DMSPt in situ (nM)'); ylabel('DMS in situ (nM)')

figure(), set(gcf,'units','centimeters','position',[10 3 33 23])
loglog(DATA(:,9),DATA(:,8),'.k','markersize',10); hold on
loglog(DATA(ilist.flagestuarine,9),DATA(ilist.flagestuarine,8),'.r','markersize',10)
loglog(DATA(ilist.sal.min,9),DATA(ilist.sal.min,8),'sr','markersize',5)
loglog(DATA(ilist.zbot.max,9),DATA(ilist.zbot.max,8),'.y','markersize',5)
loglog(DATA(ilist.flagAcidPhaeo,9),DATA(ilist.flagAcidPhaeo,8),'+m','markersize',5)
loglog(DATA(ilist.flagS2Chl,9),DATA(ilist.flagS2Chl,8),'xb','markersize',8)
loglog(DATA(ilist.dmsext,9),DATA(ilist.dmsext,8),'sb','markersize',8)
loglog(DATA(ilist.dmspt.min,9),DATA(ilist.dmspt.min,8),'+c','markersize',5)
loglog(DATA(ilist.minsumdmsp,9),DATA(ilist.minsumdmsp,8),'og','markersize',8)
axis([0.009 30 0.1 3000])
legend('all data','estuarine','low salinity','shallow bottom',...
    'acidified Phaeo','outside 1-99 quantile regr',...
    'extreme DMS','too low DMSPt','too low DMSP based on Chl',...
    'location','southeast')
xlabel('Chl in situ (mg m^{-3})'); ylabel('DMSPt in situ (nM)')


%% =========== Set values to nan acccording to selected criteria ==========
% ----------------------------- EDIT TO PLAY ------------------------------
DATA(ilist.allflags,7:9) = nan; % all quality control flags
DATA(ilist.sal.min,7:9) = nan;  % salinity   
DATA(ilist.zbot.max,7:9) = nan; % bottom depth
DATA(ilist.minsumdmsp,8) = nan; % DMSPt lower than 10% of expected from phyto classes
DATA(ilist.dmsext,7) = nan;     % extreme DMS
DATA(ilist.dmspt.min,8) = nan;  % too low DMSPt

%% Apply a conversion factor to HPLC chlorophyll 
% Most chl data measured with Turner fluorometer
% Conversion factor applied to HPLC-chl to approximate Turner-chl
% Provides marginal improvement only (check regressions in section below)
% Identify studies by contribution number (contribnum)
ilist.hplc = DATA(:,1)==135 | (DATA(:,1)>158 & (DATA(:,1)<164)) | DATA(:,1)==191 | DATA(:,1)==195;
DATA(ilist.hplc,9) = 1.39*DATA(ilist.hplc,9);

%% "Adjust" in situ chl using its regression against satellite chl (OCX algorithms)
% Useful for the development of predictive statistical models based on in
% situ chl data: removes biases inherent to comparison between in situ vs. satellite chl
% NOTE: here using regression of daily single pixel data only
% Most appropriate because it represents the best matchups
% NOTE: may be improved if non-case2 samples are removed
y = real(log10(nanmean(DATA(:,[87 88]),2)));
X=[ones(size(DATA(:,9))) real(log10(DATA(:,9)))];
[b,bint,r,rint,stats]=regress(y,X);
fprintf('\nIntercept_day_1pix = %0.2f (%0.2f__%0.2f),\nSlope_day = %0.2f (%0.2f__%0.2f),\nR^2 = %0.2f,\np = %0.1e\n',...
    b(1), bint(1,1), bint(1,2),...
    b(2), bint(2,1), bint(2,2),...
    stats(1), stats(3))
regr.adjChl = struct('b',b,'bint',bint,'r',r,'rint',rint,'stats',stats);   % save in structure   
DATA(:,9) = 10.^(regr.adjChl.b(1) + regr.adjChl.b(2)*log10(DATA(:,9)));

%% =======================================================================
%                               FINAL TIPS 
% Tricks used in the Galí et al. 2018 (Biogeosciences) paper to maximize N
% ========================================================================

%% Fill gaps of in situ SST using satellite SST to maximize N
% Satellite SST has small retrieval error compared to bio-topical
% variables. In our datasets, R^2 = 0.99 and RMSE = 0.8C (Galí et al. 2015)
DATA(isnan(DATA(:,10)),10) = DATA(isnan(DATA(:,10)),102);

%% Maximize the amount of satellite data available for algorithm tuning and validation
% by creating merged satellite variables: par, pic, kd490, and other derived variables

% PAR from 8-day MODIS-A and SeaWiFS matchups (1998 to 2012)
% Stepwise addition of 1-pixel matchups, 3x3 + 5x5 matchups, SeaWiFS climatology 
% Use of climatological data for pre-SeaWiFS measurements (<1998) does not
% degrade appreciably the relationships
parmatch = nanmean([DATA(:,55) DATA(:,60)],2); % 
parmatch(isnan(parmatch)) = nanmean([DATA(isnan(parmatch),[65 75]) DATA(isnan(parmatch),[70 80])],2); % 
parmatch(isnan(parmatch)) = DATA(isnan(parmatch),21); % SeaWiFS climatology

% PIC from satellite matchups (8-day directly)
picmatch = nanmean([DATA(:,81) DATA(:,82)],2);
picmatch(isnan(picmatch)) = nanmean([DATA(isnan(picmatch),83) DATA(isnan(picmatch),84)],2);
picmatch(isnan(picmatch)) = nanmean([DATA(isnan(picmatch),85) DATA(isnan(picmatch),86)],2);

% CHL from satellite matchups (1D and 8D stepwise)
chlmatch = nanmean([DATA(:,87) DATA(:,88)],2);
chlmatch(isnan(chlmatch)) = nanmean([DATA(isnan(chlmatch),89) DATA(isnan(chlmatch),90)],2);
chlmatch(isnan(chlmatch)) = nanmean([DATA(isnan(chlmatch),91) DATA(isnan(chlmatch),92)],2);
chlmatch(isnan(chlmatch)) = nanmean([DATA(isnan(chlmatch),51) DATA(isnan(chlmatch),56)],2);
chlmatch(isnan(chlmatch)) = nanmean([DATA(isnan(chlmatch),61) DATA(isnan(chlmatch),66)],2);
chlmatch(isnan(chlmatch)) = nanmean([DATA(isnan(chlmatch),71) DATA(isnan(chlmatch),76)],2);

% Kd from in situ data, satellite matchups (and satellite clim otherwise? removed)
Kdmatch = nanmean([DATA(:,23) DATA(:,53) DATA(:,58)],2);
Kdmatch(isnan(Kdmatch)) = nanmean([DATA(isnan(Kdmatch),63) DATA(isnan(Kdmatch),68)],2);
Kdmatch(isnan(Kdmatch)) = nanmean([DATA(isnan(Kdmatch),73) DATA(isnan(Kdmatch),78)],2);

% Fraction of surface PAR in the ML (fracmld) and mean MLD PAR (parmld)
fracmld = (1./(DATA(:,18).*Kdmatch)).*(1-exp(-DATA(:,18).*Kdmatch));
parmld = parmatch.*fracmld; % matchup par for years >= 1998, clim otherwise
parcmld = DATA(:,21).*fracmld; % fully climatological par

% Append new variables to DATA matrix and head (headers)
DATA = [DATA parmatch chlmatch picmatch Kdmatch fracmld parmld parcmld];
head = [head; {'parmatch' 'chlmatch' 'picmatch' 'Kdmatch' 'fracmld' 'parmld' 'parcmld'}'];
