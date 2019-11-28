% This main.m file is apply nonlinear granger causality framework(NGCF) on
% identify surface soil moisture-precipitation feedback, including both sign 
% and pattern distribution.

% Procedure:
% 1. Use get_terms.m to get independent terms, as period terms(i.e., inter-
% annual cycle, seasonal cycle); lagged terms(i.e., lagged P and press); spatial
% terms(i.e., lagged P and press over selected square indicated spatial impact)
%
% 2. Run random forest for P occurrence and surface soil moisture as dependent 
% terms, and independent terms metioned in 1., Addtionally, by using hybrid
% selection feature method to find the 'best' regression(to avoid overfitting 
% in some content). After applying these regression, get the residual value of 
% both surface soil moisture and precipitation.
%
% 3. Then, fit "slagperts" to residual of POCC and "offset", using only S 
% with no P on the previous day, and adding lagged atmospheric pressure as 
% independent variable. 
%
% 4. Use block bootstrap to eliminate endogeneity bias and determine signifi-
% cance of the S coefficient in the regression.(Need improved)
%
% 5. Calculate S-POCC impacts by dividing the unrestricted model by the re-
% stricted model, plotted against the seasonal S anomaly, and taking the mean 
% above and below the median of Sclim.
%
% Paramters:
% _________
% - P,S,press:
% - lon,lat:
% - startDate:
% - nAnnual,nSeasonal
% - day_lag:
% - Slen:
% - pbcrit:
%
% Attributes:
% __________
% - impact:
% - p_value:

function [impact,p_value,varargout] = main(P,S,press,...
                                           lon,lat,...
                                           startDate,...
                                           nAnnual,nSeasonal,day_lag,...
                                           Slen,...
                                           pbcrit) 
%                                      
% ***attention: only for S-P feedback and target dataset.                                  
%% class
terms = get_terms();
model = models();
%% get independent and dependent terms
% P
[depend_P_terms,lagged_P_terms,annualTerm,seasTerm,spatial_P_terms] = ...
terms.get_all_terms(P,startDate,lat,lon,Slen,nAnnual,nSeasonal,day_lag);
% press
[~,lagged_press_terms,~,~,spatial_press_terms] = ...
terms.get_all_terms(press,startDate,lat,lon,Slen,nAnnual,nSeasonal,day_lag);
% get independent terms
independ_terms = [lagged_P_terms,lagged_press_terms,...
                  annualTerm,seasTerm,...
                  spatial_P_terms,spatial_press_terms];
% get dependent terms             
POCC = 0.*depend_P_terms;
POCC(find(depend_P_terms>pbcrit))=1;
S = terms.get_all_terms(S,startDate,113-lat,lon,Slen,nAnnual,nSeasonal,day_lag);
%% 
valid = length(find(~isnan(S)==1))>= round(0.1*length(S)) && ...
        sum(POCC)>0 && sum(S)~=0;
% models
if valid 
% remove independent terms impact
[R2P,residP] = model.run_models(independ_terms,POCC,3);
[R2S,residS] = model.run_models(independ_terms,S,3);
% quatify   
    season_anomaly = terms.get_season_anomaly(nSeasonal,S,seasTerm);
    [impact,p_value] = quatify(S,residS,POCC,residP,season_anomaly);
else
    impact=[nan,nan];p_value=nan;
    R2P=nan;R2S=nan;residP=nan;residS=nan;season_anomaly=nan;
end
% optional output
if nargout>2; varargout{1} = R2P; end
if nargout>3; varargout{2} = R2S; end
if nargout>4; varargout{3} = season_anomaly; end
if nargout>5; varargout{4} = POCC; end
if nargout>6; varargout{5} = S; end
if nargout>7; varargout{6} = independ_terms; end
if nargout>8; varargout{7} = residP; end
if nargout>9; varargout{8} = residS; end

end
