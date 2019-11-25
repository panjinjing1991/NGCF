% This quaitfy method is from Tuttle and Savicci,2016[1].
% Their method derived from three steps:
% (1) regression Y = beta*X+c, and got p value.
% attention: this part remove the bootstrap method raised by [1] to save time.
% (2) Use GLM method to seperate dry days and wet days, if anomaly of target day
% is larger than seasonaly anomaly, set as wet day, vice versa.
% (3) Use mean of prediction of full model and baseline model(full model set 
% coefficient of predict variables as 0) to get the qualify value of X-Y
% feedback.

function [impact,p_value] = quatify(X,Y)
%
%
%
%
%
% qualify
if ~isnan(b(2))
    %
    PNaN = Y.*0;
    % predictions from full model
    pred_full = normcdf(1.*offsetR+...
                        (b'*[const,Y]')'+...
                        (PNaN),...
                        0,dispsave);
    % predictions from baseline model
    pred_baseline = normcdf(1.*offsetR+...
                            (b'*[const,Y.*0]')'+...
                            (PNaN),...
                            0,dispsave);
    %
    quants = quantile(season_anomaly+PNaN-nanmean(season_anomaly+PNaN),...
                      [0.1 0.25 0.5 0.75 0.9]);    
    % mean S impact in bottom/top 50% of S anomaly
    % (dividing pred_full by pred_baseline works because slagperts
    % has a mean of zero, so it does not change the mean of the prediction)
    II = find(season_anomaly + PNaN - nanmean(season_anomaly+PNaN)...
              < quants(ceil(length(quants)/2)));
    impact_lowS = nanmean(pred_full(II))./...
                  nanmean(pred_baseline(II));
    %
    II = find(season_anomaly + PNaN - nanmean(season_anomaly+PNaN)...
              > quants(ceil(length(quants)/2)));
    impact_highS = nanmean(pred_full(II))./...
                   nanmean(pred_baseline(II));
    %save mean S impact in bottom/top 50% of S anomaly
    impact = [impact_lowS,impact_highS];  
else
    p_value = NaN;
    impact = nan(2,1);    
end
end

function [b,p_value,disp] = reg_rf2glm(X,Y)
% function relate random forest regression with glm regression. 
%
%
% 
[~,predict,~,~] = RF(X,Y);
offsetR = check_rf2glm(predict);
%
const = ones(numel(Y),1);    %constant term in regression
% fit
try   
    [b,~,stats] = glmfit([const,X],Y,'binomial',...
                         'link','probit',...
                         'estdisp','on',...
                         'constant','off',...
                         'offset',offsetR);    
    disp = stats.sfit(end);    %dispersion
    p_value = stats.p(2);
catch WRONG   
    disp('error in unrestrict, uncorrected offset model')
    disp = NaN;
    p_value = NaN;
    b(2) = NaN;
end
end

function season_anomaly = get_season_anomaly(nSeasonal,X)
% seasonal anomaly by method from Tuttle and Savincci[1]

%
aic_penalty = 2;
N = numel(X);
%
season_terms = terms.get_period_terms();
% index for all possible regression
max_dec_index = 2^nSeasonal-1;
index = dec2bin(1:max_dec_index);
index = index == '1';
%
map_index = nan(2*nSeasonal,1);
for j = 1:nSeasonal
    map_index((j-1)*2+1:(j-1)*2+2) = j;
end
%
const = ones(N,1);
%
for i = 1:max_dec_index
    II = find(index(i,:)==1);
    
    JJ = zeros(length(map_index),1);
    for j = 1:length(II)
        JK = find(map_index==II(i));
        JJ(JK) = 1;
    end
    JJ = find(JJ==1);
    
    [~,dev] = glmfit([const,season_terms(:,JJ)],X,...
                     'normal',...
                     'link','identity',...
                     'estdisp','on',...
                     'constant','off');
    aic(i) = AIC(dev,aic_penalty);         
end
%find min aic
[~,ind]=nanmin(aic);
% 
II = find(index(ind,:)==1);
KK = zeros(length(map_index),1);
for i = 1:length(II)
    JK = find(map_index==II(i));
    KK(JK) = 1;
end
KK = find(KK==1);
% fit best model
[b,~,stats] = glmfit([const,season_terms(:,KK)],X,...
                     'normal',...
                     'link','identity',...
                     'estdisp','on',...
                     'constant','off');
% climatology as determined by best AIC seasonal model of S
clim_ = [const,season_terms(:,KK)]*b;
% subtract mean S climatology from S
season_anomaly = X - clim_;

end

function X = check_rf2glm(X)

% random forest changed to glm method, need to set down the offsetR
% parameters:
% __________
% X: ***X must be the prediction array of random foreset regression
%
% if x residual close to 0,corresponding offset value close to +-inf.
% thus this phenomenen should be exclude.
offsetR = norminv(X);
offsetR(find(offsetR<-1000 | offsetR>1000)) = NaN;

end 
