function [impact,p_value] = quatify(X,Y)

% if x residual close to 0,corresponding offset value close to +-inf.
% thus this phenomenen should be exclude.
offsetR = norminv(X);
offsetR(find(offsetR<-1000 | offsetR>1000)) = NaN; 
%
const = ones(numel(Y),1);    %constant term in regression
% fit
try   
    [b,~,stats] = glmfit([const,X],Y,'binomial',...
                         'link','probit',...
                         'estdisp','on',...
                         'constant','off',...
                         'offset',offsetR);    
    dispsave = stats.sfit(end);    %dispersion
    p_value = stats.p(2);
catch WRONG   
    disp('error in unrestrict, uncorrected offset model')
    dispsave = NaN;
    p_value = NaN;
    b(2) = NaN;
end
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
    quants = quantile(Sclim_aic+PNaN-nanmean(Sclim_aic+PNaN),...
                      [0.1 0.25 0.5 0.75 0.9]);    
    % mean S impact in bottom/top 50% of S anomaly
    % (dividing pred_full by pred_baseline works because slagperts
    % has a mean of zero, so it does not change the mean of the prediction)
    II = find(Sclim_aic + PNaN - nanmean(Sclim_aic+PNaN)...
              < quants(ceil(length(quants)/2)));
    impact_lowS = nanmean(pred_full(II))./...
                  nanmean(pred_baseline(II));
    %
    II = find(Sclim_aic + PNaN - nanmean(Sclim_aic+PNaN)...
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

function clim = get_season_anomaly(nSeasonal,X)
% seasonal anomaly by method from Tuttle and Savincci[1]

%
AICPEN = 2;
N = numel(X);
%
index = dec2bin(1:(2^nSeasonal-1));
index = index == '1';
%
const = ones(N,1);
for i = 1:length
    II = find(index(i,:)==1);
    
end

end
