% This quaitfy method is from Tuttle and Savicci,2016[1].
% This method derived from three steps:
% (1) regression Y = beta*X+c, and got p value.
% attention: this part remove the bootstrap method raised by [1] to save time.
% (2) Use GLM method to seperate dry days and wet days, if anomaly of target day
% is larger than seasonaly anomaly, set as wet day, vice versa.
% (3) Use mean of prediction of full model and baseline model(full model set 
% coefficient of predict variables as 0) to get the qualify value of X-Y
% feedback.
% 
% function [impact,p_value] = quatify(X,resid_X,Y,resid_Y,season_anomaly)
% %
% % qualify method in [1]
% %
% % Parameters
% % __________
% % - X/Y
% % - resid_X/resid_Y
% % - season_anomaly
% %
% % Attributes
% % - impact
% % - p_value
% 
% %
% offsetR = check_rf2glm(Y-resid_Y);
% [b,p_value,dispsave] = reg_rf2glm(resid_X,Y,offsetR);
% const = ones(numel(Y),1); 
% % qualify
% if ~isnan(b(2))
%     %
%     PNaN = Y.*0;
%     % predictions from full model
%     pred_full = normcdf(1.*offsetR+...
%                         (b'*[const,resid_X]')'+...
%                         (PNaN),...
%                         0,dispsave);
%     % predictions from baseline model
%     pred_baseline = normcdf(1.*offsetR+...
%                             (b'*[const,resid_X.*0]')'+...
%                             (PNaN),...
%                             0,dispsave);
%     %
%     quants = quantile(season_anomaly+PNaN-nanmean(season_anomaly+PNaN),...
%                       [0.1 0.25 0.5 0.75 0.9]);    
%     % mean S impact in bottom/top 50% of S anomaly
%     % (dividing pred_full by pred_baseline works because slagperts
%     % has a mean of zero, so it does not change the mean of the prediction)
%     II = find(season_anomaly + PNaN - nanmean(season_anomaly+PNaN)...
%               < quants(ceil(length(quants)/2)));
%     impact_lowS = nanmean(pred_full(II))./...
%                   nanmean(pred_baseline(II));
%     %
%     II = find(season_anomaly + PNaN - nanmean(season_anomaly+PNaN)...
%               > quants(ceil(length(quants)/2)));
%     impact_highS = nanmean(pred_full(II))./...
%                    nanmean(pred_baseline(II));
%     %save mean S impact in bottom/top 50% of S anomaly
%     impact = [impact_lowS,impact_highS];  
% else
%     p_value = NaN;
%     impact = nan(2,1);    
% end
% end
% 
% function [b,p_value,disp] = reg_rf2glm(X,Y,offsetR)
% %
% % function relate random forest regression with glm regression. 
% %
% % Parameters
% % __________
% % - X/Y
% % - offsetR
% %
% % Attributes
% % - b
% % - p_value
% % - disp
% 
% %
% const = ones(numel(Y),1);    %constant term in regression
% % fit
% try   
%     [b,~,stats] = glmfit([const,X],Y,'binomial',...
%                          'link','probit',...
%                          'estdisp','on',...
%                          'constant','off',...
%                          'offset',offsetR);    
%     disp = stats.sfit(end);    %dispersion
%     p_value = stats.p(2);
% catch WRONG   
%     disp('error in unrestrict, uncorrected offset model')
%     disp = NaN;
%     p_value = NaN;
%     b(2) = NaN;
% end
% end
% 
% function X = check_rf2glm(X)
% 
% % random forest changed to glm method, need to set down the offsetR
% % parameters:
% % __________
% % X: ***X must be the prediction array of random foreset regression
% %
% % if x residual close to 0,corresponding offset value close to +-inf.
% % thus this phenomenen should be exclude.
% 
% offsetR = norminv(X);
% offsetR(find(offsetR<-1000 | offsetR>1000)) = NaN;
% 
% end 


%% caculate S-POCC and S-P index.



function [spImpact,p] = quatify(x,y,POCC,S,Sclim_aic,N)

% transform x residual got from best model to offset of GLM.
% random forest changed to glm method, need to set down the offsetR
offsetR = norminv(x);
% if x residual close to 0,corresponding offset value close to +-inf.
% thus this phenomenen should be exclude.
offsetR(find(offsetR<-1000 | offsetR>1000)) = NaN; 

% regression of POCC residual = beta*S residual+c.
one=ones(N,1);    %constant term in regression

try
    
    [b,~,stats,~,iter_flag,~,scale_flag]=glmfit_itersave...
    ([one,S-y],POCC,'binomial','link','probit','estdisp','on','constant','off','offset',offsetR);
    
    dispsave=stats.sfit(end);    %dispersion
    bcoefsave=b;    %regression coefficients
    sesave_all=stats.se;    %standard error
    pvalsave_all=stats.p;    %coefficient significance
    
    yPredict = glmval(b,[one,S-y],'probit','constant','off');
    sSSE = sum((yPredict-mean(POCC)).^2);
    mm = (POCC-offsetR);
    sSSE = sum((mm-mean(mm)).^2)-sum((mm-yPredict).^2);
    
    
    if any([iter_flag,scale_flag])>0
        dispsave=NaN;
        bcoefsave([1,end-4:end])=NaN;
        sesave_all([1,end-4:end])=NaN;
        pvalsave_all([1,end-4:end])=NaN;
    end
    
catch ME
    
    disp('error in unrestrict, uncorrected offset model')
    dispsave=NaN;
    bcoefsave([1,end-4:end])=NaN;
    sesave_all([1,end-4:end])=NaN;
    pvalsave_all([1,end-4:end])=NaN;
    
end



% ------------------------------------------------------------------------------------------------
% use predict value of restricted model/full model to qualify S-POCC impact.
% ------------------------------------------------------------------------------------------------

if ~isnan(b(2))
    
    p = pvalsave_all(2);
    PNaN=POCC.*0;
    
    %predictions from unrestricted model - all coefficients
    %corrected for bias
    predU_all=normcdf(1.*offsetR+(squeeze(bcoefsave)'*...
        [one,S-y]')'+(PNaN),0,dispsave);
    
    %predictions from unrestricted model, but with the soil
    %moisture term set to zero - all coefficients corrected for bias
    predU_all_sp0=normcdf(1.*offsetR+(squeeze(bcoefsave)'*...
        [one,(S-y).*0]')'+(PNaN),0,dispsave);
    
%     disp(nanmean(predU_all))
%     disp(nanmean(predU_all_sp0))
    quants=quantile(Sclim_aic+PNaN-nanmean(Sclim_aic+PNaN),[0.1 0.25 0.5 0.75 0.9]);
    
    %mean S impact in bottom/top 50% of S anomaly
    %(dividing predU_all by predU_all_sp0 works because slagperts
    %has a mean of zero, so it does not change the mean of the
    %prediction)
    II=find(Sclim_aic+PNaN-nanmean(Sclim_aic+PNaN)<quants(ceil(length(quants)/2)));
    P_prob_div_lowS_U_allR=nanmean(predU_all(II))./nanmean(predU_all_sp0(II));
%     disp(nanmean(predU_all(II)))
%     disp(nanmean(predU_all_sp0(II)))
    II=find(Sclim_aic+PNaN-nanmean(Sclim_aic+PNaN)>quants(ceil(length(quants)/2)));
    P_prob_div_highS_U_allR=nanmean(predU_all(II))./nanmean(predU_all_sp0(II));
%     disp(nanmean(predU_all(II)))
%     disp(nanmean(predU_all_sp0(II)))
    %save mean S impact in bottom/top 50% of S anomaly
    spImpact=[P_prob_div_lowS_U_allR,P_prob_div_highS_U_allR];
    
else
    
    p = NaN;
    spImpact = nan(2,1);
    
end

end

