function [impact,p_value] = quatify(X,Y)


end



function [p,spImpact,sSSE,predU_all,predU_all_sp0] = spImpactGLM(x,y,POCC,S,Sclim_aic,N)

% ------------------------------------------------------------------------------------------------
% transform x residual got from best model to offset of GLM.
% ------------------------------------------------------------------------------------------------

offsetR = norminv(x);
% if x residual close to 0,corresponding offset value close to +-inf.
% thus this phenomenen should be exclude.
offsetR(find(offsetR<-1000 | offsetR>1000)) = NaN; 

% ------------------------------------------------------------------------------------------------
% regression of POCC residual = beta*S residual+c.
% ------------------------------------------------------------------------------------------------
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



function Sclim_aic = getAnomaly(maxseasterms,seasmod,S,N,AICPEN)



maxind_S=2.^(maxseasterms)-1;
index_S=dec2bin(1:maxind_S);
index_S = index_S == '1';

MAP_IND_S=nan(1,2*maxseasterms);

for jj=1:maxseasterms
    MAP_IND_S((jj-1)*2+1:(jj-1)*2+2)=jj;
end

model_params_Sseas=zeros(1+2*maxseasterms,1);
bvals_Sseas=nan(1+2*maxseasterms,1);
pvals_Sseas=nan(1+2*maxseasterms,1);


aicyrandseas_S=nan(maxind_S,1);
likeyrandseas_S=nan(maxind_S,1);

iter_save_S=nan(maxind_S,1);
iter_flag_save_S=nan(maxind_S,1);
rank_flag_save_S=nan(maxind_S,1);
scale_flag_save_S=nan(maxind_S,1);

one=ones(N,1);

for kk=1:maxind_S
    
    clear II PP
    II=find(index_S(kk,:)==1);
    
    PP=zeros(1,length(MAP_IND_S));
    for ii=1:length(II)
        JK=find(MAP_IND_S==II(ii)); %+1 is b/c of constant
        PP(JK)=1;
    end
    PP=find(PP==1);
    
    [~,dev,~,iter,iter_flag,rank_flag,scale_flag]=...
        glmfit_itersave([one,seasmod(:,PP)],S,'normal',...
        'link','identity','estdisp','on','constant','off');
    
    lenPP=1+length(PP);
    aicyrandseas_S(kk,1)=AICPEN*(lenPP+1)-2*(-dev/2);
    
    likeyrandseas_S(kk,1)=(-dev/2);
    
    iter_save_S(kk,1)=iter;
    iter_flag_save_S(kk,1)=iter_flag;
    rank_flag_save_S(kk,1)=rank_flag;
    scale_flag_save_S(kk,1)=scale_flag;
end

II=find(iter_flag_save_S==1 | rank_flag_save_S==1 | ...
    scale_flag_save_S==1);
aicyrandseas_S(II)=NaN;

[~,bind_S]=nanmin(aicyrandseas_S);  %find min aic

clear II PP

II=find(index_S(bind_S,:)==1);

PP=zeros(1,length(MAP_IND_S));
for ii=1:length(II)
    JK=find(MAP_IND_S==II(ii));
    PP(JK)=1;
end
PP=find(PP==1);

[b_Sseas,~,stats_Sseas,iter,iter_flag,rank_flag,scale_flag]=...
    glmfit_itersave([one,seasmod(:,PP)],S,'normal',...
    'link','identity','estdisp','on','constant','off');

model_params_Sseas([1,PP+1])=1;
bvals_Sseas([1,PP+1])=b_Sseas;
pvals_Sseas([1,PP+1])=stats_Sseas.p;

if any([iter_flag,rank_flag,scale_flag])>0
    model_params_Sseas([1,PP+1])=nan(length(seasmod(:,1))+1);
    bvals_Sseas([1,PP+1])=nan(length(seasmod(:,1))+1);
    pvals_Sseas([1,PP+1])=nan(length(seasmod(:,1))+1);
end

%climatology as determined by best AIC seasonal model of S
sclim_9yr_aic=[ones(length(seasmod(:,1)),1),seasmod(:,PP)]*b_Sseas;

%subtract mean S climatology from S
Sclim_aic=S-sclim_9yr_aic;
end


end



%% caculate S-POCC and S-P index.


function [R2,devR2,indexRank] = rank...
            (inPeriod,inLaggedP,inPress,inSpatialP,inSpatialPress,y,NumTrees)

% use input terms to construm full predictor matrix, [] use for construct full matrix in loop. 
inLaggedP = [inLaggedP,inSpatialP];
inPress = [inPress,inSpatialPress];
in = {inPeriod;inLaggedP;inPress;[]};
%in = {inPeriod;inLaggedP;inLaggedS;inPress;[]};


NITER = size(in); %number of terms.

% ------------------------------------------------------------------------------------------------
% y = a(x,w,e,f...)+b; get corresponding R2; exclude x,w,e,f,etc,respectively; get R2.
% contrast these R2 value to find the impact of each terms. i.e.,x,w,e,f,etc.
% ------------------------------------------------------------------------------------------------

for ii = 1:NITER
    
    II = 1:NITER;
    II(ii) = [];
    x = cat(2,in{II});
    
    coeff =  abs(corr(x(find(~isnan(y)),:),y(find(~isnan(y)))));
    xNew = x(:,find(~isnan(coeff)));
    
    B = TreeBagger(NumTrees,xNew,y,'method','regression','Surrogate','on',...
        'OOBPredictorImportance','on','PredictorSelection','curvature');    
    R2(ii) = rSquared(predict(B,xNew),y);
    
end

    devR2 = R2(end)-R2(1:end-1);
    [~,indexRank] = sort(devR2,'descend'); % indexRank means the importance array of predictors.

end
