
function R2 = models(fun)
end


function check_terms(X)

% use pearson coefficent to select top feature.

coeff =  abs(corr(x(find(~isnan(y)),:),y(find(~isnan(y)))));
% NaN value is at end of sort array.
[~,indexDescend] = sort(coeff,'descend','MissingPlacement','last'); 

% number of NaN coefficent, which means sum of the predictor variable is 0.
IJ = find(~isnan(coeff));

% select the top AA feature by pearson coefficent.
selectLinearIndex = indexDescend(1:round(AA*length(IJ)));

% ------------------------------------------------------------------------------------------------
% use random forest to select top feature.
% ------------------------------------------------------------------------------------------------

% construct the new input predictor matrix of Random forest.
xNew = x(:,IJ);

end

function R2 = RF(X,Y,param)

if nargin==2
    % use complete predictor matrix, and get corresponding
    % predict array;Rsquared;Important array of predictor;MSE
    B = TreeBagger(NumTrees,X,Y,'method','regression','Surrogate','on',...
        'OOBPredictorImportance','on','PredictorSelection','curvature');

    importanceArray = B.OOBPermutedPredictorDeltaError; 
    [~,indexDescend] = sort(importanceArray,'descend');
    R2 = R_Squared(predict(B,X),Y);

else
    [R2,best_param] = grid_search_sv(param,3);

end
end

function R2 = ANN(X,Y,param)

if nargin==2
    net = feedforwardnet(10);
    net = train(net,X,Y);
    R2 = R_Squared(sim(net,X),Y);
else
    [R2,best_param] = grid_search_cv(param,1);
end

end

function R2 = SVR(X,Y,param)

if nargin==2
    % fit SVM regression
    svr = fitrsvm(X,Y,'Standardize',true,'kernelfunction','gaussian');
    R2 = R_Squared(svr.resubPredict, Y);
else
    [R2,best_param] = grid_search_cv(param,2);
end
% plot y and ypredict
% plot(svr.resubPredict);hold on;plot(Y);
% legend('Response in Training data','Predicted Response','location','best');
% Estimate mse and epsilon-insensitive loss by cross-validation
% cv = crossval(svr);
% mse = kfoldLoss(cv);
% epsLoss = kfoldLoss(cv,'lossfun','epsiloninsensitive');
end

function [R2,best_param] = grid_search_cv(param,fun)

% grid search cv for each input machine learning method
%
% Parameters:
% __________
% param:
% 1. BP Artificial Neural Network
%    (cell)
%    - hidlaysize1/hidlaysize2 :[15 30 70][10 20 50]
%    - max epoch:[10 20 40 90]
%    - transfer function:{'logsig' 'tansig'}
%    - train option:{'traingd' 'traingda' 'traingdm' 'traingdx'}
% 2. Support vector machine regression 
%    (array as ['linear','polynomial','rbf'])
%    - kernel_function : the name of the kerenl that is to be used, this
%                        variable must be a string    
% 3. Random forest
%    - 
% fun:
% 1. BP Artificial Neural Network
% 2. Support vector machine regression
% 3. Random forest
%
% Attributes:
% __________
% R2: 
% best_param:

if fun==1   
    % set all pairs of input param.
    [p,q,r,s,t] = ndgrid(param{1},param{2},param{3},...
                         1:length(param{4}),1:length(param{5}));
    pairs = [p(:) q(:) r(:) s(:) t(:)];
    R2_ = zeros(size(pairs,1),1);
    % loop for grid search
    for i = 1:size(pairs,1)
      % set parameters according the input param.
      setdemorandstream(672880951)
      net = patternnet([pairs(i,1) pairs(i,2)]);
      net.trainParam.epochs = pairs(i,3);
      net.trainFcn	= param{4}{pairs(i,4)};
      net.layers{2}.transferFcn = param{5}{pairs(i,5)};
%       net.divideParam.trainRatio = 0.9;
%       net.divideParam.valRatio = 0.1;
%       net.divideParam.testRatio = 0; 
      % fit BP and get R2 score.
      net = train(net,X,Y);
      R2_(i) = R_Squared(sim(net,X),Y);
    end
    % get max R2
    [R2,ind] = max(R2_);
    best_param = {pairs(ind,1) pairs(ind,2) pairs(ind,3) ...
        param{4}{pairs(ind,4)} param{5}{pairs(i,5)}};
    
elseif fun==2
    R2_ = zeros(numel(param),1);
    % fit SVM regression
    tic
    for i = 1:numel(param)
    svr = fitrsvm(...
                  X,Y,'Standardize',true, ...
                  'kernelfunction',param(i),...
                  'OptimizeHyperparameters','auto',...
                  'HyperparameterOptimizationOptions',...
                  struct('AcquisitionFunctionName',...
                  'expected-improvement-plus')); 
    % fit SVR and get R2 score.
    R2_(i) = R_Squared(svr.resubPredict, Y); 
    end
    toc
    % get max R2
    [R2,ind] = max(R2_);
    best_param = param(ind);
    
elseif fun==3
    'Need improved'   
end
end



function [yPredict,selectIndex,R2,pSSE] = combineRandomForestBestModel...
         (x,y,AA,BB,NumTrees)



% use complete predictor matrix, and get corresponding
% predict array;Rsquared;Important array of predictor;MSE
B = TreeBagger(NumTrees,xNew,y,'method','regression','Surrogate','on',...
    'OOBPredictorImportance','on','PredictorSelection','curvature');

importanceArray = B.OOBPermutedPredictorDeltaError; 
[~,indexDescend] = sort(importanceArray,'descend');

R2All = rSquared(predict(B,xNew),y);

% ------------------------------------------------------------------------------------------------
% union the result of RF and Linear.
% ------------------------------------------------------------------------------------------------

selectRfIndex = IJ(indexDescend(1:round(BB*length(IJ))));
selectIndex = union(selectLinearIndex,selectRfIndex);

% use union predictor matrix from linear and RF, and get R^2.
B = TreeBagger(NumTrees,x(:,selectIndex),y,'method','regression','Surrogate','on',...
    'OOBPredictorImportance','on','PredictorSelection','curvature');

yPredict = predict(B,x(:,selectIndex));

% %pSSE = sum((yPredict-mean(y)).^2);
pSSE = sum((y-mean(y)).^2)-sum((y-yPredict).^2);


R2Best = rSquared(predict(B,x(:,selectIndex)),y);

R2 = [R2All,R2Best];

end
