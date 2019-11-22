
function [R2,resid,varargout] = models(X,Y,fun,type,pl,pnl)
% fit X to Y by machine learning models
% and adopt two method to handle overfitting.
% [1] grid search cv.
% [2] hybrid feature selection method.

% remove nan row of X,Y
[X,Y] = check_terms(X,Y);
if nargin==3
    % get param for grid search cv.
    param = get_param(fun);
    % ANN
    if fun==1
        [R2,resid,best_param] = ANN(X,Y,param);
    % SVR    
    elseif fun==2
        [R2,resid,best_param] = SVR(X,Y,param);
    % RF
    elseif fun==3
        [R2,resid,best_param] = RF(X,Y,param);
    end
    varargout{1} = best_param;
elseif nargin==4
    % ANN
    if fun==1
        [R2,resid] = ANN(X,Y);
    % SVR    
    elseif fun==2
        [R2,resid] = SVR(X,Y);
    % RF
    elseif fun==3
        [R2,resid] = RF(X,Y);
    end
elseif nargin>4
    % get feature by hybrid feature selection method.
    select_index = select_feature(X,Y,pl,pnl);
    % ANN
    if fun==1
        [R2,resid] = ANN(X(:,select_index),Y);
    % SVR    
    elseif fun==2
        [R2,resid] = SVR(X(:,select_index),Y);
    % RF
    elseif fun==3
        [R2,resid] = RF(X(:,select_index),Y);
    end
end

end

function [X,Y] = check_terms(X,Y)

% check y nan terms
nan_index_Y = find(isnan(Y));
% check x nan terms
[nan_index_X,~] = find(isnan(X));
% 
nan_index_X = unique(nan_index_X);
nan_index = union(nan_index_X,nan_index_Y);
%
X(nan_index,:) = []; Y(nan_index) = [];

end

function param = get_param(fun)
if fun==1
    param = {[15 30 70],...
             [10 20 50],...
             [10 20 40 90],...
             {'logsig' 'tansig'},...
             {'traingd' 'traingda' 'traingdm' 'traingdx'}};
elseif fun==2
    param = {'linear','polynomial','rbf'};
elseif fun==3
end
end

function index = kfold_order(N,fold_number)
%
T = round(N/fold_number);
%
index = 1:N;
for i = 1:fold_number
    index((i-1)*T+1:i*T) = i;
    if i==fold_number
        index((fold_number-1)*T+1:end) = fold_number;
    end
end
end

function [R2,resid,varargout] = RF(X,Y,param)

if nargin==2
    % use complete predictor matrix, and get corresponding
    % predict array;Rsquared;Important array of predictor;MSE
    NumTrees = 100;
    B = TreeBagger(NumTrees,X,Y,'method','regression','Surrogate','on',...
        'OOBPredictorImportance','on','PredictorSelection','curvature');

    importanceArray = B.OOBPermutedPredictorDeltaError; 
    [~,indexDescend] = sort(importanceArray,'descend');
    resid = Y-predict(B,X);
    R2 = R_Squared(predict(B,X),Y);
    varargout{1} = indexDescend;

else
    [R2,resid,best_param] = grid_search(X,Y,param,3);
    varargout{1} = best_param;
end
end

function [R2,resid,varargout] = ANN(X,Y,param)

if nargin==2
    net = feedforwardnet(10);
    net = train(net,X,Y);
    R2 = R_Squared(sim(net,X),Y);
else
    [R2,resid,best_param] = grid_search(X,Y,param,1);
     varargout{1} = best_param;
end

end

function [R2,resid,varargout] = SVR(X,Y,param)

if nargin==2
    % fit SVM regression
    svr = fitrsvm(X,Y,'Standardize',true,'kernelfunction','gaussian');
    R2 = R_Squared(svr.resubPredict, Y);
else
    [R2,resid,best_param] = grid_search(X,Y,param,2);
     varargout{1} = best_param;
end
% plot y and ypredict
% plot(svr.resubPredict);hold on;plot(Y);
% legend('Response in Training data','Predicted Response','location','best');
% Estimate mse and epsilon-insensitive loss by cross-validation
% cv = crossval(svr);
% mse = kfoldLoss(cv);
% epsLoss = kfoldLoss(cv,'lossfun','epsiloninsensitive');
end

function [R2,resid,best_param] = grid_search(X,Y,param,fun)

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

% ANN
if fun==1   
    % set all pairs of input param.
    [p,q,r,s,t] = ndgrid(param{1},param{2},param{3},...
                       1:length(param{4}),1:length(param{5}));
    pairs = [p(:) q(:) r(:) s(:) t(:)];
    R2_ = zeros(size(pairs,1),1);
    resid_ = zeros(numel(Y),size(pairs,1));
    % loop for grid search
    for i = 1:size(pairs,1)
        % set parameters according the input param.
        setdemorandstream(672880951)
        net = patternnet([pairs(i,1) pairs(i,2)]);
        net.trainParam.epochs = pairs(i,3);
        net.layers{2}.transferFcn	= param{4}{pairs(i,4)};
        net.trainFcn = param{5}{pairs(i,5)};
        net.divideParam.trainRatio = 0.9;
        net.divideParam.valRatio = 0.1;
        net.divideParam.testRatio = 0; 
        net.trainParam.showWindow = 0;
        % fit BP and get R2 score.
        net = train(net,X',Y');
        R2_(i) = R_Squared(sim(net,X'),Y');
        resid_(:,i) = Y'-sim(net,X');
        size(resid_(:,i))
    end
    % get max R2
    [R2,ind] = max(R2_);
    best_param = {pairs(ind,1) pairs(ind,2) pairs(ind,3) ...
        param{4}{pairs(ind,4)} param{5}{pairs(i,5)}};
    resid = resid_(:,ind);
% SVR  
elseif fun==2
    R2_ = zeros(numel(param),1);
    resid_ = zeros(numel(Y),numel(param));
    % fit SVM regression
    tic
    for i = 1:numel(param)
        svr = fitrsvm(...
                      X,Y,'Standardize',true, ...
                      'kernelfunction',param{i},...
                      'OptimizeHyperparameters','auto',...
                      'HyperparameterOptimizationOptions',...
                      struct('AcquisitionFunctionName',...
                      'expected-improvement-plus')); 
        % fit SVR and get R2 score.
        R2_(i) = R_Squared(svr.resubPredict, Y);
        resid_(:,i) = Y-svr.resubPredict;
    end
    toc
    % get max R2
    [R2,ind] = max(R2_);
    best_param = param(ind);
    resid = resid_(:,ind);
% RF
elseif fun==3
    'Need improved'   
end
end

function select_index = select_feature(x,y,pl,pnl)

% use pearson coefficent to select top feature.
coeff =  abs(corr(x,y));
% NaN value is at end of sort array.
[~,linear_index] = sort(coeff,'descend','MissingPlacement','last'); 
% number of NaN coefficent, which means sum of the predictor variable is 0.
nonan = find(~isnan(coeff));
% select the top AA feature by pearson coefficent.
select_linear_index = linear_index(1:round(pl*length(nonan)));

% use random forest to select top feature.
[~,~,nonlinear_index] = RF(x(:,nonan),y);
select_nonlinear_index = nonan(nonlinear_index(1:round(pnl*length(nonan))));
% combine the selected feature
select_index = union(select_linear_index,select_nonlinear_index);

end
