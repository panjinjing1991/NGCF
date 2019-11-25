% This function give several machine learning methods for regression, i.e., 
% random forest, support vector machine, artificial neural network. and for 
% these three methods, raise two method to avoid overfitting. i.e., grid search
% [1] and hybrid feature selection method[2].

% Example: model = models();
%          model.run_models(X,Y,fun,type,pl,pnl);

% Referrence:
% [1]:
% [2]:

% Copyright(c) Li Lu, 2019

function model = models

model.run_models = @run_models;
model.RF = @RF;
model.SVR = @SVR;
model.ANN = @ANN;
model.grid_search = @grid_search;
model.select_feature = @select_feature;

end

function [R2,resid,varargout] = run_models(X,Y,fun,~,pl,pnl)
% fit X to Y by machine learning models
% and adopt two method to handle overfitting.
% [1] grid search.
% [2] hybrid feature selection method.

% Parameters:
% __________
% X/Y: target and predict matrix   
% fun:
    % 1. BP Artificial Neural Network
    % 2. Support vector machine regression
    % 3. Random forest
% pl/pnl: percent of feature selected by linear and nonlinear method.

% Attributes:
% __________
% R2:
% resid: resid array of selected regression method.
% varargout: if nargin==4, return 'best_param' of grid search.

% Example:
% 1. none 
%    run_models(X,Y,fun);
% 2. grid search
%    run_models(X,Y,fun,~);
% 3. hybrid feature selection method
%    run_models(X,Y,fun,~,pl,pnl);

% remove nan row of X,Y
[X,Y] = check_terms(X,Y);
if nargin==3
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
elseif nargin==4
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

%%
function [X,Y] = check_terms(X,Y)
% remove nan rows both in target and predict matrix

% check y nan terms
nan_index_Y = find(isnan(Y));
% check x nan terms
[nan_index_X,~] = find(isnan(X));
% selected x nan rows.
nan_index_X = unique(nan_index_X);
% union x,y nan rows.
nan_index = union(nan_index_X,nan_index_Y);
% remove nan rows in both x,y
X(nan_index,:) = []; Y(nan_index) = [];

end

function param = get_param(fun)
% get parameter cells of three machine learning methods.
% 1. BP Artificial Neural Network (cell)
%    - hidlaysize1/hidlaysize2 :[15 30 70][10 20 50]
%    - max epoch:[10 20 40 90]
%    - transfer function:{'logsig' 'tansig'}
%    - train option:{'traingd' 'traingda' 'traingdm' 'traingdx'}
% 2. Support vector machine regression (cell)
%    - kernel_function : the name of the kerenl that is to be used, this
%                        variable must be a string    
% 3. Random forest
%    - 

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
% split data into train set and test set in sequence.

% Parameters:
% __________
% N: length of data;  
% fold_number:

% Attributes:
% __________
% index: testset index of specific fold_number, 1,2,...,
%        fold_number in index means trainset index,
%        the rest index means testset. thus a index array
%        could construct fold_number pairs of trainset/testset.

% ***attention: only for those data couldn't be messed up.

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

%% machine learning regression function.
%
% Parameters
% __________
% X/Y: predict,target variables.
% param: parameters cell for grid search algorithm
%
% Attributes:
% __________
% R2: determinate coefficient
% resid: 
% varargout: 

% Example:
% 1. none
%    [R2,resid] = fun(X,Y);
% 2. grid search
%    [R2,resid,best_param] = fun(X,Y,param);

function [R2,resid,varargout] = ANN(X,Y,param)

if nargin==2
    % get net
    net = feedforwardnet(10);
    % train
    net = train(net,X,Y);
    % R2 and resid
    R2 = index(sim(net,X),Y,1);
    resid = Y-sim(net,X);
else
    [R2,resid,best_param] = grid_search(X,Y,param,1);
     varargout{1} = best_param;
end

end

function [R2,resid,varargout] = SVR(X,Y,param)

if nargin==2
    % fit SVM regression
    svr = fitrsvm(X,Y,...
                  'Standardize',true,...
                  'kernelfunction','gaussian');
    % R2 and resid
    R2 = index(svr.resubPredict,Y,1);
    resid = Y-svr.resubPredict;
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

function [R2,resid,varargout] = RF(X,Y,param)
% random forest regression

if nargin==2
    % radom forest class
    NumTrees = 100;
    B = TreeBagger(NumTrees,...
                   X,Y,...
                   'method','regression',...
                   'Surrogate','on',...
                   'OOBPredictorImportance','on',...
                   'PredictorSelection','curvature');
    % get important array
    importanceArray = B.OOBPermutedPredictorDeltaError; 
    [~,index_descend] = sort(importanceArray,'descend');
    varargout{1} = index_descend;
    % get r2 and resid
    resid = Y-predict(B,X);
    R2 = R_Squared(predict(B,X),Y);
else
    [R2,resid,best_param] = grid_search(X,Y,param,3);
    varargout{1} = best_param;
end

end

%% avoid overfitting methods
function [R2,resid,best_param] = grid_search(X,Y,param,fun)

% grid search cv for each input machine learning method
%
% Parameters:
% __________
% param:
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
        R2_(i) = index(sim(net,X'),Y',1);
        resid_(:,i) = Y'-sim(net,X');
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
        R2_(i) = index(svr.resubPredict,Y,1);
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
% hybrid feature selection method
% select useful information(i.e.,feature) by appling pearson corr and random
% forest in order, and get important array of both linear and nonlinear
% methods. This method may avoid overfitting problem and handle the caculate
% efficiency problem of nonlinear machine learning method.

% Parameters:
% __________
% pl/pnl: percent of selected highly linear and nonlinear feature.

% attention: This method is highly ideally, need to be attentioned further.


% use pearson coefficent to select top feature.
coeff =  abs(corr(x,y));
% NaN value is at end of sort array.
[~,linear_index] = sort(coeff,'descend','MissingPlacement','last'); 
% number of NaN coefficent, which means sum of the predictor variable is 0.
nonan = find(~isnan(coeff));
% select the top pl feature by pearson coefficent.
select_linear_index = linear_index(1:round(pl*length(nonan)));
% use random forest to select top feature.
[~,~,nonlinear_index] = RF(x(:,nonan),y);
select_nonlinear_index = nonan(nonlinear_index(1:round(pnl*length(nonan))));
% combine the selected feature
select_index = union(select_linear_index,select_nonlinear_index);

end

%%
function index_rank = rank(X,Y)
% give the predict and target variables, remove features in predict
% variables in order, and get the difference of R2 of full and baseline
% regression, then use this diff represents the impact of the removed variable.
% and get the rank of the variable.

% y = a(x,w,e,f...)+b; get corresponding R2; 
% exclude x,w,e,f,etc,respectively; get R2.
% contrast these R2 value to find the impact of each terms.

% default use random forest regression.

NumTrees = 100;
[n_time,n_feature] = size(X);
% remove nan
[X,Y] = check_terms(X,Y);
% get R2 using full variables in X.
R2_full = RF(X,Y);
R2 = zeros(n_feature,1);
% loop for all feature
for i = 1:n_feature
    %
    II = 1:n_feature;
    II(i) = [];
    % 
    B = TreeBagger(NumTrees,...
                   X(:,II),Y,...
                   'method','regression',...
                   'Surrogate','on',...
                   'OOBPredictorImportance','on',...
                   'PredictorSelection','curvature');    
    R2(i) = index(predict(B,X(:,II)),Y,1);  
end
    % indexRank means the importance array of predictors.
    dev_R2 = R2_full - R2;
    [~,index_rank] = sort(dev_R2,'descend'); 
end
