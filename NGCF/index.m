function index = index

    index.R2 = @determination_coefficient;
    index.MSE = @mean_squared_error;
    index.RMSE = @root_mean_squared_error;
    index.MAE = @mean_average_error;
    index.AIC = @

end

function R2 = determination_coefficient(x,y)

R2 = 1 - (sum((y-x).^2)./sum((y-mean(y)).^2));

end

function MSE = mean_squared_error(x,y)

MSE = nansum((y-x).^2)/length(y);

end

function RMSE = root_mean_squared_error(x,y)

RMSE = mean_squared_error(x,y).^0.5;

end

function MAE = mean_average_error(x,y)

MAE = sum(abs(x-y))/numel(y);

end




