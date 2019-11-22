function indicator = index(x,y,type)

if type==1
    indicator = determination_coefficient(x,y);
elseif type==2
    indicator = mean_squared_error(x,y);
elseif type==3
    indicator = root_mean_squared_error(x,y);
elseif type==4
    indicator = mean_average_error(x,y);
end

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




