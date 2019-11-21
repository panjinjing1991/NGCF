function indicator = index(x,y,type)

if type==1
    indicator = determination_coefficient(x,y);
elseif type==2
    indicator = mean_squared_error(x,y);
end
end

function R2 = determination_coefficient(x,y)

R2 = 1 - (sum((y-x).^2)./sum((y-mean(y)).^2));

end

function MSE = mean_squared_error(x,y)

MSE = nansum((y-x).^2)/length(y);

end
