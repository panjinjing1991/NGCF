function pass = causality_test(X,Y)
% use random forest to apply nonlinear granger causality test.

const = ones(numel(Y),1);
R2_full = models([X,const],Y,'origin',3);
R2_baseline = models(const,Y,'origin',3);

if R2_baseline < R2_full
    pass = 1;
else
    pass = 0;
end

end
