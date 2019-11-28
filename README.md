# NGCF
Nonlinear Granger causality framework Version 1.0

## General Description

NGCF is a causal time series MATLAB package. It's a method to identify the causal relationship between two variables (e.g.,surface soil moisture and precipitation.) 

## Difference of linear and nonlinear Granger causality

run figure1.m
![figure1](https://github.com/leelew/NGCF/blob/master/Figure/Figure1.png)

## Procedure

![Procedure](https://github.com/leelew/NGCF/blob/master/Figure/Figure2.png)

## Code structure

main.m: This main.m file is apply nonlinear granger causality framework(NGCF) on identify surface soil moisture-precipitation feedback, including both sign and pattern distribution.

1. Use get_terms.m to get independent terms, as period terms(i.e., interannual cycle, seasonal cycle); lagged terms(i.e., lagged P and press); spatial terms(i.e., lagged P and press over selected square indicated spatial impact)

2. Run random forest for P occurrence and surface soil moisture as dependent terms, and independent terms metioned in 1., Addtionally, by using hybrid selection feature method to find the 'best' regression(to avoid overfitting in some content). After applying these regression, get the residual value of both surface soil moisture and precipitation.

3. Then, fit residual of S to residual of POCC, using only S with no P on the previous day.

4. Use block bootstrap to eliminate endogeneity bias and determine significance of the S coefficient in the regression.(Need improved)

5. Calculate S-POCC impacts by dividing the unrestricted model by the restricted model, plotted against the seasonal S anomaly, and taking the mean above and below the median of seasonal anomaly.

## User Agreement

By downloading NGCF you agree with the following points: NGCF is provided without any warranty or conditions of any kind. We assume no responsibility for errors or omissions in the results and interpretations following from application of NGCF.

You commit to cite above papers in your reports or publications.

## License

MIT License

Copyright (c) 2019 Li Lu

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.




