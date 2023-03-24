# Miscellaneous
Miscellaneous helper functions I use to rapidly perform basic epidemiologic analyses

# Resmat_lm
This is a complicated helper function I've piecemeal added to throughout the years to handle all sorts of repeated linear regression analyses.
It outputs R and html (kable) tables of key regression output, a list of the model objects, and a coefficient forest plot. Supported regression types
include linear regression (lm) and binomial GLMs (glm with binomial family and link "logit" for logistic regressions or Firth's logistic regression using logistf, "log" for binomial regressions, or "probit" for probit regressions).

Inputs:
- prednames: A character vector of the column names of predictors of interest you'd like to evaluate in separate regressions. These column names can refer to
  character, factor, or numeric columns in the data frame.
- outnames: A character vector of the column names of outcome variables you'd like to evaluate in separate regressions. These column names can refer to 
  character, factor, or numeric columns in the data frame. For each outnames, lm() linear regressions will be run if the referenced column is numeric, or a
  glm() GLM with a binomial() family will be run with the link function specified in the input "binomlink". Does not yet support characters or factors with >2 unique values (i.e., it can only accomodate continous or binary outcome variables).
- covnames: Either a character vector of column names for all the covariates you'd like to include in each regression or a list of character vectors equal in length
  to the length of the unique combination of prednames and outnames. Defaults to NULL, meaning no covariates are included. This can be a list in order to control for
  separate covariates for each unique combination of prednames and outnames. These can include 'ns' natural spline terms as defined by the 'splines' package.
- Data: A data frame object of class "data.frame".
- logout: If TRUE, this will log transform each outcome variable (log base will be "logbaseout") prior to running the regression and then return both the raw
  coefficients as well as converted coefficients to percent changes using the formula percent change = (exp(raw coefficient)-1) x 100. Can be either a single
  value or a logical vector of length equal to the length of outnames. Defaults to FALSE.
- logpred: If TRUE, this will log transform each predictor variable (log base will be "logbasepred") prior to running the regression and then return the raw
  coefficients. Can be either a single value or a logical vector of length equal to the length of prednames. Defaults to FALSE.
- logbasepred: The base for the log transformation of one or more of the predictor variables. Defaults to 10.
- logbaseout: The base for the log transformation of the outcome variables. Defaults to exp(1) (i.e., natural log transformation).
- Outtitle: Defines the title used for the column of outcome variables in the table and on the plot. Defaults to "Outcome".
- Predtitle: Defines the title used for the column of predictor variables in the table and on the plot. Defaults to "Exposure".
- MyMult: The table will output raw and transformed coefficients and CIs. This multiplier is an optional scalar to multiply the raw coefficients for an estimate
  of the coefficient for something other than a 1-unit increase. Defaults to 1, and can be either a number or the character "IQR", in which case the
  raw coefficient will be multiplied by the IQR of the predictor variable.
- ixterm: A column name for an optional interaction term for the prednames predictors of interest (I abbreviate interaction as 'ix'). Can refer to either a
  character, factor, or numeric-class column in the Data data frame. All relevant interaction and main coefficients will be output. Defaults to NULL.
- Firth: If TRUE and if a given outnames variable is binary, this will perform a Firth's logistic regression (logistf) instead of a logistic regression (glm).
  Defaults to FALSE.
- binomlink: The link function for the binomial family for GLMs of binary outnames. Can be "logit", "log", or "probit". Defaults to "logit".
- robust: If TRUE, this calculates and returns robust CIs and p-values based on the vcovHC function instead of the default CIs and p-values. Defaults to TRUE.
- HCtype: This is the formula used for the vcovHC function calculation of robust CIs and p-values. Defaults to "HC0".
- predspline: If TRUE, this transforms each predictor of interest into a natural spline term using ns from the splines package. A coefficient for each degree of
  freedom is returned. Defaults to FALSE.
- predsplinedf: This specifies the number of degrees of freedom for the predictor variable natural splines if predspline is TRUE. Defaults to 3.
- plotOR: If TRUE, this will plot ORs (i.e., it will plot exp(coefficient) values) in the coefficient forest plot. Defaults to FALSE.
- plotPercChange: If TRUE, this will plot percent changes (i.e., it will plot (exp(coefficient)-1) x 100 values) in the coefficient forest plot.
  Defaults to FALSE. If logistic regressions were run and both plotOR and plotPercChange are FALSE, raw coefficients will be plotted.
- LOOCV: If TRUE, this performs leave one out cross validation using the 'caret' package. If a given model is a linear regression (lm), the RMSE will be 
  reported in the tables. If a model is a GLM (glm), the accuracy will be reported in the tables. The train function input method is
  trainControl(method="repeatedcv", number=10, repeats=50). Defaults to FALSE. 
- facetcol: This defines the number of facet_wrap columns for the coefficient plots (facetted by outcome variable). If NULL, this will equal the number of
  unique outnames values. Defaults to NULL.
- covix: This character vector of covariate column names allows for covariates to also have interaction terms with the ixterm interaction variable. If one 
  wishes to have all covariates have interaction terms with ixterm, the character vector for covix should be the same as that for covnames. Defaults to 
  NULL, meaning that there are no interactions with any of the covariates.
- ixpred: If TRUE, an interaction term is included between each of the prednames predictors of interest and ixterm. This can be set to FALSE to allow for interactions
  only with covariates as defined in covix. If covix is not NULL and ixpred is TRUE, there will be interaction terms with both the predictor of interest and the 
  covariates. Defaults to TRUE.
- extradiag: If TRUE, this will include extra diagnostic information in the tables. For linear regressions (lm), this will include a p-value for the heteroskedasticity
  of model residuals (lmtest::bptest), the number of Cook's distances > 0.5, and the maximum Cook's distance, the latter two of which define high leverage points.
  For GLMs (glm), this will include output from DHARMa package tests, including a heteroskedasticity of residuals p-value from DHARMa::testQuantiles, the 
  minimum p-value from DHARMa::testQuantiles, the number of Cook's distances > 0.5, the maximum Cook's distance, a p-value for a test of uniformity 
  (DHARMa::testUniformity), a p-value for outliers (DHARMa::testOutliers), and a p-value for dispersion (DHARMa::testDispersion). Defaults to TRUE.
- leverage.test: If TRUE, this runs all regressions while excluding "high leverage" observations as defined by the user. Defaults to FALSE.
- leverage.cutoff: This is the Cook's distance cutoff above which observations will be excluded if leverage.test is TRUE. Defaults to 0.2.
- leverage.meancutoff: If leverage.cutoff is set to NULL, this is used as a multiplier by the mean Cook's distance to get a new cutoff above which observations
  will be excluded if leverage.test is TRUE. For example, if set to 2, the cutoff will be set to the mean Cook's distance times 2. If both leverage.cutoff
  and leverage.meancutoff are set to NULL, the cutoff will be set to the mean Cook's distance times 4. Defaults to NULL.
- post.power: If TRUE, this repeatedly simulates data with similar structure to the actual data using SimMultiCorrData::rcorrvar and then uses a user-specified
  predictor of interest coefficient to calculate a simulated outcome variable (all covariate coefficients are drawn from the original model estimates), 
  runs regressions on these simulated data frames, collects the iterated predictor of interest coefficient p-values, and calculates how many of those p-values 
  are <0.05. This is a form of posterior power test for different effect sizes of the predictor of interest association with the outcome variable. Defaults to 
  FALSE.
- effect.size: This is the user-specified effect size for which the post-hoc power is tested if post.power is set to TRUE. Defaults to 0.5.
- nsim: This is the number of simulations for the post-hoc power test if post.power is set to TRUE. Defaults to 1E3.
- mice: If TRUE, this performs multiple imputation (MICE using mice::mice) on the user-specified variables and pools model estimates across the imputed datasets.
  Defaults to FALSE.
- micevars: A character vector of the names of variables to be imputed using MICE. Defaults to NULL.
- miceiter: The number of multiple imputations to perform if mice is set to TRUE. Defaults to 10.

Outputs:
- Kable: An html kable version of the table in the Matrix output.
- Matrix: A table of all the key variables. Beta, LCI, and UCI refer to the coefficient for the predictor of interest and the lower and upper 95% confidence intervals, 
  respectively. Beta.Mult, etc. are those estimates transformed in some way as specified by the user, for instance by being multipled by the IQR, transformed into 
  odds ratios or percent changes, or some other multiplication specified by the user.
- GGplot: A forest ggplot of all predictor of interest coefficients and CIs. 
- LMlist: A list of all the lm, glm, or logistf objects run for each unique combination of predictor and outcome variable names (i.e., prednames and outnames).

# InteractionCoefPlot
This produces a plot of how the coefficient between one variable in a bivariate continuous interaction and the outcome changes over levels of the other continuous variable. This is a useful visualization of bivariate interactions between two continuous variables. 

Inputs:
- model: model of class lm or glm
- data: data frame used to generate model
- pred: column name of the predictor of interest
- ixterm: column name of the interacting variable
- multiplier: multiplier for the coefficients of interest (e.g., an IQR of the predictor of interest), defaults to 1
- coeftransform: optional transformation function for the coefficients (e.g., tenfoldperc if the predictor of interest is on a log10 scale where `tenfoldperc<-function(x) ((10^x)-1)*100)`
- predname: optional new name for pred
- ixname: optional new name for ixterm
- outname: optional new name for the outcome variable
- title: optional plot title
- fillcolor: color for the shading of the 95% CI region
- autotitle: option to automatically generate a title, defaults to FALSE
- addpvallab: option to add a label showing the interaction term p-value, defaults to FALSE
- labsize: label text size if addpvallab is TRUE
- lengthout: number of levels of ixterm generated with the length.out argument for the seq() function, defaults to 50
- shadebysig: option to shade the CI region areas that don't overlap the null a darker tint of fillcolor, defaults to FALSE
- otherix: option to generate both plots, one in which pred and ixterm are flipped; defaults to FALSE
- robust: option to use HC0 robust sandwich errors for the CIs, defaults to TRUE; should be set to 
-  FALSE if the model is a logistic glm or TRUE if model is a robust Poisson model for RR estimation
- logistic: option for if model is a logistic glm, defaults to FALSE; this option sets the null line to 1, exponentiates all the coefficients to bring them to an OR scale, and changes the y label to "OR" instead of "Coefficient"

Outputs:
- MainGGplot: the interaction plot
- OtherGGplot: the other interaction plot if otherix=TRUE
- MainMatrix: the matrix of values used to generate the main plot
- OtherMatrix: the matrix of values used to generate the other plot if otherix=TRUE

Plot axis labels, fonts, etc. can be changed with standard ggplot options. For example, if the model is a robust Poisson GLM, you can set logistic=TRUE when running this function and saving it into an example object plotlist and then change the "OR" label to "RR" by running `plotlist$MainGGplot + ylab("Main Variable RR")` or something similar.

Example:
```
set.seed(1)
exdat<-as.data.frame(mvtnorm::rmvnorm(500,mean=rep(0,5)))
names(exdat)<-paste0("X",1:ncol(exdat))
exdat$X1X2<-exdat$X1*exdat$X2
exdat$yhat<- c(-1 + as.matrix(exdat) %*% c(0.5,1.5,0,-2,2,1))
set.seed(2)
exdat$y<-exdat$yhat + rnorm(500,0,5)
exlm<-lm(y~X1*X2+X3+X4+X5,exdat)
explot<-InteractionCoefPlot(exlm,exdat,"X1","X2",addpvallab=T,shadebysig=T,robust=F,otherix=T)
cowplot::plot_grid(explot$MainGGplot, explot$OtherGGplot, ncol=1)
```
