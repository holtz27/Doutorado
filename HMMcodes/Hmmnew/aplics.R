library(holtz)

theta0=list(mu0=log(var(ytrain)), phi0=0.98, sigma0=0.1)
svn_fit=holtz::svnHMM(y=ytrain-mean(ytrain), 
                      m=200, gmax=5, theta_init=theta0)

theta0=list(mu0=log(var(ytrain)), phi0=0.98, sigma0=0.1, nu0=23)
svt_fit=holtz::svtHMM(y=ytrain-mean(ytrain), 
                      m=200, gmax=5, theta_init=theta0)

theta0=list(mu0=log(var(ytrain)), phi0=0.98, sigma0=0.1, nu0=5)
svs_fit=holtz::svsHMM(y=ytrain-mean(ytrain), 
                      m=200, gmax=5, theta_init=theta0)

theta0=list(mu0=log(var(ytrain)), phi0=0.98, sigma0=0.1, nu0=21)
svvg_fit=holtz::svvgHMM(y=ytrain-mean(ytrain), 
                        m=200, gmax=5, theta_init=theta0)

summary(svn_fit)
summary(svt_fit)
summary(svs_fit)
summary(svvg_fit)

holtz::info_crit(svn_fit)
holtz::info_crit(svt_fit)
holtz::info_crit(svs_fit)
holtz::info_crit(svvg_fit)

################################################################################
################################################################################

svmn_fit=holtz::svmnHMM(y=ytrain, m=200, gmax=5)

theta0=list(mu0=log(var(ytrain)), phi0=0.98, sigma0=0.1,
            beta0=c(mean(ytrain), -0.01, -0.05), nu0=23)
svmt_fit=holtz::svmtHMM(y=ytrain, m=200, gmax=5, theta_init=theta0)

theta0=list(mu0=log(var(ytrain)), phi0=0.98, sigma0=0.1,
            beta0=c(mean(ytrain), -0.01, -0.05), nu0=5)
svms_fit=holtz::svmsHMM(y=ytrain, m=200, gmax=5, theta_init=theta0)

theta0=list(mu0=log(var(ytrain)), phi0=0.98, sigma0=0.1,
         beta0=c(mean(ytrain), -0.01, -0.05), nu0=21)
svmvg_fit=holtz::svmvgHMM(y=ytrain, m=200, gmax=5, theta_init=theta0)



summary(svmn_fit)
summary(svmt_fit)
summary(svms_fit)
summary(svmvg_fit)

holtz::info_crit(svmn_fit)
holtz::info_crit(svmt_fit)
holtz::info_crit(svms_fit)
holtz::info_crit(svmvg_fit)

svmn_log.vol=holtz::logvol(svmn_fit, plot=TRUE)
svmt_log.vol=holtz::logvol(svmt_fit, plot=TRUE)
svms_log.vol=holtz::logvol(svms_fit, plot=TRUE)
svmvg_log.vol=holtz::logvol(svmvg_fit, plot=TRUE)

svmn_fpr=holtz::svmnHMM.pseudo_residuals(svmn_fit, ytest)
svmt_fpr=holtz::svmtHMM.pseudo_residuals(svmt_fit, ytest)
svms_fpr=holtz::svmsHMM.pseudo_residuals(svms_fit, ytest)
svmvg_fpr=holtz::svmvgHMM.pseudo_residuals(svmvg_fit, ytest)

#save(svn_fit,svt_fit,svs_fit,svvg_fit,
#     svmn_fit,svmt_fit,svms_fit,svmvg_fit,
#     file='~/Documentos/svmHMM/sp500/sp500_HMM.RData')

#save(svn_fit,svt_fit,svs_fit,svvg_fit,
#     svmn_fit,svmt_fit,svms_fit,svmvg_fit,
#     file='~/Documentos/svmHMM/nikkei/nikkei_HMM.RData')

#save(svn_fit,svt_fit,svs_fit,svvg_fit,
#     svmn_fit,svmt_fit,svms_fit,svmvg_fit,
#     file='~/Documentos/svmHMM/dax/dax_HMM.RData')

#save(svn_fit,svt_fit,svs_fit,svvg_fit,
#     svmn_fit,svmt_fit,svms_fit,svmvg_fit,
#     file='~/Documentos/svmHMM/ibovespa/ibovespa_HMM.RData')



