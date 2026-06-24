.First.lib <- function(lib,pkg)
{
  # Load the dll file #
  library.dynam("dbpmr",pkg,lib)
  cat("dbpmr 0.1-1 loaded\n")

  # Declare the class structures to be used #
  setClass("run.params",representation(filename="character",
                                       no_pelagic="integer", no_benthic="integer",
                                       spatial_dim="integer",
                                       coupled_flag="logical", diff_method="integer"),
                        prototype(filename="SizeSpectrum",
                                  no_pelagic=as.integer(1), no_benthic=as.integer(0),
                                  spatial_dim=as.integer(0),
                                  coupled_flag=T, diff_method=as.integer(1)))

  setClass("grid.params",representation(mmin="numeric", mmax="numeric", mstep="numeric", moutstep="numeric",
                                        t1="numeric", tmax="numeric", tstep="numeric", toutmin="numeric", toutmax="numeric", toutstep="numeric",
                                        xmin="numeric", xmax="numeric", xstep="numeric", xoutstep="numeric",
                                        ymin="numeric", ymax="numeric", ystep="numeric", youtstep="numeric"),
                        prototype(mmin=-28, mmax=14, mstep=0.2, moutstep=1,
                                  t1=0, tmax=1, tstep=(1/365), toutmin=0, toutmax=1, toutstep=(5/365),
                                  xmin=0, xmax=1000, xstep=50, xoutstep=50,
                                  ymin=0, ymax=1000, ystep=50, youtstep=50))

  setClass("plankton.params",representation(filename="character", speciestype="character",
                                            mmin="numeric", mmax="numeric",
                                            u_0="numeric", lambda="numeric", mu_0="numeric", beta="numeric",
                                            initial_flag="logical", ts_flag="logical"),
                        prototype(filename="plankton", speciestype="plankton",
                                  mmin=-28, mmax=-14,
                                  u_0=0.01, lambda=-1, mu_0=0.2, beta=-0.25,
                                  initial_flag=F, ts_flag=F))
  
  setClass("pelagic.params",representation(filename="character", speciestype="character",
                                           mmin="numeric", mmat="numeric", mmax="numeric",
                                           A="numeric", alpha="numeric", mu_0="numeric", beta="numeric", mu_s="numeric", epsilon="numeric", u_0="numeric", lambda="numeric",
                                           K_pla="numeric", R_pla="numeric", Ex_pla="numeric", K_pel="numeric", R_pel="numeric", Ex_pel="numeric", K_ben="numeric", R_ben="numeric", Ex_ben="numeric",
                                           pref_pla="numeric", pref_pel="numeric", pref_ben="numeric",
                                           q_0="numeric", sig="numeric", trunc="numeric",
                                           prey="numeric", pred="numeric", comp="numeric",gamma_prey="numeric", gamma_pred="numeric", gamma_comp="numeric",
                                           rep_method="integer", initial_flag="logical", ts_flag="logical", fishing_flag="logical"),
                        prototype(filename="pelagic",speciestype="pelagic",
                                  mmin=-14, mmat=7, mmax=14,
                                  A=640, alpha=0.82, mu_0=0.2, beta=-0.25, mu_s=0.1, epsilon=0.1, u_0=0.01, lambda=-1,
                                  K_pla=0.2, R_pla=0.2, Ex_pla=0.3, K_pel=0.2, R_pel=0.2, Ex_pel=0.3, K_ben=0.1, R_ben=0.2, Ex_ben=0.4,
                                  pref_pla=1, pref_pel=1, pref_ben=1,
                                  q_0=2*log(10), sig=1*log(10),
                                  prey=0, pred=0, comp=0.1, gamma_prey=0.33, gamma_pred=0.33, gamma_comp=0.75,
                                  rep_method=as.integer(1), initial_flag=F, ts_flag=F, fishing_flag=F))

  setClass("benthic.params",representation(filename="character",speciestype="character",
                                           mmin="numeric", mmat="numeric", mmax="numeric",
                                           A="numeric", alpha="numeric", mu_0="numeric", beta="numeric", mu_s="numeric", epsilon="numeric", u_0="numeric", lambda="numeric",
                                           K_det="numeric", R_det="numeric", Ex_det="numeric",
                                           pref_det="numeric",
                                           rep_method="integer", initial_flag="logical", ts_flag="logical", fishing_flag="logical"),
                        prototype(filename="benthic",speciestype="benthic",
                                  mmin=-17, mmat=-16, mmax=9,
                                  A=64, alpha=0.75, mu_0=0.2, beta=-0.25, mu_s=0.1, epsilon=0.1, u_0=0.01, lambda=-0.75,
                                  K_det=0.2, R_det=0.2, Ex_det=0.2,
                                  pref_det=1,
                                  rep_method=as.integer(1), initial_flag=F, ts_flag=F, fishing_flag=F))

  setClass("detritus.params",representation(filename="character",speciestype="character",
                                            w_0="numeric",
                                            initial_flag="logical", ts_flag="logical"),
                        prototype(filename="detritus",speciestype="detritus",
                                  w_0=0.6,
                                  initial_flag=F, ts_flag=F))

  setClass("plankton.results",representation(uvals="data.frame",finaluvals="data.frame",biomass="numeric",run="run.params",grid="grid.params",species="plankton.params",trange="numeric",mrange="numeric",xrange="numeric",yrange="numeric"))
  setClass("pelagic.results",representation(uvals="data.frame",finaluvals="data.frame",growth="data.frame",mortality="data.frame",predation="data.frame",fishing="data.frame",biomass="numeric",plabio="numeric",pelbio="numeric",benbio="numeric",predbio="numeric",fishbio="numeric",reproduction="numeric",run="run.params",grid="grid.params",species="pelagic.params",trange="numeric",mrange="numeric",xrange="numeric",yrange="numeric"))
  setClass("benthic.results",representation(uvals="data.frame",finaluvals="data.frame",growth="data.frame",mortality="data.frame",predation="data.frame",fishing="data.frame",biomass="numeric",detbio="numeric",predbio="numeric",fishbio="numeric",reproduction="numeric",run="run.params",grid="grid.params",species="benthic.params",trange="numeric",mrange="numeric",xrange="numeric",yrange="numeric"))
  setClass("detritus.results",representation(biomass="numeric",detin="numeric",detout="numeric",run="run.params",grid="grid.params",species="detritus.params",trange="numeric",mrange="numeric",xrange="numeric",yrange="numeric"))

  setClass("singlefile.results",representation(data="data.frame",filetype="character",speciestype="character",run="run.params",grid="grid.params",trange="numeric",mrange="numeric",xrange="numeric",yrange="numeric"))


  setClass("timestep.data",representation(data="data.frame",spatial_dim="integer",trange="numeric",mrange="numeric",yrange="numeric",xrange="numeric"))

  #Define the method to be used to show these class structure"
  setMethod("show","run.params", function(object) {
            print(t(format(data.frame(filename=object@filename,
                                      no_pelagic=object@no_pelagic, no_benthic=object@no_benthic,
                                      spatial_dim=object@spatial_dim,
                                      coupled_flag=object@coupled_flag, diff_method=object@diff_method,
                                      row.names="Value"))),quote=F,print.gap=4)
            }
  )
                                      
  setMethod("show","grid.params",function(object) {
            print(t(format(data.frame(mmin=object@mmin, mmax=object@mmax, mstep=object@mstep, moutstep=object@moutstep,
                                      t1=object@t1, tmax=object@tmax, tstep=object@tstep, toutmin=object@toutmin, toutmax=object@toutmax, toutstep=object@toutstep,
                                      xmin=object@xmin, xmax=object@xmax, xstep=object@xstep, xout=object@xoutstep,
                                      ymin=object@ymin, ymax=object@ymax, ystep=object@ystep, yout=object@youtstep,
                                      row.names="Value"))),quote=F,print.gap=4)
            }
  )
  
  setMethod("show","plankton.params", function(object) {
            print(t(format(data.frame(filename=object@filename, speciestype=object@speciestype,
                                      mmin=object@mmin, mmax=object@mmax,
                                      mu_0=object@mu_0, beta=object@beta, u_0=object@u_0, lambda=object@lambda,
                                      initial_flag=object@initial_flag, ts_flag=object@ts_flag,
                                      row.names="Value"))),quote=F,print.gap=4)
            }
  )
  
  setMethod("show","pelagic.params", function(object) {
            print(t(format(data.frame(filename=object@filename, speciestype=object@speciestype,
                                      mmin=object@mmin, mmat=object@mmat, mmax=object@mmax,
                                      A=object@A, alpha=object@alpha, mu_0=object@mu_0, beta=object@beta, mu_s=object@mu_s, epsilon=object@epsilon, u_0=object@u_0, lambda=object@lambda,
                                      K_pla=object@K_pla, R_pla=object@R_pla, Ex_pla=object@Ex_pla, K_pel=object@K_pel, R_pel=object@R_pel, Ex_pel=object@Ex_pel, K_ben=object@K_ben, R_ben=object@R_ben, Ex_ben=object@Ex_ben,
                                      pref_pla=object@pref_pla, pref_pel=object@pref_pel, pref_ben=object@pref_ben,
                                      q_0=object@q_0, sig=object@sig, trunc=object@trunc,
                                      prey=object@prey, pred=object@pred, comp=object@comp, gamma_prey=object@gamma_prey, gamma_pred=object@gamma_pred, gamma_comp=object@gamma_comp,
                                      rep_method=object@rep_method, initial_flag=object@initial_flag, ts_flag=object@ts_flag, fishing_flag=object@fishing_flag,
                                      row.names="Value"))),quote=F,print.gap=4)
            }
  )
  
  setMethod("show","benthic.params", function(object) {
            print(t(format(data.frame(filename=object@filename, speciestype=object@speciestype,
                                      mmin=object@mmin, mmat=object@mmat, mmax=object@mmax,
                                      A=object@A, alpha=object@alpha, mu_0=object@mu_0, beta=object@beta, mu_s=object@mu_s, epsilon=object@epsilon, u_0=object@u_0, lambda=object@lambda,
                                      K_det=object@K_det, R_det=object@R_det, Ex_det=object@Ex_det,
                                      pref_det=object@pref_det,
                                      rep_method=object@rep_method, initial_flag=object@initial_flag, ts_flag=object@ts_flag, fishing_flag=object@fishing_flag,
                                      row.names="Value"))),quote=F,print.gap=4)
            }
  )



  setMethod("show","detritus.params", function(object) {
            print(t(format(data.frame(filename=object@filename, speciestype=object@speciestype,
                                      w_0=object@w_0,
                                      initial_flag=object@initial_flag, ts_flag=object@ts_flag,
                                      row.names="Value"))),quote=F,print.gap=4)
            }
  )

}