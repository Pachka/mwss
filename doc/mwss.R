## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE------------------------------------------------------------
library(VTMstocMod)

## -----------------------------------------------------------------------------
params=list(
  betaPP     = 0.35,
  betaHH     = 0.35,
  betaPH     = 0.35,
  betaHP     = 0.6,
  betaPP_AL  = 0.35,
  betaHH_AL  = 0.35,
  betaPH_AL  = 0.35,
  betaHP_AL  = 0.6,
  gamma1 = 1/5,
  gamma2 = 1/2,
  gamma3 = 1/7,
  gamma4 = 1/8,
  gamma5 = 1/10,
  gammaT = 1/2,
  gammaBSL= 1/1,
  gammaSL= 1/14,
  gammapreIso= 1/1,
  gammaIso= 1/11,
  timeExCli=1,
  timeTest=1,
  nHCWS_AL=1,
  p=70/100,
  p2=20/100,
  p2p=50/100,
  pDIC=50/100,
  pD=10/100,
  pT=20/100,
  pSL=80/100,
  PrevCom=10/100,
  sensInf=90/100,
  sensNInd=30/100
)

## ---- fig.width=7, fig.height=6-----------------------------------------------
  xstart <- startvec(wards.names.ex, pop.size_pat.ex, pop.size_hcws.ex, turnover.ex, introduction.ex)
  plot_wardsConnect(matContact.ex, pop.size_pat.ex)

## -----------------------------------------------------------------------------
  nSimulations  <- 400
  duration      <- 40
  isoWard <- TRUE
  airlockEffective      <- FALSE
  scenario              <- FALSE
  test_str              <- "none"


## ----eval = FALSE-------------------------------------------------------------
#    xstart <- startvec(wards.names.ex, pop.size_pat.ex, pop.size_hcws.ex, turnover.ex)

## ----eval = FALSE-------------------------------------------------------------
#    output <- data.generate.MW(params.ex, xstart, matContact.ex, nSimulations, duration, isoWard, airlockEffective, scenario, test_str)

## ----eval = FALSE-------------------------------------------------------------
#    output.res <- resume.output(output,test_str)
#    output.res$nSimulations
#    output.res$nWards
#    output.res$TrSec
#    output.res$incP
#    output.res$incH
#    output.res$nTestsHosp
#    output.res$nDailyTestsHosp
#
#    plotDailyTest(output.res$nDailyTestsHosp$perSim)

## ----eval = FALSE-------------------------------------------------------------
#    output.resW <- resume.output.ward(output)
#    output.resW$TrSec
#    output.resW$IncP
#    output.resW$IncH

## ----eval = FALSE-------------------------------------------------------------
#    output.res <- resume.tr(output)
#    output.res <- resume.inc(output)

## ----eval = FALSE-------------------------------------------------------------
#    params$PrevCom <- 0

## ----eval = FALSE-------------------------------------------------------------
#    introduction <- data.frame(wards.names =  "Ward2",
#                               number = 1,
#                               type="P",
#                               epiStat = "E")
#    xstart <- startvec(wards.names, pop.size_pat, pop.size_hcws, turnover, introduction)
#    xstart
#    params
#
#    nSimulations  <- 400
#    duration      <- 30 #days
#    isoWard <- TRUE
#    airlockEffective      <- TRUE
#
#    output <- data.generate.MW(params, xstart, matContact, nSimulations, duration, isoWard, airlockEffective)

## ----eval = FALSE-------------------------------------------------------------
#    transmission.res <- resume.tr(output)
#
#    transmission.res$nExt
#    transmission.res$perWard$nExt
#
#    plot_pOutbreak(output, matContact, pop.size_pat)

