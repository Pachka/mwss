#' Generate summaries of trajectories
#' @usage keyoutput(trajmwss,  scale = 0, focus = NULL, outb_Thhold = 1)
#'
#' @description \code{keyoutput} apply to a mwss object returns a list with number of simulation (\code{nSimulations}), number of wards (\code{nwards}),
#' the number of admission over the whole simulation (\code{nadm}),
#' the maximal number of patients that were simultaneously isolated over the whole simulation (\code{ISO}),
#' the maximal number of professionals that simultaneously took sick leave over the whole simulation (\code{SL}),
#' the statistical summary of infections over the whole simulation (\code{infection}),
#' the statistical summary of incidence over the whole simulation (\code{incidence}),
#' the statistical summary of tests implemented over the whole simulation (\code{tests}).
#'
#' Used with the \code{scale = 0}, the function returns statistics calculated at the facility scale.
#' Used with the \code{scale = 1}, the function returns statistics calculated for each ward.
#'
#' Using focus argument limit outputs to one of the three last list contains in the generic method.
#'
#' \code{n_extinction} is the number of extinction among all simulations. \code{perc_outbreak} is the proportion of simulations were an epidemic was observed (only with \code{scale = 1}).
#' \code{n_noso} is the average number of nosocomial infections in a given population over the simulations.
#' \code{n_intro} is the average number of infected patients introduced in a given population over the simulations.
#' \code{n_out} is the average number of professionals of a given population who was infected in community over the simulations.
#' \code{incPA, incPM, incPS, incHA, incHM} and \code{incHS} are respectively the average incidence of asymptomatic, mild symptomatic and severe symptomatic patients and professionals (H: healthcare workers) over the simulations.
#' \code{incP} and \code{incH} are the respective incidence among patients and professionals (including asymptomatics and symptomatics).
#' \code{nTestP} and \code{nTestH} are the respective number of tests performed on patients and professionals.
#' \code{n_dtestP} and \code{n_dtestH} are the respective average daily number of tests performed on patients and professionals.
#'
#' @param trajmwss mwss object. List of list produced by \code{multisim} function
#' @param scale NULL or numeric vector. 0 is the facility scale, 1 is the ward scale.
#' @param focus String. Limit the output to infection data (\code{focus = "infections"}), incidence data (\code{focus = "incidence"}) or test counter  (\code{focus = "test"}).
#' @param outb_Thhold Integer. Defines the minimal number of cases requested to define an epidemic (default is 1, only used with \code{scale = 1}).
#'
#' @importFrom stats quantile
#' @importFrom data.table ':='
#'
#' @return List
#'
#' @examples
#' data("toydata")
#' list2env(toydata,envir=.GlobalEnv)
#' gdata <- build_gdata()
#' model <- mwss(ward_names, pop_size_P, pop_size_H, nVisits, LS, gdata, tSim = 30)
#' results <- multisim(model, 5, ward_names)
#'
#' keyoutput(results,  scale = 0, focus = NULL)
#' keyoutput(results,  scale = 1, focus = NULL, outb_Thhold = 1)
#' keyoutput(results,  scale = 0, focus = "incidence")
#'
#' @export

keyoutput <- function(trajmwss,
           scale = 0,
           focus = NULL,
           outb_Thhold = 1){

    nSimulations <- length(trajmwss)
    nwards <- length(trajmwss[[1]]$node %>% unique)
    ntimepoint <- length(trajmwss[[1]]$time %>% unique)
    iterations <-
      lapply(seq(nSimulations), function(x)
        rep(x, nwards)) %>% unlist
    diterations <-
      lapply(seq(nSimulations), function(x)
        rep(x, nwards * ntimepoint)) %>% unlist

    # total number of admission per simulation
    if (scale == 0)
      nadm <- lapply(trajmwss, function(sim) {
        sim[time == max(time), sum(nadm)]
      }) %>% unlist

    if (scale == 1) {
      nadm <- lapply(trajmwss, function(sim) {
        sim[time == max(time), c("node", "nadm")]
      }) %>% do.call(rbind, .)
      nadm[, iteration := iterations]
    }

    admInf <- lapply(trajmwss, function(sim) {
      sim[time == max(time), c("node", "admInf")]
    })

    # total number of patient nosocomial infection
    n_nosoP <- lapply(trajmwss, function(sim) {
      sim[time == max(time), c("node", "infP")]
    })

    # total number of professional nosocomial infection
    n_nosoH <- lapply(trajmwss, function(sim) {
      sim[time == max(time), c("node", "infH")]
    })

    # total number of professional infection in the community
    n_outH <- lapply(trajmwss, function(sim) {
      sim[time == max(time), c("node", "infHout")]
    })

    # max simultaneous number of patients isolated
    if("ISO" %in% class(trajmwss))
    n_ISO <- trajmwss

    # max simultaneous number of patients isolated
    n_SL <- trajmwss

    ######## summarize infection
    if (scale == 0) {
      n_nosoP <- lapply(n_nosoP, function(sim)
        sim[, sum(infP)]) %>% unlist
      admInf <-
        lapply(admInf, function(sim)
          sim[, sum(admInf)]) %>% unlist
      n_nosoH <-
        lapply(n_nosoH, function(sim)
          sim[, sum(infH)]) %>% unlist
      n_outH <-
        lapply(n_outH, function(sim)
          sim[, sum(infHout)]) %>% unlist

      infection <- list(
        n_extinction = sum(n_nosoP + n_nosoH == 0),
        P = list(
          n_noso = n_nosoP,
          quantiles_noso = n_nosoP %>% quantile,
          n_intro = admInf,
          quantiles_intro = admInf %>% quantile
        ),
        H = list(
          n_noso = n_nosoH,
          quantiles_noso = n_nosoH %>% quantile,
          n_out = n_outH,
          quantiles_out = n_outH %>% quantile,
          totH = n_nosoH + n_outH,
          quantiles_totH = (n_nosoH + n_outH) %>% quantile
        )
      )

      # contact restriction for patients
      if("ISO" %in% class(trajmwss)){
        iso_cpmt <- c(sapply(
          c(sapply(
            paste("PW", c("S", "E", "EA", "ES", "IA", "IM", "IS"), sep = "_"), function(pref)
              paste(pref, c("NI", "LI", "HI"), sep = "_")
          )), function(pref)
            paste(pref, "T", sep = "_")
        )) %>% paste(., collapse = " , ")

      n_ISO %<>% lapply(., function(sim) {
        eval(parse(text=(paste("
        sim[ ,sizeISO := sum(",iso_cpmt,"), by=time]"
                               ))))
        max(sim[, sizeISO])
      }) %>% unlist

      ISO <- list(maxnISO = n_ISO,
                  quantiles = quantile(n_ISO))
      }


      # sick leave
      n_SL %<>% lapply(., function(sim) {
        sim[, nSL := sum(SL + ESL), by = time]
        max(sim[, nSL])
      }) %>% unlist

      SL <- list(maxnSL = n_SL,
                 quantiles = quantile(n_SL))
    }

    if (scale == 1) {
      n_nosoP %<>% do.call(rbind, .) %>% .[, iteration := iterations]
      admInf %<>% do.call(rbind, .) %>% .[, iteration := iterations]
      n_nosoH %<>% do.call(rbind, .) %>% .[, iteration := iterations]
      n_outH %<>% do.call(rbind, .) %>% .[, iteration := iterations]

      n_nosoP[, infH := n_nosoH$infH]
      p_outbreak <-
        n_nosoP[, (infP + infH) >= outb_Thhold, by = c("node", "iteration")][, sum(V1) /
                                                                               nSimulations, by = node]
      setnames(p_outbreak, "V1", "perc_outbreak")

      totH <- n_nosoH
      totH[, infHout := n_outH$infHout]
      totH[, infH := infH + infHout, by = c("node", "iteration")]
      totH[, infHout := NULL]

      n_nosoP[, infH := NULL]

      infection <- list(
        p_outbreak = p_outbreak,
        P = list(
          n_noso = n_nosoP[, c("iteration", "node", "infP")],
          quantiles_noso = n_nosoP[, `:=`(
            `0%` = quantile(infP, p = 0),
            `25%` = quantile(infP, p =
                               0.25),
            `50%` = quantile(infP, p =
                               0.50),
            `75%` = quantile(infP, p =
                               0.75),
            `100%` = quantile(infP, p =
                                1)
          ), by = node][, c("node", "0%", "25%", "50%", "75%", "100%")] %>% unique,
          n_intro = admInf[, c("iteration", "node", "admInf")],
          quantiles_intro = admInf[, `:=`(
            `0%` = quantile(admInf, p = 0),
            `25%` = quantile(admInf, p =
                               0.25),
            `50%` = quantile(admInf, p =
                               0.50),
            `75%` = quantile(admInf, p =
                               0.75),
            `100%` = quantile(admInf, p =
                                1)
          ), by = node][, c("node", "0%", "25%", "50%", "75%", "100%")] %>% unique
        ),
        H = list(
          n_noso = n_nosoH[, c("iteration", "node", "infH")],
          quantiles_noso = n_nosoH[, `:=`(
            `0%` = quantile(infH, p = 0),
            `25%` = quantile(infH, p =
                               0.25),
            `50%` = quantile(infH, p =
                               0.50),
            `75%` = quantile(infH, p =
                               0.75),
            `100%` = quantile(infH, p =
                                1)
          ), by = node][, c("node", "0%", "25%", "50%", "75%", "100%")] %>% unique,
          n_out = n_outH[, c("iteration", "node", "infHout")],
          quantiles_out = n_outH[, `:=`(
            `0%` = quantile(infHout, p = 0),
            `25%` = quantile(infHout, p =
                               0.25),
            `50%` = quantile(infHout, p =
                               0.50),
            `75%` = quantile(infHout, p =
                               0.75),
            `100%` = quantile(infHout, p =
                                1)
          ), by = node][, c("node", "0%", "25%", "50%", "75%", "100%")] %>% unique,
          totH = totH[, c("iteration", "node", "infH")],
          quantiles_totH = totH[, `:=`(
            `0%` = quantile(infH, p = 0),
            `25%` = quantile(infH, p =
                               0.25),
            `50%` = quantile(infH, p =
                               0.50),
            `75%` = quantile(infH, p =
                               0.75),
            `100%` = quantile(infH, p =
                                1)
          ), by = node][, c("node", "0%", "25%", "50%", "75%", "100%")] %>% unique
        )
      )

      # Contact restriction for patients

      if("ISO" %in% class(trajmwss)){
        iso_cpmt <- c(sapply(
          c(sapply(
            paste("PW", c("S", "E", "EA", "ES", "IA", "IM", "IS"), sep = "_"), function(pref)
              paste(pref, c("NI", "LI", "HI"), sep = "_")
          )), function(pref)
            paste(pref, "T", sep = "_")
        )) %>% paste(., collapse = " , ")

        n_ISO %<>% lapply(., function(sim) {
          eval(parse(text=(paste("
        sim[ ,sizeISO := sum(",iso_cpmt,"), by=c('node', 'time')]"
          ))))
          sim[, maxISO := max(sizeISO), by = node] %>% .[, c("node", "maxISO")] %>% unique
        })  %>% do.call(rbind, .)

      n_ISO[, iteration := iterations]

      quantiles_ISO <- n_ISO
      quantiles_ISO %<>%
        .[, `:=`(
          `0%` = quantile(maxISO, p = 0),
          `25%` = quantile(maxISO, p = 0.25),
          `50%` = quantile(maxISO, p = 0.50),
          `75%` = quantile(maxISO, p = 0.75),
          `100%` = quantile(maxISO, p = 1)
        ), by = node] %>% .[, c("node", "0%", "25%", "50%", "75%", "100%")] %>% unique

      ISO <- list(maxnISO = n_ISO,
                  quantiles = quantiles_ISO)
      }

      # sick leave
      n_SL %<>% lapply(., function(sim) {
        sim[, nSL := sum(SL + ESL), by = c("node", "time")]
        sim[, nSL := max(nSL), by = node]
        sim[, c("node", "nSL")] %>% unique
      }) %>% do.call(rbind, .)

      n_SL[, iteration := iterations]

      quantiles_SL <- n_SL %>%
        .[, `:=`(
          `0%` = quantile(nSL, p = 0),
          `25%` = quantile(nSL, p = 0.25),
          `50%` = quantile(nSL, p = 0.50),
          `75%` = quantile(nSL, p = 0.75),
          `100%` = quantile(nSL, p = 1)
        ), by = node] %>% .[, c("node", "0%", "25%", "50%", "75%", "100%")] %>% unique

      SL <- list(maxnSL = n_SL,
                 quantiles = quantiles_SL)
    }

    # incidence
    incidence <- lapply(trajmwss, function(sim) {
      PA <- sim[time == max(time), incPA, by = node]
      PM <- sim[time == max(time), incPM, by = node]
      PS <- sim[time == max(time), incPS, by = node]

      HA <- sim[time == max(time), incHA, by = node]
      HM <- sim[time == max(time), incHM, by = node]
      HS <- sim[time == max(time), incHS, by = node]

      merge(PA, PM) %>%
        merge(., PS) %>%
        merge(., HA) %>%
        merge(., HM) %>%
        merge(., HS)

    }) %>% do.call(rbind, .)

    incidence[, iteration := lapply(1:nSimulations, function(iter)
      rep(iter, nwards)) %>% unlist]

    if (scale == 1)
      incidence <-
      list(
        incidence = incidence[, c("iteration",
                                  "node",
                                  "incPA",
                                  "incPM",
                                  "incPS",
                                  "incHA",
                                  "incHM",
                                  "incHS")],
        quantilesP = incidence[, `:=`(
          `0%` = quantile(incPA + incPM + incPS, p = 0),
          `25%` = quantile(incPA +
                             incPM + incPS, p = 0.25),
          `50%` = quantile(incPA +
                             incPM + incPS, p = 0.50),
          `75%` = quantile(incPA +
                             incPM + incPS, p = 0.75),
          `100%` = quantile(incPA +
                              incPM + incPS, p = 1)
        ), by = node][, c("node", "0%", "25%", "50%", "75%", "100%")] %>% unique,
        quantilesH = incidence[, `:=`(
          `0%` = quantile(incHA + incHM + incHS, p = 0),
          `25%` = quantile(incHA +
                             incHM + incHS, p = 0.25),
          `50%` = quantile(incHA +
                             incHM + incHS, p = 0.50),
          `75%` = quantile(incHA +
                             incHM + incHS, p = 0.75),
          `100%` = quantile(incHA +
                              incHM + incHS, p = 1)
        ), by = node][, c("node", "0%", "25%", "50%", "75%", "100%")] %>% unique
      )

    if (scale == 0) {
      incidence[, `:=`(
        incPA = sum(incPA),
        incPM = sum(incPM),
        incPS = sum(incPS),
        incHA = sum(incHA),
        incHM = sum(incHM),
        incHS = sum(incHS)
      ), by = iterations]
      incidence %<>% .[,-1] %>% unique

      incidence <- list(
        incidence = incidence,
        quantiles_incP = incidence[, incPA + incPM + incPS] %>% quantile,
        quantiles_incH = incidence[, incHA + incHM + incHS] %>% quantile
      )
    }

    # Tests counter

    n_test <- n_dtest <- copy(trajmwss)

    n_test %<>% lapply(., function(sim) {
      sim[time == max(time), c("node", "nTestP", "nTestH")]
    }) %>% do.call(rbind, .)

    n_test[, iteration := iterations]

    # Daily test
    n_dtest %<>% lapply(.,  function(sim) {
      sim[, `:=`(n_dtestP = c(nTestP[1], diff(nTestP)),
                 n_dtestH = c(nTestH[1], diff(nTestH))), by = node]
      sim[, c("node", "time", "n_dtestP", "n_dtestH")]
    }) %>% do.call(rbind, .)

    n_dtest[, iteration := diterations]

    if (scale == 0) {
      # max ntest
      n_test[, ':='(nTestP = sum(nTestP),
                    nTestH = sum(nTestH)), by = iteration]
      n_test %<>% .[, c("iteration", "nTestP", "nTestH")] %>% unique

      # daily ntest
      n_dtest[, ':='(n_dtestP = sum(n_dtestP),
                     n_dtestH = sum(n_dtestH)), by = c("iteration", "time")]
      n_dtest %<>% .[, c("iteration", "time", "n_dtestP", "n_dtestH")] %>% unique

      tests <- list(
        ntest = n_test,
        quantilesP = quantile(n_test[, nTestP]),
        quantilesH = quantile(n_test[, nTestH]),
        ndtest = n_dtest,
        quantilesPdaily = quantile(n_dtest[, n_dtestP]),
        quantilesHdaily = quantile(n_dtest[, n_dtestH])
      )
    }

    if (scale == 1)
      tests <-
      list(
        ntest = n_test[, c("iteration", "node", "nTestP", "nTestH")],
        quantilesP = n_test[, `:=`(
          `0%` = quantile(nTestP , p = 0),
          `25%` = quantile(nTestP , p =
                             0.25),
          `50%` = quantile(nTestP , p =
                             0.50),
          `75%` = quantile(nTestP , p =
                             0.75),
          `100%` = quantile(nTestP , p =
                              1)
        ), by = node][, c("node", "0%", "25%", "50%", "75%", "100%")] %>% unique,
        quantilesH = n_test[, `:=`(
          `0%` = quantile(nTestH , p = 0),
          `25%` = quantile(nTestH , p =
                             0.25),
          `50%` = quantile(nTestH , p =
                             0.50),
          `75%` = quantile(nTestH , p =
                             0.75),
          `100%` = quantile(nTestH , p =
                              1)
        ), by = node][, c("node", "0%", "25%", "50%", "75%", "100%")] %>% unique,
        ndtest = n_dtest[, c("iteration", "node", "time", "n_dtestP", "n_dtestH")],
        quantilesPdaily = n_dtest[, `:=`(
          `0%` = quantile(n_dtestP, p = 0),
          `25%` = quantile(n_dtestP , p =
                             0.25),
          `50%` = quantile(n_dtestP , p =
                             0.50),
          `75%` = quantile(n_dtestP , p =
                             0.75),
          `100%` = quantile(n_dtestP , p =
                              1)
        ), by = node][, c("node", "0%", "25%", "50%", "75%", "100%")] %>% unique,
        quantilesHdaily = n_dtest[, `:=`(
          `0%` = quantile(n_dtestH , p = 0),
          `25%` = quantile(n_dtestH , p =
                             0.25),
          `50%` = quantile(n_dtestH , p =
                             0.50),
          `75%` = quantile(n_dtestH , p =
                             0.75),
          `100%` = quantile(n_dtestH , p =
                              1)
        ), by = node][, c("node", "0%", "25%", "50%", "75%", "100%")] %>% unique
      )


    if (is.null(focus)){
      if("ISO" %in% class(trajmwss))
      x <-  list(
        nSimulations = nSimulations,
        nwards = nwards,
        nadm = nadm,
        ISO = ISO,
        SL = SL,
        infection = infection,
        incidence = incidence,
        tests = tests
      ) else
        x <-  list(
          nSimulations = nSimulations,
          nwards = nwards,
          nadm = nadm,
          SL = SL,
          infection = infection,
          incidence = incidence,
          tests = tests
        )
        }
    else
      if (focus == "infections")
        x <-  infection
    else
      if (focus == "incidence")
        x <- incidence
    else
      if (focus == "test")
        x <- tests
    else
      warning(
        paste(
          "keyoutput does not know how to handle focus: ",
          focus,
          "and can only be used with focus NULL, \"infections\", \"incidence\", and \"test\""
        )
      )

    return(x)
  }
