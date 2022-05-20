#' Plot the cumulative incidence
#'
#' @description The function \code{plot_incidence} returns a barplot representing the average cumulative incidence over simulations at different scales and for different populations
#'
#' @usage plot_incidence(
#'           trajmwss, scale = 0,
#'           pop = FALSE, iter = FALSE,
#'           ward = FALSE, display_sd = TRUE)
#'
#' @param trajmwss List of data.table. Epidemiological trajectories simulated by the function \code{mwss::multisim}
#' @param scale Integer. Indicated the scale for statistics: 0/ whole facility, 1/ wards.
#' @param pop String. Define a subpopulation for statistics ("P":Cumulative incidence among patients, "Psymp": Cumulative incidence among patients (only considering symptomatics), "H": Cumulative incidence among professionals, "Hsymp": Cumulative incidence among professionals (only considering symptomatics); default is FALSE: considering all individuals).
#' @param iter Integer. Define a specific simulation. Default is FALSE, providing average values over all simulations.
#' @param ward String. Define a specific ward. Default is FALSE, considering all wards. At the facility scale '0' no ward can be specified.
#' @param display_sd Logical. Determine standard deviation display. Default is TRUE.
#'
#' @importFrom data.table setnames
#' @import ggplot2
#' @importFrom stats sd
#' @importFrom data.table ':='
#'
#' @return Cumulative incidence plot.
#'
#' @examples
#'
#' data("toydata")
#' list2env(toydata,envir=.GlobalEnv)
#' gdata <- build_gdata()
#'
#' model <- mwss(ward_names, pop_size_P, pop_size_H, nVisits, LS, gdata, tSim = 30)
#' results <- multisim(model, 5, ward_names)
#'
#' plot_incidence(trajmwss = results)
#'
#' @export

plot_incidence <- function(trajmwss,
                           scale = 0,
                           pop = FALSE,
                           iter = FALSE,
                           ward = FALSE,
                           display_sd = TRUE) {
  # checks
  if (!isFALSE(ward) & scale == 0) {
    scale = 1
  }

  if (!isFALSE(ward) & !ward %in% unique(trajmwss[[1]]$node))
    stop(paste(
      "If not FALSE, ward argument must be a string among:",
      paste(unique(trajmwss[[1]]$node), collapse = ", ")
    ))

  if (!isFALSE(iter) &
      (!is.numeric(iter) | !(iter %in% seq(length(trajmwss)))))
    stop(paste(
      "If not FALSE, iter argument must be a number between 1 and",
      length(trajmwss)
    ))

    if (is.numeric(iter))
      trajmwss %<>% .[iter]

  if (scale == 1) {
    if (ward %in% trajmwss[[1]]$node)
      trajmwss %<>% lapply(., function(sim)
        sim %<>% .[node == ward])

    trajmwss <- lapply(seq(length(trajmwss)), function(sim) {
      trajmwss[[sim]][, iteration := sim]
      trajmwss[[sim]]
    }) %>% do.call(rbind, .)

    trajmwss[, ':='(
      incP = sum(incPA + incPM + incPS),
      incPsymp = sum(incPM + incPS),
      incH = sum(incHA + incHM + incHS),
      incHsymp = sum(incHM + incHS),
      inc = sum(incPA + incPM + incPS + incHA + incHM + incHS)
    ), by = c("iteration", "time", "node")]

    trajmwss[, ':='(
      mean_incP = mean(incP),
      sd_incP = sd(incP),
      mean_incPsymp = mean(incPsymp),
      sd_incPsymp = sd(incPsymp),
      mean_incH = mean(incH),
      sd_incH = sd(incH),
      mean_incHsymp = mean(incHsymp),
      sd_incHsymp = sd(incHsymp),
      mean_inc = mean(inc),
      sd_inc = sd(inc)
    ), by = c("time", "node")]

    if (pop == "P") {
      setnames(trajmwss, "mean_incP", "meanVal")
      setnames(trajmwss, "sd_incP", "sdVal")
      ylabel = "Cumulative incidence among patients"
    } else
      if (pop == "Psymp") {
        setnames(trajmwss, "mean_incPsymp", "meanVal")
        setnames(trajmwss, "sd_incPsymp", "sdVal")
        ylabel = "Cumulative incidence among patients (only considering symptomatics)"
      } else
        if (pop == "H") {
          setnames(trajmwss, "mean_incH", "meanVal")
          setnames(trajmwss, "sd_incH", "sdVal")
          ylabel = "Cumulative incidence among professionals"
        } else
          if (pop == "Hsymp") {
            setnames(trajmwss, "mean_incHsymp", "meanVal")
            setnames(trajmwss, "sd_incHsymp", "sdVal")
            ylabel = "Cumulative incidence among professionals (only considering symptomatics)"
          } else
            if (isFALSE(pop)) {
              setnames(trajmwss, "mean_inc", "meanVal")
              setnames(trajmwss, "sd_inc", "sdVal")
              ylabel = "Cumulative incidence (patients + professionals)"
            }

    data <-
      trajmwss[, c("time", "node", "meanVal", "sdVal")] %>% unique

    p <- ggplot(data = data,
                aes(
                  x = time,
                  y = meanVal,
                  group = node %>% as.factor,
                  fill = node %>% as.factor
                )) +
      geom_bar(stat = "identity",
               color = "grey",
               width = 1) +
      facet_grid(node %>% as.factor ~ .) +
      xlab("Time (day)") +
      ylab(ylabel) +
      theme(legend.position = "none")

    if (isTRUE(display_sd) & is.logical(iter))
      p <-
      p +  geom_errorbar(
        aes(
          ymin = ifelse(meanVal - sdVal >= 0, meanVal - sdVal, 0),
          ymax = meanVal + sdVal
        ),
        colour = "grey",
        alpha = 0.75,
        width = .2,
        position = position_dodge(.9)
      )
  }

  if (scale == 0) {
    trajmwss <- lapply(seq(length(trajmwss)), function(sim) {
      trajmwss[[sim]][, iteration := sim]
      trajmwss[[sim]]
    })

    trajmwss %<>% do.call(rbind, .)

    trajmwss[, ':='(
      incP = sum(incPA + incPM + incPS),
      incPsymp = sum(incPM + incPS),
      incH = sum(incHA + incHM + incHS),
      incHsymp = sum(incHM + incHS),
      inc = sum(incPA + incPM + incPS + incHA + incHM + incHS)
    ), by = c("iteration", "time")]

    trajmwss[, ':='(
      mean_incP = mean(incP),
      sd_incP = sd(incP),
      mean_incPsymp = mean(incPsymp),
      sd_incPsymp = sd(incPsymp),
      mean_incH = mean(incH),
      sd_incH = sd(incH),
      mean_incHsymp = mean(incHsymp),
      sd_incHsymp = sd(incHsymp),
      mean_inc = mean(inc),
      sd_inc = sd(inc)
    ), by = c("time")]

    if (pop == "P") {
      setnames(trajmwss, "mean_incP", "meanVal")
      setnames(trajmwss, "sd_incP", "sdVal")
      ylabel = "Cumulative incidence among patients"
    } else
      if (pop == "Psymp") {
        setnames(trajmwss, "mean_incPsymp", "meanVal")
        setnames(trajmwss, "sd_incPsymp", "sdVal")
        ylabel = "Cumulative incidence among patients (only considering symptomatics)"
      } else
        if (pop == "H") {
          setnames(trajmwss, "mean_incH", "meanVal")
          setnames(trajmwss, "sd_incH", "sdVal")
          ylabel = "Cumulative incidence among professionals"
        } else
          if (pop == "Hsymp") {
            setnames(trajmwss, "mean_incHsymp", "meanVal")
            setnames(trajmwss, "sd_incHsymp", "sdVal")
            ylabel = "Cumulative incidence among professionals (only considering symptomatics)"
          } else
            if (isFALSE(pop)) {
              setnames(trajmwss, "mean_inc", "meanVal")
              setnames(trajmwss, "sd_inc", "sdVal")
              ylabel = "Cumulative incidence (patients + professionals)"
            }

    data <- trajmwss[, c("time", "meanVal", "sdVal")] %>% unique

    p <- ggplot(data, aes(x = time, y = meanVal)) +
      geom_bar(stat = "identity",
               color = "grey",
               width = 1) +
      xlab("Time (day)") +
      ylab(ylabel) +
      theme(legend.position = "none")

    if (isTRUE(display_sd) & is.logical(iter))
      p <-
      p +  geom_errorbar(
        aes(
          ymin = ifelse(meanVal - sdVal >= 0, meanVal - sdVal, 0),
          ymax = meanVal + sdVal
        ),
        colour = "grey",
        alpha = 0.75,
        width = .2,
        position = position_dodge(.9)
      )
  }

  return(p)
}
