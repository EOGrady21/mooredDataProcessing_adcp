####bin by bin plot###
#'
#'Bin by bin plot
#'
#'
#'use to plot each "bin" of any chosen variable (u, v, error, echo intensity)
#'
#'
#' @param v

plotBin <- function(v){
  for(i in 1:length(v[1, ]))
    plot(v[,i], xlab = "time", ylab = "m/sec", main = (paste( "bin", i)), type = 'l')
  }
