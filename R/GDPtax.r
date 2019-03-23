#' GDP per capita with country groups and profit tax
#'
#' GDP per capita with country groups and profit tax
#'
#'
#' The data includes the GDP per capita (current US$) for years 1980-2017.
#' The data is World Bank national accounts data, and OECD National Accounts data files.
#' The data have been downloaded from the webpage
#' \url{https://datamarket.com/data/set/15c9/gdp-per-capita-current-us#!ds=15c9!hd1&display=line},
#' distributed under the CC-BY 4.0 (\url{https://datacatalog.worldbank.org/public-licenses#cc-by}).
#' From the same webpage the profit tax in 2010 was downloaded.
#' Furthermore, different country groups were formed from
#' countries for which the GDP was available for 1980-2017 and profit tax for 2010:
#' \itemize{
#' \item Group 1 (Major Advanced Economies (G7)): "Canada", "France", "Germany", "Italy", Japan"
#' \item Group 2 (Euro Area excluding G7): "Austria", "Belgium", "Cyprus", "Finland", "Greece", "Ireland",
#' "Luxembourg", "Netherlands", "Portugal", "Spain"
#' \item Group 3 (Other Advanced Economies (Advanced Economies excluding G7 and Euro Area)):
#' "Australia", "Denmark", "Iceland", "Norway", "Sweden", "Switzerland"
#' \item Group 4 (Emerging and Developing Asia): "Bangladesh", "Bhutan", "China", "Fiji", "India",
#' "Indonesia", "Malaysia", "Nepal", "Philippines", "Thailand", "Vanuatu"
#' }
#'
#' @format A list of a \code{curve_set} object with components \code{r} and \code{obs}
#' containing the years of observations and the GDP curves, i.e. the observed values of GDP in
#' those years. Each column of \code{obs} contains the GDP for the years for a particular
#' country (seen as column names). The country grouping is given in the list component \code{Group}
#' and the profit tax in \code{Profittax}.
#'
#' @usage data(GDPtax)
#' @references
#'
#' @keywords datasets
#' @keywords curves
#' @name GDPtax
#' @docType data
#' @seealso \code{\link{graph.fglm}}
#' @examples
#' data(GDPtax)
#' GDPcset <- GDPtax$curve_set
#' # Plot data in groups
#' subs <- function(group, ...) {
#'   cset <- GDPcset
#'   cset$obs <- GDPcset$obs[, GDPtax$Group == group]
#'   plot(cset, ...)
#' }
#' par(mfrow=c(2,2))
#' for(i in 1:4) subs(i, main=paste("Group ", i, sep=""), ylab="GDP")
#'
#' # See example analysis in ?graph.fglm
NULL
