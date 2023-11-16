# Functions to plot operation codes

library(data.table)
library(ggplot2)
library(codeCollection)

# NOMESCO operation codes to data.table
NOM <- as.data.table(TPKoodit)

#' Create barplot of operation code counts with understandable names to codes
#' 
#' @param x A vector of operation codes
#' @param lang Language of code descriptions (fi, en, sv)
#' @param mapping Mapping of NOMESCO codes to descriptions
#' @returns A ggplot object
#' @examples
#' add(1, 1)
ggplot_geom_bar_operations <- function(x, lang) {
    X <- data.table(CODE = x)
    X <- X[, .N, by = .(CODE)]
    X <- merge(X, NOM, by.x = "CODE", by.y = "Koodi")[order(N)]

    ggplot(X) +
        geom_bar(stat = "identity", mapping = aes(x = CODE, y = N)) +
        theme_minimal() +
        ggtitle("Operations by code") +
        scale_x_discrete(limits = X[, CODE],
                        breaks = X[, CODE],
                        labels = X[, get(lang)]) +
        coord_flip()
}
