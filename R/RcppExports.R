# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

linprogPD2C <- function(x, A, b, epsilon, pdtol = 1e-3, pdmaxiter = 50L) {
    .Call('clime2_linprogPD2C', PACKAGE = 'clime2', x, A, b, epsilon, pdtol, pdmaxiter)
}

