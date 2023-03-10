// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// vusC_U
double vusC_U(NumericVector tt1, NumericVector tt2, NumericVector tt3);
RcppExport SEXP _emplikROCS_vusC_U(SEXP tt1SEXP, SEXP tt2SEXP, SEXP tt3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tt1(tt1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tt2(tt2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tt3(tt3SEXP);
    rcpp_result_gen = Rcpp::wrap(vusC_U(tt1, tt2, tt3));
    return rcpp_result_gen;
END_RCPP
}
// vusC_varEL
NumericVector vusC_varEL(NumericVector tt1, NumericVector tt2, NumericVector tt3);
RcppExport SEXP _emplikROCS_vusC_varEL(SEXP tt1SEXP, SEXP tt2SEXP, SEXP tt3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tt1(tt1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tt2(tt2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tt3(tt3SEXP);
    rcpp_result_gen = Rcpp::wrap(vusC_varEL(tt1, tt2, tt3));
    return rcpp_result_gen;
END_RCPP
}
// vusC_full_core
List vusC_full_core(NumericVector tt1, NumericVector tt2, NumericVector tt3);
RcppExport SEXP _emplikROCS_vusC_full_core(SEXP tt1SEXP, SEXP tt2SEXP, SEXP tt3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tt1(tt1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tt2(tt2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tt3(tt3SEXP);
    rcpp_result_gen = Rcpp::wrap(vusC_full_core(tt1, tt2, tt3));
    return rcpp_result_gen;
END_RCPP
}
// place_U
NumericVector place_U(NumericVector tt1, NumericVector tt2, NumericVector tt3);
RcppExport SEXP _emplikROCS_place_U(SEXP tt1SEXP, SEXP tt2SEXP, SEXP tt3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tt1(tt1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tt2(tt2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tt3(tt3SEXP);
    rcpp_result_gen = Rcpp::wrap(place_U(tt1, tt2, tt3));
    return rcpp_result_gen;
END_RCPP
}
// vusC_ties
double vusC_ties(NumericVector tt1, NumericVector tt2, NumericVector tt3);
RcppExport SEXP _emplikROCS_vusC_ties(SEXP tt1SEXP, SEXP tt2SEXP, SEXP tt3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tt1(tt1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tt2(tt2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tt3(tt3SEXP);
    rcpp_result_gen = Rcpp::wrap(vusC_ties(tt1, tt2, tt3));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_emplikROCS_vusC_U", (DL_FUNC) &_emplikROCS_vusC_U, 3},
    {"_emplikROCS_vusC_varEL", (DL_FUNC) &_emplikROCS_vusC_varEL, 3},
    {"_emplikROCS_vusC_full_core", (DL_FUNC) &_emplikROCS_vusC_full_core, 3},
    {"_emplikROCS_place_U", (DL_FUNC) &_emplikROCS_place_U, 3},
    {"_emplikROCS_vusC_ties", (DL_FUNC) &_emplikROCS_vusC_ties, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_emplikROCS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
