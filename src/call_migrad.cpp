#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <cassert>
// #include <iostream> // DEBUG

#include <BhatF_helper.hpp>


/// Set Fortran globals (common-block variables) for the optimization problem
extern "C" void ftn_set_variables (const int* nvars_total,
                                   const char x_labels[],
                                   const double x_ini[],
                                   const int is_fixed[],
                                   const double x_est[],
                                   const double x_min[],
                                   const double x_max[]);

/// Call Fortran migrad() function to solve the optimization problem
extern "C" void ftn_migrad(int* nvars,
                           double x_final[],
                           double* y_final,
                           int* ncalls,
                           int* status);


// FIXME:ToDo: Encapsulate "what should be called after what" in some class ctor?

namespace migrad_caller {
    using namespace BhatF_caller;
    
    static inline double objective_fn_stub(const std::vector<double>&, int) {
        // Rcpp::Rcerr << "The objective function is left uninitialized!\nAborting...";
        Rcpp::Function r_stop("stop");
        r_stop("The objective function is left uninitialized!\nStopping...");
        return 0;
    } 

    // FIXME!!!  This should be assigned prior to calling call_migrad()
    //           ToDo: Encapsulate it in some class ctor to enforce?
    // Warning: the optimization is not re-enterable because of this static object!
    // Note: the optimization is not re-enterable anyway because of Fortran common blocks.
    objective_fn_t ObjectiveFn {objective_fn_stub};
}

/// Called from Fortran when it needs to evaluate the objective function
extern "C"
void cpp_callback(double* f,
                  const int* iflag, const double u[], const int* npar)
{
    // ToDo: better to use a `span`
    *f = migrad_caller::ObjectiveFn(std::vector<double>(u, u+*npar), *iflag);
}


/// Called from Fortran when it needs to print a message
extern "C"
void message_callback(const char* message, const int* msglen)
{
    std::string msg(message, *msglen);
    // Rcpp::Rcout << msg << "\n";
    Rcpp::Function r_message("message");
    r_message(msg);
}


//' dfp
//'
//' Call Davidon-Fletcher-Powell optimization Fortran routine
//'
//' @param vars List with components:
//' 'label' (character vector),
//' 'est' (initial guess vector),
//' 'low' (lower bounds vector),
//' 'upp' (upper bounds vector)
//' 'fixed' (optional, logicals vector, TRUE if variable to be held constant)
//' @param fn R Function or a pointer to Rcpp function to be optimized
//' @export
// [[Rcpp::export]]
Rcpp::List dfp(const Rcpp::List vars, Rcpp::RObject fn) {
    using Rcpp::as;
    using Rcpp::is;
    using Rcpp::NumericVector;
    using Rcpp::Function;

    if (is<Function>(fn)) {
        migrad_caller::ObjectiveFn = [fn](const std::vector<double>& x, int) { return as<double>(as<Function>(fn)(x)); };
    }
    else if (TYPEOF(fn)==EXTPTRSXP) {
        migrad_caller::xptr_t xptr(fn);
        migrad_caller::ObjectiveFn = *xptr;
    } else {
        throw std::runtime_error("The argument `fn` must be an R function"
                                 " or a pointer to a native function");
    }

    
    // it is not said anywhere that Rcpp::NumericVector is a continuous storage!
    // so, we are copying the R vectors to C++ vectors
    std::vector<double> x_est(as<NumericVector>(vars["est"]).cbegin(), as<NumericVector>(vars["est"]).cend());
    std::vector<double> x_ini(x_est);
    std::vector<double> x_min(as<NumericVector>(vars["low"]).cbegin(), as<NumericVector>(vars["low"]).cend());
    std::vector<double> x_max(as<NumericVector>(vars["upp"]).cbegin(), as<NumericVector>(vars["upp"]).cend());
    int nvars = x_est.size(); // FIXME: use `auto` and check for MAX_INT

    if (x_min.size()!=static_cast<unsigned>(nvars) || x_max.size()!=static_cast<unsigned>(nvars)) {
        throw std::runtime_error("Inconsistent number of values"); // FIXME: better error message
    }

    std::vector<int> is_fixed(nvars);
    if (vars.containsElementNamed("fixed")) {
        Rcpp::LogicalVector r_is_fixed=vars["fixed"];
        if (r_is_fixed.size()!=nvars) {
            throw std::runtime_error("Inconsistent number of values defining fixed variables"); // FIXME: better error message
        }
        std::transform(r_is_fixed.cbegin(), r_is_fixed.cend(),
                       is_fixed.begin(),
                       [](int x) { return !!x; } );
    }
    // const int npar = std::count(is_fixed.begin(), is_fixed.end(), 0);

    Rcpp::CharacterVector r_labels = vars["label"];
    if (r_labels.size()!=nvars) {
        throw std::runtime_error("Inconsistent number of variable labels"); // FIXME: better error message
    }
    std::vector<char> f_labels(nvars*10, ' ');
    for (auto i=0; i<nvars; ++i) {
        std::copy(r_labels[i].begin(), r_labels[i].end(), &f_labels[i*10]);
    }

    ftn_set_variables(&nvars, f_labels.data(), x_ini.data(), is_fixed.data(), x_est.data(), x_min.data(), x_max.data());


    std::vector<double> x_final(nvars);
    double y_final {};
    int status=-1;
    int ncalls=0;
    ftn_migrad(&nvars, x_final.data(), &y_final, &ncalls, &status);

    using Rcpp::_;
    return Rcpp::List::create(_("fmin")=y_final,
                              _("est")=x_final,
                              _("status")=status,
                              _("nfcn")=ncalls,
                              _("label")=r_labels);
}
