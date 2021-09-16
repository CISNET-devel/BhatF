#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <cassert>
// #include <iostream> // DEBUG


extern "C" void ftn_set_variables (const int* nvars_total,
                                   const char x_labels[],
                                   const double x_ini[],
                                   const int is_fixed[],
                                   const double x_est[],
                                   const double x_min[],
                                   const double x_max[]);


extern "C" void ftn_migrad(int* nvars,
                           double x_final[],
                           double* y_final,
                           int* ncalls,
                           int* status);


using objective_fn_t = std::function<double(const std::vector<double>& u, int iflag)>;

// Test function to be minimized
double test_fn(const std::vector<double>& u, int /*iflag*/) {
    static const double data1[] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,5,5,5,5,
        5,5,5,5,5,5,5,5,5,5,5,5,5,
        5,5,5,5,5,5,5,5,5,5,5,5,5,
        5,5,5,5,5,5,5,5,5,5,5,5,5,
        5,5,5,5,5,5,5,10,10,10,10,10,10,
        10,10,10,10,10,10,10,10,10,10,10,10,10,
        10,10,10,10,10,10,10,10,10,10,10,10,10,
        10,10,10,10,10,10,10,10,10,10,10,10,10,
        10,10,10,10,10
    };

    static const double data2[] = {
        23,23,15,22,21,13,20,19,22,25,
        18,19,15, 21,20,19,19,18,21,25,
        31,18,19,20,20,17,19,19,21,27,
        30,17,14,24,21,22,22,13,23,15,
        19,24,18,30,16,21,15,28,20,31,
        38,38,32,22,31,38,26,42,32,26,
        29,28,31,22,30,26,41,41,33,33,
        39,35,13,31,30,28,18,23,33,27,
        40,26,28,27,32,34,28,32,29,34,
        23,32,26,29,21,32,31,31,41,26,
        67,60,60,44,59,75,38,38,45,56,
        46,57,55,62,55,56,73,52,50,57,
        52,50,61,69,65,38,52,62,59,65,
        68,68,64,54,64,54,52,71,59,46,
        51,59,60,45,49,41,56,61,59,48,
        68,71,60,78,68,63,75,73,57,70,
        64,82,73,75,79,59,83,79,77,47,
        66,60,89,64,90,67,78,72,65,81,
        82,59,75,66,63,80,77,67,59,66,
        80,66,74,71,69,85,70,70,66,63
    };
    
    const auto sz = sizeof(data1)/sizeof(*data1);
    assert((sz==sizeof(data2)/sizeof(*data2)) && "Arrays must have the same size");
    Rcpp::NumericVector d1(data1, data1+sz);
    Rcpp::NumericVector d2(data2, data2+sz);

    Rcpp::NumericVector uu(u.begin(), u.end());
    Rcpp::NumericVector g = uu(0)*(1+d1*uu(1)*(1-uu(2)*d1));
    return sum(g - d2*log(g));
}


// FIXME!!!  This should be assigned when calling call_migrad()
objective_fn_t ObjectiveFn {test_fn} ;   // FIXME: bare global --- hide it in a class or namespace
                                        // Warning: the whole thing is not re-enterable because of the global!

extern "C"
void cpp_callback(double* f,
                  const int* iflag, const double u[], const int* npar)
{
    *f = ObjectiveFn(std::vector<double>(u, u+*npar), *iflag);
}


extern "C"
void message_callback(const char* message, const int* msglen)
{
    std::string msg(message, *msglen);
    // Rcpp::Rcout << msg << "\n";
    Rcpp::Function r_message("message");
    r_message(msg);
}


//' call_migrad
//'
//' Call migrad subroutine
//'
//' @param vars List with components:
//' 'label' (character vector),
//' 'est' (initial guess vector),
//' 'low' (lower bounds vector),
//' 'upp' (upper bounds vector)
//' 'fixed' (optional, logicals vector, TRUE if variable to be held constant)
//' @param fn R Function to be optimized
//' @export
//[[Rcpp::export]]
Rcpp::List call_migrad(const Rcpp::List vars, Rcpp::Function fn) {
    using Rcpp::as;
    using Rcpp::NumericVector;
    
    ObjectiveFn = [fn](const std::vector<double>& x, int) { return as<double>(fn(x)); };
    
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
