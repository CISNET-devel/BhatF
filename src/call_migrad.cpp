#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream> // DEBUG

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


extern "C" void ftn_set_variables (const int* nvars_total,
                                   const char x_labels[],
                                   const double x_ini[],
                                   const int is_fixed[],
                                   const double x_est[],
                                   const double x_min[],
                                   const double x_max[]);


extern "C" void ftn_migrad(int* nvars,
                           double x_final[],
                           double* y_final);


extern "C"
void cpp_callback(double* f,
                  const int* /*iflag*/, const double u[], const int* npar)
{
    // std::cout << "DEBUG: cpp_callback() is called." << std::endl;
    const auto sz = sizeof(data1)/sizeof(*data1);
    assert((sz==sizeof(data2)/sizeof(*data2)) && "Arrays must have the same size");
    Rcpp::NumericVector d1(data1, data1+sz);
    Rcpp::NumericVector d2(data2, data2+sz);

    Rcpp::NumericVector uu(u, u + *npar);
    Rcpp::NumericVector g = uu(0)*(1+d1*uu(1)*(1-uu(2)*d1));
    *f = sum(g - d2*log(g));
    // std::cout << "DEBUG: cpp_callback() is returning *f=" << *f << std::endl;
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
//' @param vars Matrix of variables, 1 row/variable containing is_fixed,initial,estimated,min,max
//' @export
//[[Rcpp::export]]
std::vector<double> call_migrad(const Rcpp::NumericMatrix vars) {
    int nvars = vars.nrow(); // fixme: use `auto` and check for MAX_INT
    if (vars.ncol()!=5) {
        throw std::runtime_error("Variable matrix should have 5 columns: "
                                 "is_fixed,initial,estimated,min,max");
    }
    // it is not said anywhere that Rcpp::NumericVector is a continuous storage!
    // so, we are copying the column vectors to C++ vectors
    std::vector<double> x_ini(vars.column(1).cbegin(), vars.column(1).cend());
    std::vector<double> x_est(vars.column(2).cbegin(), vars.column(2).cend());
    std::vector<double> x_min(vars.column(3).cbegin(), vars.column(3).cend());
    std::vector<double> x_max(vars.column(4).cbegin(), vars.column(4).cend());

    std::vector<int> is_fixed(nvars);
    std::transform(vars.column(0).cbegin(), vars.column(0).cend(),
                   is_fixed.begin(),
                   [](double x) { return x!=0; } );

    // const int npar = std::count(is_fixed.begin(), is_fixed.end(), 0);

    Rcpp::CharacterVector r_labels = rownames(vars);
    std::vector<char> f_labels(nvars*10, ' ');
    for (auto i=0; i<nvars; ++i) {
        std::copy(r_labels[i].begin(), r_labels[i].end(), &f_labels[i*10]);
    }

    // std::cout << "DEBUG: calling ftn_set_variables()" << std::endl;
    ftn_set_variables(&nvars, f_labels.data(), x_ini.data(), is_fixed.data(), x_est.data(), x_min.data(), x_max.data());


    std::vector<double> yx_final(nvars+1);
    
    // std::cout << "DEBUG: calling ftn_migrad()" << std::endl;
    ftn_migrad(&nvars, yx_final.data()+1, &yx_final[0]);
    
    // std::cout << "DEBUG: back from migrad(), ymin=" << yx_final[0] <<  std::endl;

    // std::cout << "DEBUG: back from cpp_callback()" << std::endl;
    return yx_final;
}
