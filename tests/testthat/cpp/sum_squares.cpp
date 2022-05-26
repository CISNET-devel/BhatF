#include <Rcpp.h>
#include "../../../inst/include/BhatF_helper.hpp"


// Returns 1.25+Σ(xᵢ-1)²
double test_sum_squares_cpp(const std::vector<double>& u, int /*iflag*/) {
    double result = 1.25;
    for (const auto& x: u) { result += (x-1)*(x-1); }
    return result;
}

// Returns the address of the C++ function
//[[Rcpp::export]]
BhatF_caller::xptr_t get_sum_squares_xptr() {
    using namespace BhatF_caller;
    return fptr2xptr(test_sum_squares_cpp);
}
