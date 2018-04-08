////////////////////////////////////////////////////////////////////////////////
//                  _
//     ________ ___(_) ___
//    / __/ _  | __| |/   \
//    \__ \(_) | | | | Y Y |
//    /___/\_|_|_| |_|_|_|_|
//
//
// sarim is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// sarim is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with sarim. If not, see <http://www.gnu.org/licenses/>.
//
// This file contains:
// -------------------
//
//   Definition for calculation of the log-likelihood
//
// Written by:
//  Christopher Küster
//
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "log_likelihood.hpp"
#include "iwls.hpp"


// Compute the loglikelihood depending on the family and link function 
Eigen::VectorXd loglike (const Eigen::VectorXd & y, 
                         const Eigen::VectorXd & eta, 
                         const std::string & family,
                         const std::string & link,
                         const double & Ntrials) {
    
    // Intialise output vectors
    Eigen::VectorXd out, h_tmp;
    int nrows = eta.rows();
    
    Eigen::VectorXd h = response_function(eta, link);
    
    if (family == "poisson") {
        out = ( y.array() * eta.array() - h.array() ).colwise().sum();
    }
    
    if (family == "binomial") {
        
        /*  Try to make log-likelihood independent of starting values in gamma
         *  for (int i = 0; i < nrows; ++i) {
         *      if (h(i) == 1) h(i) = 1 - 1e-10;
         *      if (h(i) == 0) h(i) = 1e-10;
         *  }
         */
        out = ( y.array() * h.array().log() + (Ntrials - y.array()).array() * (1 - h.array() ).array().log() ).colwise().sum();
    }
    
    return out;
}




