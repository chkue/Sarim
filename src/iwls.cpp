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
//   Definition for iterative weighted least square functions
//
// Written by:
//   Christopher Küster
//
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "iwls.hpp"
#include "response_function.hpp"


// structure to calculate the working dependent variable and weight matrix
IWLS compute (const Eigen::VectorXd & y, 
              const Eigen::VectorXd & eta, 
              const std::string & fam, 
              const std::string & link,
              const double & Ntrials) {
    // make structure for multiple output
    IWLS iwls;
    
    int n = y.rows(); // number of observations
    // initialisation 
    Eigen::VectorXd y_t, W_tmp, h; 
    Eigen::SparseMatrix<double> W_t(n, n);
    
    
    if (fam == "poisson") {
        if (link == "log") {
            // calcualte response function which also the weight function
            h = response_function(eta, link);
            y_t = eta + (y - h).cwiseQuotient(h);
            
            // insert weights in diagonal matrix
            for (int i = 0; i < n; ++i) {
                W_t.insert(i, i) = h(i);
            }
            
            // compress sparse matrix due to better calculation for Eigen
            W_t.makeCompressed();
        }
    }
    
    
    if (fam == "binomial") {
        if (link == "logit") {
            // calcualte response function 
            // and attempt simple trick to preserve from numerical problems
            h = response_function(eta, link);
            Eigen::VectorXd tmp2 = response_function(-eta, link);
            for (int i = 0; i < n; ++i) {
                if (h(i) == 1)  h(i) = 1 - 1e-7;
                if (h(i) == 0)  h(i) = 1e-7;
                if (tmp2(i) == 1) tmp2(i) = 1 - 1e-7;
                if (tmp2(i) == 0) tmp2(i) = 1e-7;
            };
            
            // working observations
            y_t = eta.array() + (y - Ntrials * h).array() * ( 1 / (tmp2.array() * (1 - tmp2.array()) * Ntrials) ).array();
            
            // weights for weight matrix
            Eigen::VectorXd w = ( (Ntrials * tmp2.array() * (1 - tmp2.array()) ).array().square() ).array() / (Ntrials * h.array() * (1 - h.array()) ).array(); 
            for (int i = 0; i < n; ++i) {
                W_t.insert(i, i) = w(i);
            }
            
            // compress sparse matrix due to better calculation for Eigen
            W_t.makeCompressed();
        }
    }
    
    iwls.y_tilde = y_t;
    iwls.W = W_t;
    return iwls;
}
