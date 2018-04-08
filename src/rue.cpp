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
//   Definition for Rue (2001) algorithm to sample form gaussian distribution
//
// Written by:
//  Christopher Küster
//
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "rue.hpp"
#include "misc.hpp"

// Rue-algorithm
//
// first decompose Q by Cholesky, i.e. Q = LL'
// a sample x ~ N(0, Q^{-1}) can then be obtained by sampling z from N(0, I)
// and solving 

RueSolv algorithm (const Eigen::SparseMatrix<double> & Q, const Eigen::VectorXd & b) {
    
    RueSolv rue_solver;
    // number of rows in Q, necessary for sampling from z
    int k = Q.rows();
    
    // initialise vectors for solving system
    Eigen::VectorXd x(k), w(k), mu(k);
    
    // sample z from N(0, I)
    Eigen::VectorXd z = random_gauss(k);
    
    // compute Cholesky decomposition of Q and save as "Lt"
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > lltOfQ(Q); 
    Eigen::SparseMatrix<double> Lt = lltOfQ.matrixU();
    Lt.makeCompressed(); // compression required for SparseQR
    
    // solving "L'x = z" to get a sample from N(0, Q^{-1})
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > Ltsolve(Lt);
    x = Ltsolve.solve(z);
    
    // solving "Lw = b"
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > Lsolve(Lt.transpose());
    w = Lsolve.solve(b);
    
    // solving "L'mu = w"
    mu = Ltsolve.solve(w);
    
    rue_solver.x = x;
    rue_solver.mu = mu;
    return rue_solver;
};