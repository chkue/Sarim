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
//   Declaration for iterative sampling from a gaussian distribution by Lanczos
//
// Written by:
//  Christopher Küster
//
////////////////////////////////////////////////////////////////////////////////

#ifndef LANCZOS_H_
#define LANCZOS_H_

#include <RcppEigen.h>


Eigen::MatrixXd tri_mat_creation (const int & n, 
                                  const Eigen::VectorXd & a, 
                                  const Eigen::VectorXd & b);


struct Lanczos 
{
    Eigen::VectorXd x;
    Eigen::VectorXd error;
    int Iteration;
};


Lanczos algorithm (const Eigen::SparseMatrix<double> & Q, 
                   const int & m, 
                   const Eigen::SparseMatrix<double> & F1, 
                   const Eigen::SparseMatrix<double> & F2, 
                   const double & thr);

#endif // LANCZOS_H_