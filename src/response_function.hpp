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
//   Response function for logit or log-function
//
// Written by:
//   Christopher Küster
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RESPONSEFUNC_H_
#define RESPONSEFUNC_H_

#include <RcppEigen.h>

Eigen::VectorXd response_function (const Eigen::VectorXd & eta, 
                                   const std::string & link);


#endif // IWLS_H_