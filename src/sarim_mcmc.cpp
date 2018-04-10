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
//   Main function to sample with a Metropolis-Hasting step 
//
// Written by:
//  Christopher Küster
//
////////////////////////////////////////////////////////////////////////////////



// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

#include <RcppEigen.h>
#include "rue.hpp"
#include "lanczos.hpp"
#include "misc.hpp"
#include "iwls.hpp"
#include "log_likelihood.hpp"

//' MCMC sampler with Metropolis-Hasting step in Gibbs-sampler for use in sarim()-function
//' 
//' This generates samples for the coefficients using a Metropolis-Hasting step in
//' the Gibbs-sampling for the precision parameters and further assumed that the 
//' response variable is either binomial or poisson. So family = "binomial"/"poisson" with
//' link = "logit"/"log" is choosen respectively.
//' 
//' @param y Response variable, given as a vector, use as.numeric() if error occur.
//' @param eta_first First calculation of eta = Z_1 * gamma_1 + ... + Z_p * gamma_p
//' @param Z Design matrices of the covariates, given as list with sparse matrix, 
//'          use e.g. as(matrix, "dgCMatrix") from library(Matrix). Use for example the
//'          useful sx()-function for smoothing.
//' @param K Structure/penalty matrices for the coefficients, given as list, 
//'          also with sparse matrix. Can be choosen in sx()-function.
//' @param gamma List of coefficient, given as vector. Row-length need to be the 
//'          same as the columns of Z. Per starting default from uniform distribution 
//'          is sampled, but a specific starting value can be given, using the sx()-function
//'          in the formular, e.g. y ~ sx(x1, gamma = c(rep(1, 5))).
//' @param kappa Start value for kappa, given as a list and double/float value.
//' @param kappa_values Coefficients for kappa, given as list within as vector c(kappa_a, kappa_b).
//' @param solver List of the solvers ("rue" or "lanczos") for sampling from a gaussian distribution, i.e. gamma ~ N(eta, Q)
//'          with Q as precision matrix. Can be choosen in sx()-function.
//' @param family Can currently only be "binomial" or "poisson"
//' @param link Link need to be chooses correspondingly to the family function. If
//'          family = "binomial", please choose link = "logit" or if family = "poisson", 
//'          please choose link = "log". Currently no other link function is present.
//' @param nIter Number of iterations for MCMC-algorithm
//' @param Ntrials Number of trails, only interesting for binomial distribution
//' @param m Number of maximal Lanczos-iterations
//' @param thr threshold when the Lanczos-algorithm should stop
//' 
//' @return Return a list of values:
//' "coef_results" = result of the estimated coefficient, output given as matrix;
//' "kappa_result" = result of the estimated kappa (precision) parameters, output given as vector
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List sarim_mcmc(const Eigen::Map<Eigen::VectorXd> & y,
                      const Eigen::Map<Eigen::VectorXd> & eta_first,
                      const Rcpp::List & Z,
                      const Rcpp::List & K,
                      const Rcpp::List & gamma,
                      const Rcpp::List & ka_start,
                      const Rcpp::List & ka_values,
                      const Rcpp::List & solver,
                      const Rcpp::String & family,
                      const Rcpp::String & link,
                      const int & nIter,
                      const double & Ntrials = 1,
                      const int & m = 50,
                      const double & thr = 0.0001,
                      const bool & display_progress = true,
                      const Rcpp::String & constraint = "No") {
    
    
    // display progress bar
    Progress pro(nIter, display_progress);
    // number of covariates (p), get from R:List Z
    int p = Z.size();
    // number of observations
    int n = y.rows();
    
    // calculate rank of structure matrices and save in vector K_rk
    Eigen::VectorXd K_rk(p);
    for (int i = 0; i < p; ++i) {
        Eigen::SparseMatrix<double> K_tmp = K(i);
        K_rk(i) = rank_calculation(K_tmp);
    }
    
    
    // generate list for coefficients 
    Rcpp::List coef_results(p);
    Rcpp::List mu_results(p);
    // generate list for potential scale reduction factor
    // and corresponding matrices for coefficient-result-output and psrf-output
    for (int i = 0; i < p; ++i) {
        Eigen::VectorXd gamma_tmp = gamma(i);
        int k_size = gamma_tmp.size();
        Eigen::MatrixXd gamma_results(k_size, 1); // initialise with one column and add +1column iteratively
        gamma_results.col(0) = gamma_tmp;
        coef_results[i] = gamma_results;
        mu_results[i] = gamma_results;
    };
    
    
    Eigen::VectorXd tmp(nIter + 1);
    
    // generate list for kappa
    Rcpp::List kappa_results(p);
    // and corresponding vector for result-output
    for (int i = 0; i < p; ++i) {
        double kappa_tmp = ka_start(i);
        Eigen::VectorXd kappa_results_tmp = tmp;
        kappa_results_tmp(0) = kappa_tmp;
        kappa_results[i] = kappa_results_tmp;
    };
    
    // save eta from first time calculation
    Eigen::VectorXd eta = eta_first;
    // eta_tmp for calculating eta_{-k}
    Eigen::VectorXd eta_tmp = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd eta_tmp_proposal = eta_tmp;
    
    /* too slow?!
     for (int i = 0; i < p; ++i) {
     Eigen::SparseMatrix<double> Z_tmp = Z(i);
     Eigen::VectorXd gamma_tmp = gamma(i);
     eta = eta + Z_tmp * gamma_tmp;
     };
     */
    
    
    // initialise values, required for computation
    Eigen::SparseMatrix<double> Z_k, K_k, Q, W, ZtW;
    Eigen::MatrixXd gamma_matrix, gamma_current;
    
    Eigen::VectorXd b, ka_vector, y_tilde, gamma_proposal, mu, mu_tmp, ka_tmp, x;
    
    Eigen::VectorXd ll, ll_proposal;
    Eigen::VectorXd proposal_c_given_p, proposal_p_given_c, prior_c, prior_p;
    Eigen::VectorXd alpha, u;
    
    
    // Initialise working observations and weights
    for (int k = 0; k < p; ++k) {
        Z_k = Z(k);
        gamma_matrix = coef_results(k);
        ka_vector = kappa_results(k);
        K_k = K(k);
        
        IWLS iwls = compute(y, eta, Ntrials, family, link);
        W = iwls.W;
        y_tilde = iwls.y_tilde;
        
        ZtW = Z_k.transpose() * W;
        b = ZtW * (y_tilde - (eta - Z_k * gamma_matrix.col(0)));
        
        Q = ZtW * Z_k + ka_vector(0) * K_k;
        
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > lltOfQ(Q); 
        Eigen::VectorXd mu = lltOfQ.solve(b);
        
        mu_results[k] = mu;
    }
    
    // loop for number of iterations "nIter" 
    // C++ begin in 0 !!! important for loops
    for (int n_mcmc = 1; n_mcmc <= nIter; ++n_mcmc) {
        //if (Progress::check_abort() )
        //    return -1.0;
        pro.increment();
        
        // loop for covariates p 
        for (int k = 0; k < p; ++k) {
            
            // initialise required matrices z and K and coefficient gamma and kappa
            // set ._k for required calculations
            std::string solv = solver(k);
            Z_k = Z(k);
            K_k = K(k);
            gamma_matrix = coef_results(k);
            gamma_current = gamma_matrix.col(n_mcmc - 1);
            ka_vector = kappa_results(k);
            mu = mu_results(k);
            
            // calculate log-likelihood
            ll = loglike(y, eta, family, link, Ntrials);
            
            // calculate eta_{-k} = eta + Z_k * (mu - gamma) 
            eta_tmp = eta + Z_k * (mu - gamma_current);
            
            // calculate Q = Z'WZ + kappa * K -- and b = Z'_{k} W * (y_tilde - eta_{-k})
            IWLS iwls = compute(y, eta_tmp, Ntrials, family, link);
            W = iwls.W;
            y_tilde = iwls.y_tilde;
            
            // calculate one time for efficiency Z'W
            ZtW = Z_k.transpose() * W;
            
            // calculate Q (precision) matrix 
            Q = ZtW * Z_k + ka_vector(n_mcmc - 1) * K_k;
            
            // calculate b;  b =  Z'_{k} W * (y_tilde - eta_{-k})
            b = ZtW * (y_tilde - (eta_tmp - Z_k * mu) ); 
            
            
            // sample from gaussian distribution
            // depending on solver 'rue' or 'lanczos'
            if (solv == "rue") {
                RueSolv rue_solver = algorithm(Q, b);
                Eigen::VectorXd x = rue_solver.x;
                mu_tmp = rue_solver.mu;
                
                gamma_proposal = x + mu_tmp;
                
            } else {
                Eigen::SparseMatrix<double> M;
                M = ichol(Q);
                Eigen::SparseMatrix<double> Mt = M.transpose();
                Lanczos lanczos_solver = algorithm(Q, m, M, Mt, thr);
                Eigen::VectorXd x = lanczos_solver.x;
                if (lanczos_solver.Iteration >= m) {
                    Rcpp::Rcout << "Number of lanczos-iteration <= " << m <<" in iteration " 
                                << n_mcmc << ". Set m higher (>" << m << ")" << std::endl;
                }
                Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, 
                                         Eigen::Lower|Eigen::Upper, 
                                         Eigen::IncompleteCholesky<double> > iccg(Q);
                mu_tmp = iccg.solve(b);
                
                gamma_proposal = x + mu_tmp;
                
            };
            // compute
            eta_tmp_proposal = eta_tmp + Z_k * (gamma_proposal - mu);
            
            ll_proposal = loglike(y, eta_tmp_proposal, family, link, Ntrials);
            
            // compute p(ga_c | \mu^p, Q^p)
            proposal_c_given_p = -0.5 * ( (gamma_current - mu).transpose() * Q * (gamma_current - mu) );
            
            // compute p(ga_p | \mu^c, Q^c)  
            proposal_p_given_c = -0.5 * ( (gamma_proposal - mu).transpose() * Q * (gamma_proposal - mu) );
            
            // compute p(ga_c | \ka^c)
            prior_c = - (ka_vector.row(n_mcmc - 1) * 0.5) * gamma_current.transpose() * K_k * gamma_current;
            
            // compute p(ga_p | \ka^c)
            prior_p = - (ka_vector.row(n_mcmc - 1) * 0.5) * gamma_proposal.transpose() * K_k * gamma_proposal;
            
            // acceptance probability
            alpha = (ll_proposal + prior_p + proposal_c_given_p) - (ll + prior_c + proposal_p_given_c); 
            
            u = (random_uniform(1)).array().log();
            
            gamma_matrix.conservativeResize(gamma_matrix.rows(), gamma_matrix.cols() + 1);
            if (alpha(0) > u(0)) {
                gamma_matrix.col(n_mcmc) = gamma_proposal;
                coef_results[k] = gamma_matrix;
                mu_results[k] = mu_tmp;
                eta = eta_tmp_proposal;
            } else {
                gamma_matrix.col(n_mcmc) = gamma_current;
                coef_results[k] = gamma_matrix;
                mu_results[k] = mu;
                eta = eta_tmp + Z_k * gamma_current; 
            }
            
            
            
            // update kappa_k by sampling from  ~ Ga(kappa_alpha + rk(K_i)/2  , 
            // kappa_beta + gamma_k ' * K_k * gamma_k / 2)
            ka_tmp = ka_values(k);
            double ka_alpha = ka_tmp.coeff(0, 0) + 0.5 * K_rk(k);
            double ka_beta;
            ka_beta = ka_tmp.coeff(1, 0) + 
                0.5 * (gamma_matrix.col(n_mcmc)).transpose() * K_k * gamma_matrix.col(n_mcmc);
            ka_vector.row(n_mcmc) = random_gamma(1, ka_alpha, 1/ka_beta);
            kappa_results[k] = ka_vector;
        };
        
    };
    
    
    
    // return lists for gamma, kappa and sigma and rank(K) for control
    return Rcpp::List::create(Rcpp::Named("coef_results") = coef_results, 
                              Rcpp::Named("kappa_results") = kappa_results);
}
