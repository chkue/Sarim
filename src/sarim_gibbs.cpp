// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

#include <RcppEigen.h>
#include "rue.hpp"
#include "lanczos.hpp"
#include "misc.hpp"

//' Gibbs sampler for use in sarim()-function
//' 
//' This generates samples for the coefficients using a gibbs sampler, 
//' assumed the response variable is a gaussian distribution, so family = "gaussian" and
//' link = "identity" is choosen in sarim()-function
//' 
//' @param y Response variable, given as a vector, use as.numeric() if error occur.
//' @param Z Design matrices of the covariates, given as list with sparse matrix, 
//'          use e.g. as(matrix, "dgCMatrix") from library(Matrix). Use for example the
//'          useful sx()-function for smoothing.
//' @param K Structure/penalty matrices for the coefficients, given as list, 
//'          also with sparse matrix. Can be choosen in sx()-function.
//' @param K_rank List of ranks of the structure/penalty matrix.
//' @param gamma List of coefficient, given as vector. Row-length need to be the 
//'          same as the columns of Z. Per starting default from uniform distribution 
//'          is sampled, but a specific starting value can be given, using the sx()-function
//'          in the formular, e.g. y ~ sx(x1, gamma = c(rep(1, 5))).
//' @param kappa Start value for kappa, given as a list and double/float value.
//' @param kappa_values Coefficients for kappa, given as list within as vector c(kappa_a, kappa_b).
//' @param solver List of the solvers ("rue" or "lanczos") for sampling from a 
//'          gaussian distribution, i.e. gamma ~ N(eta, Q)
//'          with Q as precision matrix. Can be choosen in sx()-function.
//' @param sigma Start value for sigma given as numeric value, variance of response y. 
//' @param sigma_values Values for sigma ~ IG(sigma_a, sigma_b) given as a vector c(sigma_a, sigma_b), 
//'          similar to kappa_values. Can be also choosen in sx()-function.
//' @param nIter Number of iterations for MCMC-algorithm
//' @param m Number of maximal Lanczos-iterations
//' @param thr threshold when the Lanczos-algorithm should stop
//' 
//' @return Return a list of values:
//' "coef_results" = result of the estimated coefficient, output given as matrix; 
//' "kappa_result" = result of the estimated kappa (precision) parameters, output given as vector;
//' "sigma_results" = result of the sigma value, output given as vector
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List sarim_gibbs(const Eigen::Map<Eigen::VectorXd> & y,
                       const Rcpp::List & Z,
                       const Rcpp::List & K,
                       const Rcpp::List & K_rank,
                       const Rcpp::List & gamma,
                       const Rcpp::List & ka_start,
                       const Rcpp::List & ka_values,
                       const Rcpp::List & solver,
                       const double & sigma,
                       const Eigen::Map<Eigen::VectorXd> & sigma_values,
                       const int & nIter,
                       const int & m,
                       const double & thr,
                       const bool & display_progress = true) {
    
    
    // display progress bar
    Progress pro(nIter, display_progress);
    // number of covariates (p), get from R:List Z
    int p = Z.size();
    // number of observations
    int n = y.rows();
    
    // generate list for coefficients 
    Rcpp::List coef_results(p);
    // and corresponding matrices for coefficient-result-output
    for (int i = 0; i < p; ++i) {
        Eigen::VectorXd gamma_tmp = gamma(i);
        int k_size = gamma_tmp.size();
        Eigen::MatrixXd gamma_results = Eigen::MatrixXd::Zero(k_size, nIter + 1);
        gamma_results.col(0) = gamma_tmp;
        coef_results[i] = gamma_results;
    };
    
    // generate list for kappa
    Rcpp::List kappa_results(p);
    // and corresponding vector for result-output
    for (int i = 0; i < p; ++i) {
        double kappa_tmp = ka_start(i);
        Eigen::VectorXd kappa_results_tmp(nIter + 1);
        kappa_results_tmp(0) = kappa_tmp;
        kappa_results[i] = kappa_results_tmp;
    };
    
    // calculate eta the first time, eta = Z_1 * gamma_1 + ... + Z_p * gamma_p
    Eigen::VectorXd eta = Eigen::VectorXd::Zero(n);
    // eta_tmp for calculating eta_{-k}
    Eigen::VectorXd eta_tmp = eta;
    
    
    for (int i = 0; i < p; ++i) {
        Eigen::SparseMatrix<double> Z_tmp = Z(i);
        Eigen::VectorXd gamma_tmp = gamma(i);
        eta = eta + Z_tmp * gamma_tmp;
    };
    
    // initialise results vector for sigma and set first value to starting value
    Eigen::VectorXd sigma_results(nIter + 1);
    sigma_results(0) = sigma;
    
    // initialise values, required for computation
    Eigen::SparseMatrix<double> Z_k, K_k, Q, M, Mt;
    Eigen::MatrixXd gamma_matrix;
    Eigen::VectorXd K_rk, gamma_matrix_tmp, ka_vector, ka_tmp, b;
    Eigen::VectorXd x, mu_tmp, gamma_proposal;
    
    
    // loop for number of iterations "nIter" 
    // C++ begin in 0 !!! important for loops
    for (int n_mcmc = 1; n_mcmc <= nIter; ++n_mcmc) {
        //if (Progress::check_abort() )
        //    return -1.0;
        pro.increment();
        
        // loop for covariates p 
        for (int k = 0; k < p; ++k) {
            
            // initialise required matrices Z and K and coefficient gamma and kappa
            // set ._tmp for required calculations
            std::string solv = solver(k);
            Z_k = Z(k);
            K_k = K(k);
            K_rk = K_rank(k);
            gamma_matrix = coef_results(k);
            ka_vector = kappa_results(k);
            
            // calculate eta_{-k} = eta - Z_k * gamma
            eta_tmp = eta - Z_k * gamma_matrix.col(n_mcmc - 1);
            
            // calculate Q (precision) matrix ; 
            Q = (1/sigma_results(n_mcmc - 1)) * Z_k.transpose() * Z_k + 
                ka_vector(n_mcmc - 1) * K_k;
            
            // calculate b;  b = (1/var(y)) * Z'_{k} * (y - eta_{-k})
            b = 1/sigma_results(n_mcmc - 1) * Z_k.transpose() * (y - eta_tmp); 
            
            // sample from gaussian distribution
            // depending on solver 'rue' or 'lanczos'
            if (solv == "rue") {
                RueSolv rue_solver = algorithm(Q, b);
                x = rue_solver.x;
                mu_tmp = rue_solver.mu;
                
                gamma_proposal = x + mu_tmp;
                
            } else {
                M = ichol(Q);
                Mt = M.transpose();
                Lanczos lanczos = algorithm(Q, m, M, Mt, thr);
                
                if (lanczos.Iteration >= m) {
                    Rcpp::Rcout << "Number of lanczos-iteration <= " << m <<" in iteration " 
                                << n_mcmc << ". Set m higher (>" << m << ")" << std::endl;
                }
                
                x = lanczos.x;
                
                Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, 
                                         Eigen::Lower|Eigen::Upper, 
                                         Eigen::IncompleteCholesky<double> > iccg(Q);
                mu_tmp = iccg.solve(b);
                
                gamma_proposal = x + mu_tmp;
            };
            
            gamma_matrix.col(n_mcmc) = gamma_proposal;
            coef_results[k] = gamma_matrix;
            eta = eta_tmp + Z_k * gamma_matrix.col(n_mcmc);
            
            // update kappa_k by sampling from  ~ Ga(kappa_alpha + rk(K_i)/2  , 
            // kappa_beta + gamma_k ' * K_k * gamma_k / 2)
            ka_tmp = ka_values(k);
            double ka_alpha = ka_tmp.coeff(0, 0) + 0.5 * K_rk(0);
            double ka_beta;
            ka_beta = ka_tmp.coeff(1, 0) + 
                0.5 * (gamma_matrix.col(n_mcmc)).transpose() * K_k * gamma_matrix.col(n_mcmc);
            ka_vector.row(n_mcmc) = random_gamma(1, ka_alpha, 1/ka_beta);
            kappa_results[k] = ka_vector;
        };
        
        // updata sigma, which is "~ IG(sigma_a + n/2, 1/2 * (y - eta)'(y - eta))"
        double sigma_alpha = sigma_values(0) + n/2;
        double sigma_beta = sigma_values(1) + 0.5 * (y - eta).transpose() * (y - eta);
        sigma_results.row(n_mcmc) = random_gamma(1, sigma_alpha, 1/sigma_beta).cwiseInverse();
    };
    
    
    
    // return lists for gamma, kappa and sigma and rank(K) for control
    return Rcpp::List::create(Rcpp::Named("coef_results") = coef_results, 
                              Rcpp::Named("kappa_results") = kappa_results,
                              Rcpp::Named("sigma_results") = sigma_results);
}