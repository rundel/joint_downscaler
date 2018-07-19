// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(local_include)]]

#include <RcppArmadillo.h>

#include <assert.hpp>
#include <distance.hpp>


arma::mat zero_trunc(arma::mat const& v)
{
    return (v+arma::abs(v))/2.0;
}

// [[Rcpp::export]]
Rcpp::List downscaler_joint_predict_phi(arma::mat Q_P, arma::mat s_M, arma::mat s_P,
                                        arma::mat phi, arma::vec xi,
                                        arma::vec eta2, arma::mat sigma2,
                                        arma::mat beta0, arma::mat beta1,
                                        arma::mat beta0s_1, arma::mat beta0s_2,
                                        arma::mat beta0s_3, arma::mat beta0s_4,
                                        arma::mat beta0s_5, arma::mat beta0s_o, 
                                        std::string trans)
{
    int n_P = s_P.n_rows;
    int n_iter = phi.n_rows;
    int n_species = 5;
    int other = 5;

    arma::mat log_Q_P;
    if (trans=="log")
        log_Q_P = arma::log(Q_P);

    // Data checks

    RT_ASSERT(phi.n_rows == xi.n_rows,"");
    RT_ASSERT(phi.n_rows == eta2.n_rows,"");
    RT_ASSERT(phi.n_rows == sigma2.n_rows,"");
    RT_ASSERT(phi.n_rows == beta0.n_rows,"");
    RT_ASSERT(phi.n_rows == beta1.n_rows,"");

    RT_ASSERT(phi.n_rows == beta0s_1.n_rows,"");
    RT_ASSERT(phi.n_rows == beta0s_2.n_rows,"");
    RT_ASSERT(phi.n_rows == beta0s_3.n_rows,"");
    RT_ASSERT(phi.n_rows == beta0s_4.n_rows,"");
    RT_ASSERT(phi.n_rows == beta0s_5.n_rows,"");
    RT_ASSERT(phi.n_rows == beta0s_o.n_rows,"");

    RT_ASSERT(phi.n_cols == n_species,"");
    RT_ASSERT(sigma2.n_cols == n_species,"");
    RT_ASSERT(beta0.n_cols == n_species+1,"");
    RT_ASSERT(beta1.n_cols == n_species+1,"");

    RT_ASSERT(beta0s_2.n_cols == beta0s_1.n_cols,"");
    RT_ASSERT(beta0s_3.n_cols == beta0s_1.n_cols,"");
    RT_ASSERT(beta0s_4.n_cols == beta0s_1.n_cols,"");
    RT_ASSERT(beta0s_5.n_cols == beta0s_1.n_cols,"");
    RT_ASSERT(beta0s_o.n_cols == beta0s_1.n_cols,"");

    // beta0s

    std::vector< arma::mat > beta0s(n_species+1);

    beta0s[0] = beta0s_1;
    beta0s[1] = beta0s_2;
    beta0s[2] = beta0s_3;
    beta0s[3] = beta0s_4;
    beta0s[4] = beta0s_5;
    beta0s[5] = beta0s_o;


    // Distance Matrices
    
    arma::mat d_M = circle_dist( s_M );
    arma::mat d_P = circle_dist( s_P );
    arma::mat d_M_P = circle_dist( s_M, s_P );

    
    // Output

    arma::mat Z_1_P(n_P, n_iter);
    arma::mat Z_2_P(n_P, n_iter);
    arma::mat Z_3_P(n_P, n_iter);
    arma::mat Z_4_P(n_P, n_iter);
    arma::mat Z_5_P(n_P, n_iter);
    arma::mat Z_o_P(n_P, n_iter);
    arma::mat Z_tot_P(n_P, n_iter);

    arma::mat beta0s_1_P(n_P, n_iter);
    arma::mat beta0s_2_P(n_P, n_iter);
    arma::mat beta0s_3_P(n_P, n_iter);
    arma::mat beta0s_4_P(n_P, n_iter);
    arma::mat beta0s_5_P(n_P, n_iter);
    arma::mat beta0s_o_P(n_P, n_iter);


    std::vector<arma::mat> beta0s_mu(n_species+1);
    std::vector<arma::mat> beta0s_S_U(n_species+1);

    for(int j=0; j != n_iter; ++j)
    {
        //////
        // Covariance Matrices
        /////    

        for(int i = 0; i != n_species+1; ++i)
        {
            double cur_phi;
            if (i == n_species)
            {
                cur_phi = xi(j);
            }
            else
            {
                cur_phi = phi(j,i);
            }

            arma::mat Sigma_M_inv = arma::inv_sympd( arma::exp(-cur_phi * d_M) );
            arma::mat Sigma_P = arma::exp(-cur_phi * d_P);
            arma::mat Sigma_M_P = arma::exp(-cur_phi * d_M_P);

            beta0s_mu[i] = Sigma_M_P.t() * Sigma_M_inv;
            beta0s_S_U[i] = arma::chol(Sigma_P - (beta0s_mu[i] * Sigma_M_P));
        }

        ////
        // Predict \beta_0(s)
        ////

        arma::mat beta0s_P(n_P, n_species+1);   

        for(int i = 0; i != n_species; ++i)
        {
            beta0s_P.col(i) = beta0s_mu[i] * beta0s[i].row(j).t()
                              + sqrt(sigma2(j,i)) * beta0s_S_U[i].t() * arma::randn<arma::vec>(n_P);
        }
        beta0s_P.col(other) = beta0s_mu[other] * beta0s[other].row(j).t()
                              + sqrt(eta2(j)) * beta0s_S_U[other].t() * arma::randn<arma::vec>(n_P);


        beta0s_1_P.col(j) = beta0s_P.col(0);
        beta0s_2_P.col(j) = beta0s_P.col(1);
        beta0s_3_P.col(j) = beta0s_P.col(2);
        beta0s_4_P.col(j) = beta0s_P.col(3);
        beta0s_5_P.col(j) = beta0s_P.col(4);
        beta0s_o_P.col(j) = beta0s_P.col(5);

        
        ////
        // Predict Z_i
        ////

        for(int i = 0; i != n_species+1; ++i)
        {
            arma::vec Z_i_P;
            if (trans == "log")
                Z_i_P = arma::exp(beta0(j,i) + beta0s_P.col(i) + beta1(j,i) * log_Q_P.col(i));
            else if (trans == "max")
                Z_i_P = zero_trunc( beta0(j,i) + beta0s_P.col(i) + beta1(j,i) * Q_P.col(i) );
            else
                Z_i_P = beta0(j,i) + beta0s_P.col(i) + beta1(j,i) * Q_P.col(i);
        

            if (i==0) Z_1_P.col(j) = Z_i_P;
            if (i==1) Z_2_P.col(j) = Z_i_P;
            if (i==2) Z_3_P.col(j) = Z_i_P;
            if (i==3) Z_4_P.col(j) = Z_i_P;
            if (i==4) Z_5_P.col(j) = Z_i_P;
            if (i==5) Z_o_P.col(j) = Z_i_P;
        }
    }
    
    Z_tot_P = Z_o_P + Z_1_P + Z_2_P + Z_3_P + Z_4_P + Z_5_P;


    return Rcpp::List::create(
                Rcpp::Named("Z_1")   = Z_1_P,
                Rcpp::Named("Z_2")   = Z_2_P,
                Rcpp::Named("Z_3")   = Z_3_P,
                Rcpp::Named("Z_4")   = Z_4_P,
                Rcpp::Named("Z_5")   = Z_5_P,
                Rcpp::Named("Z_o")   = Z_o_P,
                Rcpp::Named("Z_tot") = Z_tot_P,
                Rcpp::Named("beta0s_1") = beta0s_1_P,
                Rcpp::Named("beta0s_2") = beta0s_2_P,
                Rcpp::Named("beta0s_3") = beta0s_3_P,
                Rcpp::Named("beta0s_4") = beta0s_4_P,
                Rcpp::Named("beta0s_5") = beta0s_5_P,
                Rcpp::Named("beta0s_o") = beta0s_o_P 
           );      
}

