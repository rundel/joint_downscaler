// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(local_include)]]


#include <RcppArmadillo.h>

#include <assert.hpp>
#include <adapt_mcmc.hpp>
#include <distance.hpp>


arma::mat zero_trunc(arma::mat const& v)
{
    return (v+arma::abs(v))/2.0;
}



// [[Rcpp::export]]
Rcpp::List downscaler_joint(arma::mat C, arma::mat I, arma::mat F,
                            arma::mat Q_C, arma::mat Q_I, arma::mat Q_F, arma::mat Q_P,
                            arma::mat s_C, arma::mat s_I, arma::mat s_F, arma::mat s_P,
                            arma::vec phi, double xi, 
                            arma::vec beta0, arma::vec beta1, 
                            int niter, int nthin, int nburn,
                            bool predict = false,
                            bool verbose = false,
                            bool profile = false,
                            bool fix_phi = true,
                            bool fix_beta1 = false,
                            int bisp1 = 0, int bisp2 = 0
                           )
{
    bool debug = false;

    if(debug)
    {
        niter=10;
        nburn=10;
        nthin=1;
    }

    int nr_C = s_C.n_rows;
    int nr_I = s_I.n_rows;
    int nr_F = s_F.n_rows;
    int nr_P = s_P.n_rows;

    int nr_CI = nr_C + nr_I;
    int nr_CIF = nr_C + nr_I + nr_F;

    int n_species = 5;
    int tot_index = 5;
    int other_index = 5;

    RT_ASSERT(C.n_cols == n_species+1, "Dimension mismatch with CSN.");
    RT_ASSERT(I.n_cols == n_species+1, "Dimension mismatch with IMPROVE.");
    RT_ASSERT(phi.n_elem == n_species, "Dimension mismatch with phi.");
    
    bool use_bivariate = bisp1 && bisp2;

    if (use_bivariate)
    {
        bisp1--;
        bisp2--;

        RT_ASSERT(bisp1 >= 0 && bisp1 < n_species, "Bivariate species 1 out of range.");
        RT_ASSERT(bisp2 >= 0 && bisp2 < n_species, "Bivariate species 2 out of range.");
    }

    // Data

    arma::vec C_tot = C.col(tot_index);
    arma::vec I_tot = I.col(tot_index);
    arma::vec F_tot = F.col(0);
    
    arma::vec CI_tot = arma::join_cols(C_tot,I_tot);
    arma::vec CIF_tot = arma::join_cols(CI_tot,F_tot);

    arma::mat CI = arma::join_cols(C,I);

    arma::mat Q_CI = arma::join_cols(Q_C, Q_I);
    arma::mat Q = arma::join_cols(Q_CI, Q_F);

    // Distance Matrices

    arma::mat s_CI = arma::join_cols(s_C, s_I);
    arma::mat s_CIF = arma::join_cols(s_CI, s_F);

    arma::mat d_CI = circle_dist( s_CI );
    
    arma::mat d_CI_F = circle_dist( s_CI, s_F );
    arma::mat d_F = circle_dist( s_F );

    arma::mat d_CIF = circle_dist( s_CIF );

    arma::mat d_CIF_P = circle_dist( s_CIF, s_P );
    arma::mat d_P = circle_dist( s_P );

    double max_d = d_CIF.max();

    // Priors
    double b2_0 = 500; // prior variance for \beta_0
    double b2_1 = 500; // prior variance for \beta_1

    double c2_A11 = 500;
    double c2_A22 = 500;
    double c2_A21 = 500;

    arma::mat tau2_i_b = 2.0 * arma::ones<arma::mat>(2, n_species);
    double tau2_i_a = 2.0;

    arma::vec tau2_tot_b = 2.0 * arma::ones<arma::vec>(3);
    double tau2_tot_a = 2.0;    

    
    double eta2_a = 2;
    double eta2_b = 2;

    double sigma2_a = 2.0;
    double sigma2_b = 2.0;
    

    // Model state

    if (fix_beta1)
    {
        for(int i = 0; i != n_species; ++i)
        {
            beta1(i) = 0;
        }
    }

    arma::vec w0, w1;

    if(use_bivariate)
    {
        w0 = arma::zeros<arma::vec>(nr_CIF);
        w1 = arma::zeros<arma::vec>(nr_CIF);
    }

    double A11 = 1.0, 
           A21 = 0.0,
           A22 = 1.0;

    arma::mat beta0s = arma::zeros<arma::mat>(nr_CIF, n_species+1);
    
    if (use_bivariate)
    {
        beta0s.col(bisp1) = A11 * w0;
        beta0s.col(bisp2) = A21 * w0 + A22 * w1;
    }

    arma::mat Z_i(nr_CIF,n_species);
    arma::mat Z_i_star(nr_CIF,n_species);

    for(int i = 0; i != n_species; ++i)
    {   
        Z_i_star.col(i) = beta0(i) + beta0s.col(i) + beta1(i) * Q.col(i);
    }  
    Z_i = zero_trunc(Z_i_star); 

    arma::vec Z_o_star = beta0(other_index) + beta0s.col(other_index) + beta1(other_index) * Q.col(other_index);
    arma::vec Z_o = zero_trunc(Z_o_star); 

    arma::mat tau2_i = arma::ones<arma::mat>(2,n_species);
    arma::mat tau2_tot = arma::ones<arma::vec>(3);

    arma::mat tau2_i_vec(nr_CI,n_species);
    for(int i = 0; i != n_species; ++i)
    {   
        tau2_i_vec.col(i) = arma::join_cols( arma::ones<arma::vec>(nr_C) * tau2_i(0,i),
                                             arma::ones<arma::vec>(nr_I) * tau2_i(1,i) );
    }

    arma::vec tau2_tot_vec = arma::join_cols( arma::ones<arma::vec>(nr_C) * tau2_tot(0),
                             arma::join_cols( arma::ones<arma::vec>(nr_I) * tau2_tot(1),
                                              arma::ones<arma::vec>(nr_F) * tau2_tot(2) ));


    std::vector<arma::mat> Sigma_CIF_inv(n_species);
    std::vector<double> Sigma_CIF_log_det(n_species);
    
    for(int i = 0; i != n_species; ++i)
    {
        arma::mat Sigma_CIF = arma::exp(-phi(i) * d_CIF);
        arma::mat Sigma_CIF_U = arma::chol(Sigma_CIF);
        arma::mat Sigma_CIF_U_inv = arma::inv(arma::trimatu(Sigma_CIF_U));
        
        Sigma_CIF_inv[i] = Sigma_CIF_U_inv * Sigma_CIF_U_inv.t();
        Sigma_CIF_log_det[i] = 2*arma::accu(arma::log(Sigma_CIF_U.diag()));
    }
    
    arma::vec sigma2 = 1.0 * arma::ones<arma::vec>(n_species);

    arma::mat Sigma_o_CIF = arma::exp(-xi * d_CIF);
    
    arma::mat Sigma_o_CIF_U = arma::chol(Sigma_o_CIF);
    arma::mat Sigma_o_CIF_U_inv = arma::inv(arma::trimatu(Sigma_o_CIF_U));
    
    arma::mat Sigma_o_CIF_inv = Sigma_o_CIF_U_inv * Sigma_o_CIF_U_inv.t();
    double Sigma_o_CIF_log_det = 2*arma::accu(arma::log(Sigma_o_CIF_U.diag()));

    double eta2 = 1.0;



    // MCMC Adaptation

    arma::vec beta0_tuning = 0.01 * arma::ones<arma::vec>(n_species+1);
    vihola_ind_adapt beta0_amcmc(nburn*nthin, 0.45, 0.5, beta0_tuning);
    arma::uvec accept_beta0 = arma::zeros<arma::uvec>(n_species+1);

    arma::vec beta1_tuning = 0.01 * arma::ones<arma::vec>(n_species+1);
    vihola_ind_adapt beta1_amcmc(nburn*nthin, 0.45, 0.5, beta1_tuning);
    arma::uvec accept_beta1 = arma::zeros<arma::uvec>(n_species+1);

    arma::vec beta0s_tuning = 0.01 * arma::ones<arma::vec>(nr_CIF);
    std::vector<vihola_ind_adapt> beta0s_amcmc(n_species+1, vihola_ind_adapt(nburn*nthin, 0.45, 0.5, beta0s_tuning));
    arma::umat accept_beta0s = arma::zeros<arma::umat>(nr_CIF, n_species+1);
    
    arma::vec A_tuning = 0.1 * arma::ones<arma::vec>(3);
    vihola_ind_adapt A_amcmc(nburn*nthin, 0.45, 0.5, A_tuning);
    arma::uvec accept_A = arma::zeros<arma::uvec>(3);

    arma::vec phi_tuning = 0.001 * arma::ones<arma::vec>(n_species+1);
    vihola_ind_adapt phi_amcmc(nburn*nthin, 0.45, 0.5, phi_tuning);
    arma::uvec accept_phi = arma::zeros<arma::uvec>(n_species+1);

    // Prediction
    arma::mat beta0s_P = arma::zeros<arma::mat>(nr_P, n_species+1);

    arma::vec w0_P, w1_P;

    if(use_bivariate)
    {
        w0_P.resize(nr_P);
        w1_P.resize(nr_P);
    }

    arma::mat Z_i_star_P = arma::mat(nr_P,n_species);
    arma::mat Z_i_P = arma::mat(nr_P,n_species);
    arma::vec Z_o_P = arma::zeros<arma::vec>(nr_P); 
    arma::vec Z_tot_P = arma::zeros<arma::vec>(nr_P); 

    std::vector<arma::mat> Sigma_P(n_species);
    std::vector<arma::mat> Sigma_CIF_P(n_species);
    
    for(int i = 0; i != n_species; ++i)
    {
        Sigma_P[i] = arma::exp(-phi(i) * d_P);
        Sigma_CIF_P[i] = arma::exp(-phi(i) * d_CIF_P);
    }

    arma::mat Sigma_o_P = arma::exp(-xi * d_P);
    arma::mat Sigma_o_CIF_P = arma::exp(-xi * d_CIF_P);


    // Posterior Samples

    int n_post = (niter == nburn) ? niter : niter-nburn;

    arma::mat post_beta0(n_post, n_species+1);
    arma::mat post_beta1(n_post, n_species+1);

    arma::mat post_beta0s_1(n_post, nr_CIF);
    arma::mat post_beta0s_2(n_post, nr_CIF);
    arma::mat post_beta0s_3(n_post, nr_CIF);
    arma::mat post_beta0s_4(n_post, nr_CIF);
    arma::mat post_beta0s_5(n_post, nr_CIF);
    arma::mat post_beta0s_6(n_post, nr_CIF);

    arma::mat post_w0, post_w1;
    arma::mat post_A;
    arma::mat post_w0_P, post_w1_P;

    
    if(use_bivariate)
    {
        post_w0.resize(n_post, nr_CIF);
        post_w1.resize(n_post, nr_CIF);

        post_A.resize(n_post, 3);

        post_w0_P.resize(n_post, nr_P);
        post_w1_P.resize(n_post, nr_P);        
    }

    arma::mat post_Z_o(n_post, nr_CIF);
    arma::mat post_Z_1(n_post, nr_CIF);
    arma::mat post_Z_2(n_post, nr_CIF);
    arma::mat post_Z_3(n_post, nr_CIF);
    arma::mat post_Z_4(n_post, nr_CIF);
    arma::mat post_Z_5(n_post, nr_CIF);

    arma::mat post_tau2_1(n_post, 2);
    arma::mat post_tau2_2(n_post, 2);
    arma::mat post_tau2_3(n_post, 2);
    arma::mat post_tau2_4(n_post, 2);
    arma::mat post_tau2_5(n_post, 2);

    arma::mat post_tau2_tot(n_post, 3);

    arma::mat post_sigma2(n_post, 5);

    arma::vec post_eta2(n_post);

    arma::mat post_phi(n_post, n_species);
    arma::mat post_xi(n_post, 1);

    arma::mat post_beta0s_1_P(n_post, nr_P);
    arma::mat post_beta0s_2_P(n_post, nr_P);
    arma::mat post_beta0s_3_P(n_post, nr_P);
    arma::mat post_beta0s_4_P(n_post, nr_P);
    arma::mat post_beta0s_5_P(n_post, nr_P);
    arma::mat post_beta0s_6_P(n_post, nr_P);


    arma::mat post_Z_o_P(n_post, nr_P);
    arma::mat post_Z_1_P(n_post, nr_P);
    arma::mat post_Z_2_P(n_post, nr_P);
    arma::mat post_Z_3_P(n_post, nr_P);
    arma::mat post_Z_4_P(n_post, nr_P);
    arma::mat post_Z_5_P(n_post, nr_P);

    arma::mat post_Z_tot_P(n_post, nr_P);


    // MCMC

    for(int iter = 0; iter != niter; ++iter)
    {
        for(int thin = 0; thin != nthin; ++thin)
        {
            ////
            // Update \beta_0
            ////

            {
                arma::vec jump = beta0_amcmc.get_jump();
                arma::vec alpha(n_species+1);


                // 5 Species
                for(int i = 0; i != n_species; ++i)
                {   
                    arma::vec n_D = CIF_tot - Z_o;

                    for(int j = 0; j != n_species; ++j)
                    {
                        if(i != j)
                            n_D -= Z_i.col(j);
                    }

                    double beta0_prop = beta0(i) + jump(i);

                    arma::vec Z_i_star_prop = beta0_prop + beta0s.col(i) + beta1(i) * Q.col(i);
                    arma::vec Z_i_prop = zero_trunc(Z_i_star_prop);

                    double loglik_cur = -0.5 * pow(beta0(i),2) / b2_0
                                        -0.5 * arma::accu(arma::square(CI.col(i) - Z_i.submat(0,i,nr_CI-1,i)) / tau2_i_vec.col(i))
                                        -0.5 * arma::accu(arma::square(n_D - Z_i.col(i)) / tau2_tot_vec);
                    

                    double loglik_prop = -0.5 * pow(beta0_prop,2) / b2_0
                                         -0.5 * arma::accu(arma::square(CI.col(i) - Z_i_prop.subvec(0,nr_CI-1)) / tau2_i_vec.col(i))
                                         -0.5 * arma::accu(arma::square(n_D - Z_i_prop) / tau2_tot_vec);

                    alpha(i) = std::min(1.0, exp(loglik_prop-loglik_cur));

                    if (Rcpp::runif(1)[0] <= alpha(i))
                    {
                        beta0(i) = beta0_prop;
                        
                        Z_i_star.col(i) = Z_i_star_prop;
                        Z_i.col(i) = Z_i_prop; 
                        
                        accept_beta0(i)++;
                    }
                }

                // Other
                arma::vec n_D = CIF_tot;

                for(int j = 0; j != n_species; ++j)
                {
                    n_D -= Z_i.col(j);
                }

                double beta0_prop = beta0(other_index) + jump(other_index);

                arma::vec Z_o_star_prop = beta0_prop + beta0s.col(other_index) + beta1(other_index) * Q.col(other_index);
                arma::vec Z_o_prop = zero_trunc(Z_o_star_prop);

                double loglik_cur = -0.5 * pow(beta0(other_index),2) / b2_0
                                    -0.5 * arma::accu(arma::square(n_D - Z_o) / tau2_tot_vec);
                

                double loglik_prop = -0.5 * pow(beta0_prop,2) / b2_0
                                     -0.5 * arma::accu(arma::square(n_D - Z_o_prop) / tau2_tot_vec);

                alpha(other_index) = std::min(1.0, exp(loglik_prop-loglik_cur));

                if (Rcpp::runif(1)[0] <= alpha(other_index))
                {
                    beta0(other_index) = beta0_prop;
                    
                    Z_o_star = Z_o_star_prop;
                    Z_o      = Z_o_prop; 
                    
                    accept_beta0(other_index)++;
                }


                beta0_amcmc.update(iter*nthin+thin, alpha);
            }


            ////
            // Update w0(s)
            ////

            
            if (use_bivariate)
            {
                arma::vec n = CIF_tot - Z_o;

                for(int j = 0; j != n_species; ++j)
                {
                    if(j != bisp1 && j != bisp2) 
                        n -= Z_i.col(j);
                }

                arma::vec jump = beta0s_amcmc[bisp1].get_jump();
                arma::vec alpha(nr_CIF);
                
                for(int j = 0; j != nr_CIF; ++j)
                {
                    double w0_prop = w0(j) + jump(j);
                    
                    double beta0s_sp1_prop = A11 * w0_prop;
                    double beta0s_sp2_prop = A21 * w0_prop + A22 * w1(j);
                    
                    double Z_i_sp1_star_prop = beta0(bisp1) + beta0s_sp1_prop + beta1(bisp1) * Q(j,bisp1);
                    double Z_i_sp2_star_prop = beta0(bisp2) + beta0s_sp2_prop + beta1(bisp2) * Q(j,bisp2);
                    
                    double Z_i_sp1_prop = std::max(0.0, Z_i_sp1_star_prop);
                    double Z_i_sp2_prop = std::max(0.0, Z_i_sp2_star_prop);

                    double loglik_ratio = -0.5 * (   jump(j) * jump(j) * Sigma_CIF_inv[bisp1](j,j)
                                                   + 2 * jump(j) * arma::as_scalar(Sigma_CIF_inv[bisp1].row(j) * w0)
                                                   + pow(n(j) - Z_i_sp1_prop - Z_i_sp2_prop, 2) / tau2_tot_vec(j)
                                                   - pow(n(j) - Z_i(j,bisp1) - Z_i(j,bisp2), 2) / tau2_tot_vec(j)
                                                 );

                    if (j < nr_CI)
                        loglik_ratio += -0.5 * ( pow(CI(j,bisp1)-Z_i_sp1_prop,2) - pow(CI(j,bisp1)-Z_i(j,bisp1),2) ) / tau2_i_vec(j,bisp1)
                                        -0.5 * ( pow(CI(j,bisp2)-Z_i_sp2_prop,2) - pow(CI(j,bisp2)-Z_i(j,bisp2),2) ) / tau2_i_vec(j,bisp2);


                    if (debug)
                    {
                        arma::vec w0_prop_vec = w0;
                        w0_prop_vec(j) += jump(j);

                        arma::vec Z_i_sp1_prop_vec = zero_trunc(beta0(bisp1) + A11 * w0_prop_vec + beta1(bisp1) * Q.col(bisp1));
                        arma::vec Z_i_sp2_prop_vec = zero_trunc(beta0(bisp2) + A21 * w0_prop_vec + A22 * w1 + beta1(bisp2) * Q.col(bisp2));

                        double logpost_cur = arma::as_scalar(
                                                  -0.5 * w0.t() * Sigma_CIF_inv[bisp1] * w0
                                                + -0.5 * arma::square(CI.col(bisp1) - Z_i.submat(0,bisp1,nr_CI-1,bisp1)).t() * (1/tau2_i_vec.col(bisp1))
                                                + -0.5 * arma::square(CI.col(bisp2) - Z_i.submat(0,bisp2,nr_CI-1,bisp2)).t() * (1/tau2_i_vec.col(bisp2))
                                                + -0.5 * arma::square(n - Z_i.col(bisp1) - Z_i.col(bisp2)).t() * (1/tau2_tot_vec)
                                             );

                        double logpost_prop = arma::as_scalar(
                                                  -0.5 * w0_prop_vec.t() * Sigma_CIF_inv[bisp1] * w0_prop_vec
                                                + -0.5 * arma::square(CI.col(bisp1) - Z_i_sp1_prop_vec.rows(0,nr_CI-1)).t() * (1/tau2_i_vec.col(bisp1)) 
                                                + -0.5 * arma::square(CI.col(bisp2) - Z_i_sp2_prop_vec.rows(0,nr_CI-1)).t() * (1/tau2_i_vec.col(bisp2)) 
                                                + -0.5 * arma::square( n - Z_i_sp1_prop_vec - Z_i_sp2_prop_vec).t() * (1/tau2_tot_vec)
                                              );


                        if (abs((logpost_prop - logpost_cur) - loglik_ratio) > 1e-6)
                        {   
                            Rcpp::Rcout << "Short calc: " << loglik_ratio << "\n";
                            Rcpp::Rcout << "Long  calc: " << logpost_prop - logpost_cur << "\n";

                            RT_ASSERT(false, "Metropolis ratio error (w0).");
                        }
                    }

                    alpha(j) = std::min(1.0, exp(loglik_ratio));

                    if (Rcpp::runif(1)[0] <= alpha(j))
                    {
                        w0(j) = w0_prop;

                        beta0s(j,bisp1) = beta0s_sp1_prop;
                        beta0s(j,bisp2) = beta0s_sp2_prop;
                        
                        Z_i_star(j,bisp1) = Z_i_sp1_star_prop;
                        Z_i_star(j,bisp2) = Z_i_sp2_star_prop;
                        
                        Z_i(j,bisp1) = Z_i_sp1_prop;
                        Z_i(j,bisp2) = Z_i_sp2_prop;

                        accept_beta0s(j,bisp1)++;
                    }
                }

                beta0s_amcmc[bisp1].update(iter*nthin+thin, alpha);
            }

            ////
            // Update w1(s)
            ////

            if (use_bivariate)
            {
                arma::vec n = CIF_tot - Z_o;

                for(int j = 0; j != n_species; ++j)
                {
                    if(j != bisp2) 
                        n -= Z_i.col(j);
                }

                arma::vec jump = beta0s_amcmc[bisp2].get_jump();
                arma::vec alpha(nr_CIF);
                
                for(int j = 0; j != nr_CIF; ++j)
                {
                    double w1_prop = w1(j) + jump(j);
                    double beta0s_prop = A21 * w0(j) + A22 * w1_prop;
                    double Z_i_star_prop = beta0(bisp2) + beta0s_prop + beta1(bisp2) * Q(j,bisp2);
                    double Z_i_prop = std::max(0.0, Z_i_star_prop);

                    double loglik_ratio = -0.5 * (   jump(j) * jump(j) * Sigma_CIF_inv[bisp2](j,j)
                                                   + 2 * jump(j) * arma::as_scalar(Sigma_CIF_inv[bisp2].row(j) * w1)
                                                   + pow(n(j) - Z_i_prop,2) / tau2_tot_vec(j)
                                                   - pow(n(j) - Z_i(j,bisp2),2) / tau2_tot_vec(j)
                                                 );

                    if (j < nr_CI)
                        loglik_ratio += -0.5 * (  pow(CI(j,bisp2)-Z_i_prop,2)
                                                - pow(CI(j,bisp2)-Z_i(j,bisp2),2) ) / tau2_i_vec(j,bisp2);


                    if (debug)
                    {
                        double logpost_cur = arma::as_scalar(
                                                  -0.5 * w1.t() * Sigma_CIF_inv[bisp2] * w1
                                                + -0.5 * arma::square(CI.col(bisp2) - Z_i.submat(0,bisp2,nr_CI-1,bisp2)).t() * (1/tau2_i_vec.col(bisp2))
                                                + -0.5 * arma::square(n - Z_i.col(bisp2)).t() * (1/tau2_tot_vec)
                                             );

                        arma::vec w1_prop_vec = w1;
                        w1_prop_vec(j) += jump(j);
                        
                        arma::vec Z_i_prop_vec = zero_trunc(beta0(bisp2) + A21 * w0 + A22 * w1_prop_vec + beta1(bisp2) * Q.col(bisp2));


                        double logpost_prop = arma::as_scalar(
                                                  -0.5 * w1_prop_vec.t() * Sigma_CIF_inv[bisp2] * w1_prop_vec
                                                + -0.5 * arma::square(CI.col(bisp2) - Z_i_prop_vec.rows(0,nr_CI-1)).t() * (1/tau2_i_vec.col(bisp2))
                                                + -0.5 * arma::square(n - Z_i_prop_vec).t() * (1/tau2_tot_vec)
                                              );


                        if (abs((logpost_prop - logpost_cur) - loglik_ratio) > 1e-6)
                        {   
                            Rcpp::Rcout << "Short calc: " << loglik_ratio << "\n";
                            Rcpp::Rcout << "Long  calc: " << logpost_prop - logpost_cur << "\n";

                            RT_ASSERT(false, "Metropolis ratio error (w1).");
                        }
                    }

                    alpha(j) = std::min(1.0, exp(loglik_ratio));

                    if (Rcpp::runif(1)[0] <= alpha(j))
                    {
                        w1(j) = w1_prop;
                        beta0s(j,bisp2) = beta0s_prop;
                        Z_i_star(j,bisp2) = Z_i_star_prop;
                        Z_i(j,bisp2) = Z_i_prop;

                        accept_beta0s(j,bisp2)++;
                    }
                }

                beta0s_amcmc[bisp2].update(iter*nthin+thin, alpha);
            }
            

            ////
            // Update \beta_0(s)
            ////


            // 5 Species
            for(int i = 0; i != n_species; ++i)
            {
                if (use_bivariate && (i == bisp1 || i == bisp2))
                    continue;

                arma::vec n = CIF_tot - Z_o;

                for(int j = 0; j != n_species; ++j)
                {
                    if(i != j) 
                        n -= Z_i.col(j);
                }

                arma::mat C = Sigma_CIF_inv[i] / sigma2(i);

                arma::vec jump = beta0s_amcmc[i].get_jump();
                arma::vec alpha(nr_CIF);
                
                for(int j = 0; j != nr_CIF; ++j)
                {
                    double beta0s_prop = beta0s(j,i) + jump(j);
                    double Z_i_star_prop = beta0(i) + beta0s_prop + beta1(i) * Q(j,i);
                    double Z_i_prop = std::max(0.0, Z_i_star_prop);

                    double loglik_ratio = -0.5 * (   jump(j) * jump(j) * C(j,j)
                                                   + 2 * jump(j) * arma::as_scalar(C.row(j) * beta0s.col(i))
                                                   + pow(n(j) - Z_i_prop,2) / tau2_tot_vec(j)
                                                   - pow(n(j) - Z_i(j,i),2) / tau2_tot_vec(j)
                                                 );

                    if (j < nr_CI)
                        loglik_ratio += -0.5 * (  pow(CI(j,i)-Z_i_prop,2)
                                                - pow(CI(j,i)-Z_i(j,i),2) ) / tau2_i_vec(j,i);

                    alpha(j) = std::min(1.0, exp(loglik_ratio));

                    if (Rcpp::runif(1)[0] <= alpha(j))
                    {
                        beta0s(j,i) = beta0s_prop;
                        
                        Z_i_star(j,i) = Z_i_star_prop;
                        Z_i(j,i) = Z_i_prop;

                        accept_beta0s(j,i)++;
                    }
                }

                beta0s_amcmc[i].update(iter*nthin+thin, alpha);
            }

            // Other
            {
                arma::vec n = CIF_tot;

                for(int j = 0; j != n_species; ++j)
                {
                    n -= Z_i.col(j);
                }

                arma::mat C = Sigma_o_CIF_inv / eta2;

                arma::vec jump = beta0s_amcmc[other_index].get_jump();
                arma::vec alpha(nr_CIF);
                
                for(int j = 0; j != nr_CIF; ++j)
                {
                    double beta0s_prop = beta0s(j,other_index) + jump(j);
                    
                    double Z_o_star_prop = beta0(other_index) + beta0s_prop + beta1(other_index) * Q(j,other_index);
                    double Z_o_prop = std::max(0.0, Z_o_star_prop);

                    double loglik_ratio = -0.5 * (   jump(j) * jump(j) * C(j,j)
                                                   + 2 * jump(j) * arma::as_scalar(C.row(j) * beta0s.col(other_index))
                                                   + pow(n(j) - Z_o_prop,2) / tau2_tot_vec(j)
                                                   - pow(n(j) - Z_o(j),2) / tau2_tot_vec(j)
                                                 );

                    alpha(j) = std::min(1.0, exp(loglik_ratio));

                    if (Rcpp::runif(1)[0] <= alpha(j))
                    {
                        beta0s(j,other_index) = beta0s_prop;
                        
                        Z_o_star(j) = Z_o_star_prop;
                        Z_o(j) = Z_o_prop;

                        accept_beta0s(j,other_index)++;
                    }
                }

                beta0s_amcmc[other_index].update(iter*nthin+thin, alpha);
            }


            ////
            // Update A_11, A_21, A_22
            ////

            if (use_bivariate)
            {
                arma::vec alpha(3);
                arma::vec jump = A_amcmc.get_jump();

                double A11_prop = A11 + jump(0);
                double A21_prop = A21 + jump(1);
                double A22_prop = A22 + jump(2);

                if (A11_prop > 0)
                {
                    arma::vec n = CIF_tot - Z_o;

                    for(int i = 0; i != n_species; ++i)
                    {
                        if (i != bisp1)
                            n -= Z_i.col(i);
                    }

                    arma::vec Z_i_prop = zero_trunc(beta0(bisp1) + A11_prop * w0 + beta1(bisp1) * Q.col(bisp1));

                    double log_post_cur = -0.5 * arma::accu( arma::square(CI.col(bisp1) - Z_i.submat(0,bisp1,nr_CI-1,bisp1)) / tau2_i_vec.col(bisp1))
                                          -0.5 * arma::accu( arma::square(n - Z_i.col(bisp1)) / tau2_tot_vec )
                                          -0.5 * pow(log(A11),2) / c2_A11 - log(A11);

                    double log_post_prop = -0.5 * arma::accu( arma::square(CI.col(bisp1) - Z_i_prop.rows(0,nr_CI-1)) / tau2_i_vec.col(bisp1))
                                           -0.5 * arma::accu( arma::square(n - Z_i_prop) / tau2_tot_vec )
                                           -0.5 * pow(log(A11_prop),2) / c2_A11 - log(A11_prop);

                    alpha(0) = std::min(1.0, exp(log_post_prop - log_post_cur));

                    if (Rcpp::runif(1)[0] <= alpha(0))
                    {
                        A11 = A11_prop;
                        accept_A(0)++;

                        arma::vec old_beta0s = beta0s.col(bisp1);
                        beta0s.col(bisp1) = A11 * w0;
                        Z_i_star.col(bisp1) += beta0s.col(bisp1) - old_beta0s;
                        Z_i.col(bisp1) = zero_trunc(Z_i_star.col(bisp1));
                    }
                }
                else
                {
                    alpha(0) = 0.0;
                }


                if (A21_prop > 0)
                {
                    arma::vec n = CIF_tot - Z_o;

                    for(int i = 0; i != n_species; ++i)
                    {
                        if (i != bisp2)
                            n -= Z_i.col(i);
                    }

                    arma::vec Z_i_prop = zero_trunc(beta0(bisp2) + A21_prop * w0 + A22 * w1 + beta1(bisp2) * Q.col(bisp2));

                    double log_post_cur = -0.5 * arma::accu( arma::square(CI.col(bisp2) - Z_i.submat(0,bisp2,nr_CI-1,bisp2)) / tau2_i_vec.col(bisp2))
                                          -0.5 * arma::accu( arma::square(n - Z_i.col(bisp2)) / tau2_tot_vec )
                                          -0.5 * pow(A21,2) / c2_A21;

                    double log_post_prop = -0.5 * arma::accu( arma::square(CI.col(bisp2) - Z_i_prop.rows(0,nr_CI-1)) / tau2_i_vec.col(bisp2))
                                           -0.5 * arma::accu( arma::square(n - Z_i_prop) / tau2_tot_vec )
                                           -0.5 * pow(A21_prop,2) / c2_A21;



                    alpha(1) = std::min(1.0, exp(log_post_prop - log_post_cur));

                    if (Rcpp::runif(1)[0] <= alpha(1))
                    {
                        A21 = A21_prop;
                        accept_A(1)++;

                        arma::vec old_beta0s = beta0s.col(bisp2);
                        beta0s.col(bisp2) = A21 * w0 + A22 * w1;
                        Z_i_star.col(bisp2) += beta0s.col(bisp2) - old_beta0s;
                        Z_i.col(bisp2) = zero_trunc(Z_i_star.col(bisp2));
                    }
                }
                else
                {
                    alpha(1) = 0.0;
                }

                if (A22_prop > 0)
                {
                    arma::vec n = CIF_tot - Z_o;

                    for(int i = 0; i != n_species; ++i)
                    {
                        if (i != bisp2)
                            n -= Z_i.col(i);
                    }

                    arma::vec Z_i_prop = zero_trunc(beta0(bisp2) + A21 * w0 + A22_prop * w1 + beta1(bisp2) * Q.col(bisp2));

                    double log_post_cur = -0.5 * arma::accu( arma::square(CI.col(bisp2) - Z_i.submat(0,bisp2,nr_CI-1,bisp2)) / tau2_i_vec.col(bisp2))
                                          -0.5 * arma::accu( arma::square(n - Z_i.col(bisp2)) / tau2_tot_vec )
                                          -0.5 * pow(log(A22),2) / c2_A22 - log(A22);

                    double log_post_prop = -0.5 * arma::accu( arma::square(CI.col(bisp2) - Z_i_prop.rows(0,nr_CI-1)) / tau2_i_vec.col(bisp2))
                                           -0.5 * arma::accu( arma::square(n - Z_i_prop) / tau2_tot_vec )
                                           -0.5 * pow(log(A22_prop),2) / c2_A22 - log(A22_prop);



                    alpha(2) = std::min(1.0, exp(log_post_prop - log_post_cur));

                    if (Rcpp::runif(1)[0] <= alpha(2))
                    {
                        A22 = A22_prop;
                        accept_A(2)++;

                        arma::vec old_beta0s = beta0s.col(bisp2);
                        beta0s.col(bisp2) = A21 * w0 + A22 * w1;
                        Z_i_star.col(bisp2) += beta0s.col(bisp2) - old_beta0s;
                        Z_i.col(bisp2) = zero_trunc(Z_i_star.col(bisp2));
                    }
                }
                else
                {
                    alpha(2) = 0.0;
                }

                A_amcmc.update(iter*nthin+thin, alpha);
            }

            ////
            // Update \beta_1
            ////


            if (!fix_beta1)
            {
                arma::vec jump = beta1_amcmc.get_jump();
                arma::vec alpha(n_species+1);

                // 5 Species
                for(int i = 0; i != n_species; ++i)
                {
                    arma::vec n_D = CIF_tot - Z_o;

                    for(int j = 0; j != n_species; ++j)
                    {
                        if(i != j)
                            n_D -= Z_i.col(j);
                    }

                    double beta1_prop = beta1(i) + jump(i);

                    arma::vec Z_i_star_prop = beta0(i) + beta0s.col(i) + beta1_prop * Q.col(i);
                    arma::vec Z_i_prop = zero_trunc(Z_i_star_prop);

                    double loglik_cur = -0.5 * pow(beta1(i),2) / b2_1
                                        -0.5 * arma::accu(arma::square(CI.col(i) - Z_i.submat(0,i,nr_CI-1,i)) / tau2_i_vec.col(i))
                                        -0.5 * arma::accu(arma::square(n_D - Z_i.col(i)) / tau2_tot_vec);
                    
                    double loglik_prop = -0.5 * pow(beta1_prop,2) / b2_1
                                         -0.5 * arma::accu(arma::square(CI.col(i) - Z_i_prop.subvec(0,nr_CI-1)) / tau2_i_vec.col(i))
                                         -0.5 * arma::accu(arma::square(n_D - Z_i_prop) / tau2_tot_vec);

                    alpha(i) = std::min(1.0, exp(loglik_prop-loglik_cur));
                    
                    if (Rcpp::runif(1)[0] <= alpha(i))
                    {
                        beta1(i) = beta1_prop;
                        
                        Z_i_star.col(i) = Z_i_star_prop; 
                        Z_i.col(i) = Z_i_prop; 
                        
                        accept_beta1(i)++;
                    }  
                }
                
                // Other

                arma::vec n_D = CIF_tot;

                for(int j = 0; j != n_species; ++j)
                {
                    n_D -= Z_i.col(j);
                }

                double beta1_prop = beta1(other_index) + jump(other_index);

                arma::vec Z_o_star_prop = beta0(other_index) + beta0s.col(other_index) + beta1_prop * Q.col(other_index);
                arma::vec Z_o_prop = zero_trunc(Z_o_star_prop);

                double loglik_cur = -0.5 * pow(beta1(other_index),2) / b2_1
                                    -0.5 * arma::accu(arma::square(n_D - Z_o) / tau2_tot_vec);
                
                double loglik_prop = -0.5 * pow(beta1_prop,2) / b2_1
                                     -0.5 * arma::accu(arma::square(n_D - Z_o_prop) / tau2_tot_vec);

                alpha(other_index) = std::min(1.0, exp(loglik_prop-loglik_cur));
                
                if (Rcpp::runif(1)[0] <= alpha(other_index))
                {
                    beta1(other_index) = beta1_prop;
                    
                    Z_o_star = Z_o_star_prop; 
                    Z_o = Z_o_prop; 
                    
                    accept_beta1(other_index)++;
                }

                beta1_amcmc.update(iter*nthin+thin, alpha);
            }            
            
            if (debug)
            {
                arma::mat tmp = Z_i;
                for(int i = 0; i != n_species; ++i)
                {       
                    tmp.col(i) = beta0(i) + beta0s.col(i) + beta1(i) * Q.col(i);
                }

                arma::mat diff = arma::abs(Z_i - zero_trunc(tmp));
                if (arma::accu(diff) > 1e-6)
                {
                    Rcpp::Rcout << arma::accu(diff) << "\n";
                    RT_ASSERT(false, "Debug failure - beta1 - Z_i status.");
                }

                if (use_bivariate)
                {
                    tmp.col(bisp1) += A11 * w0 - beta0s.col(bisp1);
                    tmp.col(bisp2) += A21 * w0 + A22 * w1 - beta0s.col(bisp2);
                
                    arma::mat diff = arma::abs(Z_i - zero_trunc(tmp));
                    if (arma::accu(diff) > 1e-6)
                    {
                        Rcpp::Rcout << arma::accu(diff) << "\n";
                        RT_ASSERT(false, "Debug failure - beta1 - Z_i status.");
                    }
                }
            }
            

            ////
            // Update \sigma^2
            ////
            for(int i = 0; i != n_species; ++i)
            {   
                if (use_bivariate && (i == bisp1 || i == bisp2))
                {
                    sigma2(i) = 0.0/0.0; // NaN
                    continue;
                }

                double a = sigma2_a + (nr_CIF)/2;
                double b = sigma2_b + arma::as_scalar(beta0s.col(i).t() * Sigma_CIF_inv[i] * beta0s.col(i)) / 2;

                sigma2(i) = 1.0 / ::Rf_rgamma(a, 1/b);
            }

            ////
            // Update \eta^2
            ////
            {   
                double a = eta2_a + (nr_CIF)/2;
                double b = eta2_b + arma::as_scalar( beta0s.col(other_index).t() * Sigma_o_CIF_inv * beta0s.col(other_index) ) / 2;

                eta2 = 1.0 / ::Rf_rgamma(a, 1/b);
            }

            ////
            // Update \phi & \xi
            ////

            if (!fix_phi)
            {
                arma::vec alpha(n_species+1);

                double phi_min = 3.0 / (max_d/2);
                double phi_max = 3.0 / 1.0;


                for(int i = 0; i != n_species; ++i)
                {
                    double prop_phi = phi(i) + phi_amcmc.get_jump()(i);

                    if (prop_phi >= phi_min && prop_phi <= phi_max)
                    {
                        arma::mat prop_Sigma = arma::exp(-prop_phi * d_CIF);
                        arma::mat prop_Sigma_U = arma::chol(prop_Sigma);
                        arma::mat prop_Sigma_U_inv = arma::inv(arma::trimatu(prop_Sigma_U));
                        arma::mat prop_Sigma_inv = prop_Sigma_U_inv * prop_Sigma_U_inv.t();
                        double prop_Sigma_log_det = 2*arma::accu(arma::log(prop_Sigma_U.diag()));

                        //Rcpp::Rcout << i << " : " << arma::accu(prop_Sigma_inv - arma::inv(arma::sympd(prop_Sigma))) << "\n";

                        double cur_log_prob, prop_log_prob;

                        
                        if (use_bivariate)
                        {
                            if (i == bisp1) 
                            {
                                cur_log_prob = - 0.5 * arma::as_scalar(w0.t() * Sigma_CIF_inv[i] * w0)
                                               - 0.5 * Sigma_CIF_log_det[i];

                                prop_log_prob = - 0.5 * arma::as_scalar(w0.t() * prop_Sigma_inv * w0)
                                                - 0.5 * prop_Sigma_log_det;
                            }
                            if (i == bisp2)
                            {
                                cur_log_prob = - 0.5 * arma::as_scalar(w1.t() * Sigma_CIF_inv[i] * w1)
                                               - 0.5 * Sigma_CIF_log_det[i];

                                prop_log_prob = - 0.5 * arma::as_scalar(w1.t() * prop_Sigma_inv * w1)
                                                - 0.5 * prop_Sigma_log_det;  
                            }
                        }
                        else
                        {
                            cur_log_prob = - 0.5 * arma::as_scalar(beta0s.col(i).t() * Sigma_CIF_inv[i] * beta0s.col(i)) / sigma2(i)
                                           - 0.5 * Sigma_CIF_log_det[i];

                            prop_log_prob = - 0.5 * arma::as_scalar(beta0s.col(i).t() * prop_Sigma_inv * beta0s.col(i)) / sigma2(i)
                                            - 0.5 * prop_Sigma_log_det;
                        }


                        alpha[i] = std::min(1.0, exp(prop_log_prob - cur_log_prob));
                                                
                        if (Rcpp::runif(1)[0] <= alpha[i])
                        {
                            phi(i) = prop_phi;
                            Sigma_CIF_inv[i] = prop_Sigma_inv;
                            Sigma_CIF_log_det[i] = prop_Sigma_log_det;
                            
                            accept_phi[i]++;
                        }
                    }
                    else
                    {
                        alpha[i] = 0.0;
                    }


                }


                double xi_min = phi_min;
                double xi_max = phi_max;

                double prop_xi = xi + phi_amcmc.get_jump()(other_index);

                if (prop_xi >= xi_min && prop_xi <= xi_max)
                {
                    double cur_log_prob = - 0.5 * arma::as_scalar(  beta0s.col(other_index).t() 
                                                                  * Sigma_o_CIF_inv 
                                                                  * beta0s.col(other_index)     ) / eta2
                                          - 0.5 * Sigma_o_CIF_log_det;


                    arma::mat prop_Sigma_o = arma::exp(-prop_xi * d_CIF);
                    arma::mat prop_Sigma_o_U = arma::chol(prop_Sigma_o);
                    arma::mat prop_Sigma_o_U_inv = arma::inv(arma::trimatu(prop_Sigma_o_U));
                    arma::mat prop_Sigma_o_inv = prop_Sigma_o_U_inv * prop_Sigma_o_U_inv.t();
                    double prop_Sigma_o_log_det = 2*arma::accu(arma::log(prop_Sigma_o_U.diag()));


                    double prop_log_prob = - 0.5 * arma::as_scalar(  beta0s.col(other_index).t() 
                                                                   * prop_Sigma_o_inv 
                                                                   * beta0s.col(other_index)    ) / eta2
                                           - 0.5 * prop_Sigma_o_log_det;

                    alpha[other_index] = std::min(1.0, exp(prop_log_prob - cur_log_prob));
                    
                    //Rcpp::Rcout << "Cur  : " << cur_log_prob << "\n"
                    //            << "Prop : " << prop_log_prob << "\n"
                    //            << "alpha: " << alpha[6] << "\n";

                    if (Rcpp::runif(1)[0] <= alpha[other_index])
                    {
                        xi = prop_xi;
                        Sigma_o_CIF_inv = prop_Sigma_o_inv;
                        Sigma_o_CIF_log_det = prop_Sigma_o_log_det;
                        
                        accept_phi[other_index]++;
                    }
                }
                else
                {
                    alpha[other_index] = 0.0;
                }

                //Rcpp::Rcout << alpha.t() << "\n";
                phi_amcmc.update(iter*nthin+thin,alpha);
            }


            ////
            // Update \tau^2_i
            ////
            for(int i = 0; i != n_species; ++i)
            {   
                arma::vec sq_diff = arma::square(CI.col(i) - Z_i.submat(0,i,nr_CI-1,i));

                double b0 = tau2_i_b(0,i) + 0.5 * arma::accu( sq_diff.subvec(0,nr_C-1) );
                double b1 = tau2_i_b(1,i) + 0.5 * arma::accu( sq_diff.subvec(nr_C,nr_CI-1) );
                
                double a0 = tau2_i_a + 0.5 * nr_C;
                double a1 = tau2_i_a + 0.5 * nr_I;

                
                tau2_i(0,i) = 1.0 / ::Rf_rgamma(a0, 1/b0);
                tau2_i(1,i) = 1.0 / ::Rf_rgamma(a1, 1/b1);


                tau2_i_vec.col(i) = arma::join_cols( arma::ones<arma::vec>(nr_C) * tau2_i(0,i),
                                                     arma::ones<arma::vec>(nr_I) * tau2_i(1,i) );
            }


            ////
            // Update tau^2_tot
            ////
            {   
                arma::vec m = CIF_tot - Z_o;
                
                for(int j = 0; j != n_species; ++j)
                {
                    m -= Z_i.col(j);
                }

                double b0 = tau2_tot_b(0) + 0.5 * arma::accu( arma::square(m.subvec(   0 , nr_C  -1)) );
                double b1 = tau2_tot_b(1) + 0.5 * arma::accu( arma::square(m.subvec(nr_C , nr_CI -1)) );
                double b2 = tau2_tot_b(2) + 0.5 * arma::accu( arma::square(m.subvec(nr_CI, nr_CIF-1)) );

                double a0 = tau2_tot_a + 0.5 * nr_C;
                double a1 = tau2_tot_a + 0.5 * nr_I;
                double a2 = tau2_tot_a + 0.5 * nr_F;
                    
                tau2_tot(0) = 1.0 / ::Rf_rgamma(a0, 1/b0);
                tau2_tot(1) = 1.0 / ::Rf_rgamma(a1, 1/b1);
                tau2_tot(2) = 1.0 / ::Rf_rgamma(a2, 1/b2);

                tau2_tot_vec = arma::join_cols( arma::ones<arma::vec>(nr_C) * tau2_tot(0),
                               arma::join_cols( arma::ones<arma::vec>(nr_I) * tau2_tot(1),
                                                arma::ones<arma::vec>(nr_F) * tau2_tot(2) ) );

            }


            ////
            // Output Status 
            ////
            int s = iter*nthin+thin;
            if ((verbose && (s+1) % ( std::max(niter*nthin/10,1) ) == 0) || iter+1 == niter)
            {
                Rcpp::Rcout << "beta0 accept rates       :";
                for(int i=0; i!=n_species+1; ++i)
                    Rcpp::Rcout << " " << ( (double) accept_beta0(i))  / ( (double) s + 1.0);
                Rcpp::Rcout << "\n";

                Rcpp::Rcout << "beta0s (avg) accept rates:";
                for(int i=0; i!=n_species+1; ++i)
                {
                    if (use_bivariate && (i == bisp1 || i == bisp2)) {
                        Rcpp::Rcout << " " << "---";
                    } else {
                        Rcpp::Rcout << " " << arma::mean(accept_beta0s.col(i))  / ( (double) s + 1.0);
                    }
                }
                Rcpp::Rcout << "\n";
                
                Rcpp::Rcout << "beta1 accept rates       :";
                for(int i=0; i!=n_species+1; ++i)
                    Rcpp::Rcout << " " << ( (double) accept_beta1(i))  / ( (double) s + 1.0);
                Rcpp::Rcout << "\n";

                Rcpp::Rcout << "  eta^2 = " << eta2 << "\n"
                            << "  sigma^2 = " << sigma2.t();

                
                if (use_bivariate)
                {
                    Rcpp::Rcout << "w0 and w1 (avg) accept rates:";
                    Rcpp::Rcout << " " << arma::mean(accept_beta0s.col(bisp1))  / ( (double) s + 1.0);
                    Rcpp::Rcout << " " << arma::mean(accept_beta0s.col(bisp2))  / ( (double) s + 1.0);
                    Rcpp::Rcout << "\n";

                    Rcpp::Rcout << "A accept rate:"
                                << " " << ((double) accept_A(0)) / ( (double) s + 1.0) 
                                << " " << ((double) accept_A(1)) / ( (double) s + 1.0) 
                                << " " << ((double) accept_A(2)) / ( (double) s + 1.0) 
                                << "\n";
                           
                    Rcpp::Rcout << " A11 = " << A11 << "\n"
                                << " A21 = " << A21 << "\n"
                                << " A22 = " << A22 << "\n";
                }
                
                if (!fix_phi)
                {
                    Rcpp::Rcout << "phi/xi accept rates: "
                                << arma::conv_to< arma::vec >::from(accept_phi).t()  / ( (double) s + 1.0);
                    
                    Rcpp::Rcout << " eff. range =";
                    for(int i=0; i!=n_species; i++)
                        Rcpp::Rcout << " " << floor(3/phi(i));

                    Rcpp::Rcout << " " << floor(3/xi) << "\n";
                }


                Rcpp::Rcout << "\n";
            }
        }

        if (iter >= nburn || nburn == niter)
        {
            if (predict)
            {
                ////
                // Update cross covariance matrices if needed
                ////

                if (!fix_phi)
                {
                    for(int i = 0; i != n_species; ++i)
                    {
                        Sigma_P[i] = arma::exp(-phi(i) * d_P);
                        Sigma_CIF_P[i] = arma::exp(-phi(i) * d_CIF_P);
                    }

                    Sigma_o_P = arma::exp(-xi * d_P);
                    Sigma_o_CIF_P = arma::exp(-xi * d_CIF_P);
                }


                ////
                // Update Predict w0 and w1
                ////

                if (use_bivariate)
                {

                    {
                        arma::mat tmp = Sigma_CIF_P[bisp1].t() * Sigma_CIF_inv[bisp1];

                        arma::mat S = Sigma_P[bisp1] - (tmp * Sigma_CIF_P[bisp1]);
                        arma::vec mu = tmp * beta0s.col(bisp1);

                        arma::mat S_U;
                        if (!arma::chol(S_U, S))
                            throw std::runtime_error("Prediction w0_P Cholesky failed!");

                        w0_P = mu + S_U.t() * arma::randn<arma::vec>(nr_P);
                    }

                    {
                        arma::mat tmp = Sigma_CIF_P[bisp2].t() * Sigma_CIF_inv[bisp2];

                        arma::mat S = Sigma_P[bisp2] - (tmp * Sigma_CIF_P[bisp2]);
                        arma::vec mu = tmp * beta0s.col(bisp2);

                        arma::mat S_U;
                        if (!arma::chol(S_U, S))
                            throw std::runtime_error("Prediction w1_P Cholesky failed!");

                        w1_P = mu + S_U.t() * arma::randn<arma::vec>(nr_P);
                    }


                    beta0s_P.col(bisp1) = A11 * w0_P;
                    beta0s_P.col(bisp2) = A21 * w0_P + A22 * w1_P;
                }


                ////
                // Update Predict \beta_0(s)
                ////

                // 5 Species
                for(int i = 0; i != n_species; ++i)
                {
                    if (use_bivariate && (i == bisp1 || i == bisp2))
                    {
                        continue;
                    }

                    arma::mat S = sigma2(i)*Sigma_P[i] - sigma2(i)*(Sigma_CIF_P[i].t() * Sigma_CIF_inv[i] * Sigma_CIF_P[i]);
                    arma::vec mu = Sigma_CIF_P[i].t() * Sigma_CIF_inv[i] * beta0s.col(i);

                    arma::mat S_U;
                    if (!arma::chol(S_U, S))
                        throw std::runtime_error("Prediction beta_0(s) Cholesky failed!");

                    beta0s_P.col(i) = mu + S_U.t() * arma::randn<arma::vec>(nr_P);
                }

                // Other
                {
                    arma::mat S = eta2*Sigma_o_P - eta2*(Sigma_o_CIF_P.t() * Sigma_o_CIF_inv * Sigma_o_CIF_P);
                    arma::vec mu = Sigma_o_CIF_P.t() * Sigma_o_CIF_inv * beta0s.col(other_index);

                    arma::mat S_U;
                    if (!arma::chol(S_U, S))
                        throw std::runtime_error("Prediction beta_0(s) Cholesky failed!");

                    beta0s_P.col(other_index) = mu + S_U.t() * arma::randn<arma::vec>(nr_P);
                }


                ////
                // Update Predict Z_i
                ////
                for(int i = 0; i != n_species; ++i)
                {
                    Z_i_star_P.col(i) = beta0(i) + beta0s_P.col(i) + beta1(i) * Q_P.col(i);
                }
                Z_i_P = zero_trunc(Z_i_star_P);

                Z_o_P = zero_trunc(beta0(other_index) + beta0s_P.col(other_index) + beta1(other_index) * Q_P.col(other_index));


                ////
                // Update Predict Z_tot
                ////
                {
                    Z_tot_P = Z_o_P;
                    for(int i = 0; i != n_species; ++i)
                        Z_tot_P += Z_i_P.col(i);
                }
            }

            // Save posterior samples

            int i = (niter==nburn) ? iter : iter - nburn;

            post_beta0.row(i) = beta0.t();
            
            post_beta0s_1.row(i) = beta0s.col(0).t();
            post_beta0s_2.row(i) = beta0s.col(1).t();
            post_beta0s_3.row(i) = beta0s.col(2).t();
            post_beta0s_4.row(i) = beta0s.col(3).t();
            post_beta0s_5.row(i) = beta0s.col(4).t();
            post_beta0s_6.row(i) = beta0s.col(5).t();

            post_beta1.row(i) = beta1.t();

            post_Z_o.row(i) = Z_o.t();
            post_Z_1.row(i) = Z_i.col(0).t();
            post_Z_2.row(i) = Z_i.col(1).t();
            post_Z_3.row(i) = Z_i.col(2).t();
            post_Z_4.row(i) = Z_i.col(3).t();
            post_Z_5.row(i) = Z_i.col(4).t();

            post_tau2_1.row(i) = tau2_i.col(0).t();
            post_tau2_2.row(i) = tau2_i.col(1).t();
            post_tau2_3.row(i) = tau2_i.col(2).t();
            post_tau2_4.row(i) = tau2_i.col(3).t();
            post_tau2_5.row(i) = tau2_i.col(4).t();

            post_tau2_tot.row(i) = tau2_tot.t();

            post_sigma2.row(i) = sigma2.t();
            post_eta2(i) = eta2;

            post_beta0s_1_P.row(i) = beta0s_P.col(0).t();
            post_beta0s_2_P.row(i) = beta0s_P.col(1).t();
            post_beta0s_3_P.row(i) = beta0s_P.col(2).t();
            post_beta0s_4_P.row(i) = beta0s_P.col(3).t();
            post_beta0s_5_P.row(i) = beta0s_P.col(4).t();
            post_beta0s_6_P.row(i) = beta0s_P.col(5).t();

            if (use_bivariate)
            {
                post_w0.row(i) = w0.t();
                post_w1.row(i) = w1.t();
                
                post_A(i,0) = A11;
                post_A(i,1) = A21;
                post_A(i,2) = A22;

                post_w0_P.row(i) = w0_P.t();
                post_w1_P.row(i) = w1_P.t();
            }

            if (!fix_phi)
            {
                post_phi.row(i) = phi.t();
                post_xi(i,0) = xi;
            }


            post_Z_o_P.row(i) = Z_o_P.t();
            post_Z_1_P.row(i) = Z_i_P.col(0).t();
            post_Z_2_P.row(i) = Z_i_P.col(1).t();
            post_Z_3_P.row(i) = Z_i_P.col(2).t();
            post_Z_4_P.row(i) = Z_i_P.col(3).t();
            post_Z_5_P.row(i) = Z_i_P.col(4).t();

            post_Z_tot_P.row(i) = Z_tot_P.t();
        }
    }

    Rcpp::List settings = Rcpp::List::create(
                            Rcpp::Named("C") = C,
                            Rcpp::Named("I") = I,
                            Rcpp::Named("F") = F,
                            Rcpp::Named("Q_C") = Q_C,
                            Rcpp::Named("Q_I") = Q_I,
                            Rcpp::Named("Q_F") = Q_F,
                            Rcpp::Named("Q_P") = Q_P,
                            Rcpp::Named("s_C") = s_C,
                            Rcpp::Named("s_I") = s_I,
                            Rcpp::Named("s_F") = s_F,
                            Rcpp::Named("s_P") = s_P,
                            Rcpp::Named("phi") = phi,
                            Rcpp::Named("xi")  = xi,
                            Rcpp::Named("niter") = niter,
                            Rcpp::Named("nthin") = nthin,
                            Rcpp::Named("nburn") = nburn,
                            Rcpp::Named("predict") = predict,
                            Rcpp::Named("verbose") = verbose
                          );

    Rcpp::List pred = Rcpp::List::create(
                        Rcpp::Named("beta0s_1_P") = post_beta0s_1_P,
                        Rcpp::Named("beta0s_2_P") = post_beta0s_2_P,
                        Rcpp::Named("beta0s_3_P") = post_beta0s_3_P,
                        Rcpp::Named("beta0s_4_P") = post_beta0s_4_P,
                        Rcpp::Named("beta0s_5_P") = post_beta0s_5_P,
                        Rcpp::Named("beta0s_6_P") = post_beta0s_6_P,
                        Rcpp::Named("Z_o_P") = post_Z_o_P,
                        Rcpp::Named("Z_1_P") = post_Z_1_P,
                        Rcpp::Named("Z_2_P") = post_Z_2_P,
                        Rcpp::Named("Z_3_P") = post_Z_3_P,
                        Rcpp::Named("Z_4_P") = post_Z_4_P,
                        Rcpp::Named("Z_5_P") = post_Z_5_P,
                        Rcpp::Named("Z_tot_P") = post_Z_tot_P
                      );
    
    

    Rcpp::List beta = Rcpp::List::create(
                        Rcpp::Named("beta0") = post_beta0,
                        Rcpp::Named("beta0s_1") = post_beta0s_1,
                        Rcpp::Named("beta0s_2") = post_beta0s_2,
                        Rcpp::Named("beta0s_3") = post_beta0s_3,
                        Rcpp::Named("beta0s_4") = post_beta0s_4,
                        Rcpp::Named("beta0s_5") = post_beta0s_5,
                        Rcpp::Named("beta0s_6") = post_beta0s_6,
                        Rcpp::Named("beta1") = post_beta1
                      );

    Rcpp::List Z = Rcpp::List::create(
                        Rcpp::Named("Z_o") = post_Z_o,
                        Rcpp::Named("Z_1") = post_Z_1,
                        Rcpp::Named("Z_2") = post_Z_2,
                        Rcpp::Named("Z_3") = post_Z_3,
                        Rcpp::Named("Z_4") = post_Z_4,
                        Rcpp::Named("Z_5") = post_Z_5
                   );

    Rcpp::List tau2 = Rcpp::List::create(
                        Rcpp::Named("tau2_1") = post_tau2_1,
                        Rcpp::Named("tau2_2") = post_tau2_2,
                        Rcpp::Named("tau2_3") = post_tau2_3,
                        Rcpp::Named("tau2_4") = post_tau2_4,
                        Rcpp::Named("tau2_5") = post_tau2_5,
                        Rcpp::Named("tau2_tot") = post_tau2_tot
                      );

    Rcpp::List post = Rcpp::List::create(
                        Rcpp::Named("beta") = beta,
                        Rcpp::Named("Z") = Z,
                        Rcpp::Named("tau2") = tau2,
                        Rcpp::Named("sigma2") = post_sigma2,
                        Rcpp::Named("eta2") = post_eta2
                      );

    Rcpp::List accept = Rcpp::List::create(
                            Rcpp::Named("accept_beta0") = accept_beta0,
                            Rcpp::Named("accept_beta0s") = accept_beta0s,
                            Rcpp::Named("accept_beta1") = accept_beta1
                        );

    if (use_bivariate)
    {
        post["A"] = post_A;
        post["w0"] = post_w0;
        post["w1"] = post_w1;

        pred["w0"] = post_w0_P;
        pred["w1"] = post_w1_P;  

        settings["bisp1"] = bisp1;             
        settings["bisp2"] = bisp2;

        accept["accept_A"] = accept_A;
    }

    if (!fix_phi)
    {
        post["phi"] = post_phi;  
        post["xi"] = post_xi;  
        accept["accept_phi"] = accept_phi;
    }

    return Rcpp::List::create(
                Rcpp::Named("settings") = settings,
                Rcpp::Named("post")     = post,
                Rcpp::Named("pred")     = pred,
                Rcpp::Named("accept") = accept
           );
}

