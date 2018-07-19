library(RcppArmadillo)
library(inline)


euclid = cxxfunction( signature(X_r = "numeric",Y_r = "numeric"),
'
    arma::mat X = Rcpp::as<arma::mat>(X_r);
    int n = X.n_rows;

    arma::mat Y = Rcpp::as<arma::mat>(Y_r);
    int m = Y.n_rows;

    if (X.n_cols!=Y.n_cols) throw("Dimension mismatch between coordinates.");

    arma::mat D(n,m);
    for(int i=0; i!=n; ++i)
    {
        for(int j=0; j!=m; ++j)
        {
           D.at(i,j) = sqrt(arma::accu(arma::square(X.row(i)-Y.row(j))));
        }
    }

    return Rcpp::wrap(D);
', plugin = "RcppArmadillo" )

euclid_sym = cxxfunction( signature(X_r = "numeric"),
'
    arma::mat X = Rcpp::as<arma::mat>(X_r);
    int n = X.n_rows;

    arma::mat D(n,n);
    D.diag() = arma::zeros<arma::vec>(n);

    if (n != 1)
    {
        for(int i=1; i<n; ++i)
        {
            for(int j=0; j<i; ++j)
            {
               D.at(i,j) = sqrt(arma::accu(arma::square(X.row(i)-X.row(j))));
               D.at(j,i) = D.at(i,j);
            }
        }
    }

    return Rcpp::wrap(D);
', plugin = "RcppArmadillo" )


sp_dist = function(x, y=NULL, method="euclidean", ...)
{
    args = list(...)

    if (!is.na(pmatch(method, "euclidian"))) 
        method = "euclidean"
    METHODS = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
    
    method = pmatch(method, METHODS)
    if (is.na(method)) stop("invalid distance method")
    if (method == -1)  stop("ambiguous distance method")
    method = METHODS[method]

    x = as.matrix(x)
    if (!missing(y)) {
        y = as.matrix(y)
        if (ncol(x) != ncol(y)) 
            stop("Dimensionality of x and y mismatch")
    }

    if (method == "euclidean")
    {
        if (missing(y)) return( euclid_sym(x) )
        else            return( euclid(x, y) )
    }
    else
    {
        stop("Unsupported distance function.")
    }
}