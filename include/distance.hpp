#ifndef DISTANCE_HPP
#define DISTANCE_HPP


template <typename T> inline int sign(T val) 
{
    return (T(0) < val) - (val < T(0));
}


arma::mat circle_dist(arma::mat const& x1, bool miles = true) 
{

    double R = 6378.388;
    if (miles) 
        R = 3963.34;

    arma::vec coslat1 = arma::cos( (x1.col(1) * arma::datum::pi)/180 );
    arma::vec sinlat1 = arma::sin( (x1.col(1) * arma::datum::pi)/180 );
    arma::vec coslon1 = arma::cos( (x1.col(0) * arma::datum::pi)/180 );
    arma::vec sinlon1 = arma::sin( (x1.col(0) * arma::datum::pi)/180 );
    
    arma::mat m1(x1.n_rows, 3);

    m1.col(0) = coslat1 % coslon1;
    m1.col(1) = coslat1 % sinlon1;
    m1.col(2) = sinlat1;

    arma::mat pp = arma::acos(m1 * m1.t());

    for(int i = 0; i != pp.n_rows; ++i)
    {
        for(int j = 0; j != pp.n_cols; ++j)
        {
            if (abs(pp(i,j)) > 1)
                pp(i,j) = 1.0 * sign(pp(i,j));
        }   

        pp(i,i) = 0.0;
    }

    return R * pp;
}

arma::mat circle_dist(arma::mat const& x1, arma::mat const& x2, bool miles = true) 
{
    double R = 6378.388;
    if (miles) 
        R = 3963.34;

    arma::vec coslat1 = arma::cos( (x1.col(1) * arma::datum::pi)/180 );
    arma::vec sinlat1 = arma::sin( (x1.col(1) * arma::datum::pi)/180 );
    arma::vec coslon1 = arma::cos( (x1.col(0) * arma::datum::pi)/180 );
    arma::vec sinlon1 = arma::sin( (x1.col(0) * arma::datum::pi)/180 );
    
    arma::vec coslat2 = arma::cos( (x2.col(1) * arma::datum::pi)/180 );
    arma::vec sinlat2 = arma::sin( (x2.col(1) * arma::datum::pi)/180 );
    arma::vec coslon2 = arma::cos( (x2.col(0) * arma::datum::pi)/180 );
    arma::vec sinlon2 = arma::sin( (x2.col(0) * arma::datum::pi)/180 );

    arma::mat m1(x1.n_rows, 3);

    m1.col(0) = coslat1 % coslon1;
    m1.col(1) = coslat1 % sinlon1;
    m1.col(2) = sinlat1;

    arma::mat m2(x2.n_rows, 3);

    m2.col(0) = coslat2 % coslon2;
    m2.col(1) = coslat2 % sinlon2;
    m2.col(2) = sinlat2;

    arma::mat pp = arma::acos(m1 * m2.t());

    for(int i = 0; i != pp.n_rows; ++i)
    {
        for(int j = 0; j != pp.n_cols; ++j)
        {
            if (abs(pp(i,j)) > 1)
                pp(i,j) = 1.0 * sign(pp(i,j));
        }   
    }

    return R * pp;
}

#endif