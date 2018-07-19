library(Rcpp)
library(fields)

save = TRUE

verbose = TRUE
predict = TRUE
bivar   = FALSE
fix_phi = FALSE
fix_beta1 = TRUE


niter = as.integer(5000)
nthin = as.integer(10)
nburn = as.integer(4000)

weeks = 1:52

data_dir = path.expand(file.path(getwd(),"Data"))

if (fix_beta1) {
    res_dir = "output/krige"
} else {
    res_dir = "output/joint"
}

if (bivar)
{
    bisp1 = 1 #sulfate
    bisp2 = 3 #ammonium

    res_dir = paste0(res_dir,"_biv")
}

if (fix_phi) {
    res_dir = paste0(res_dir,"_phi")
} else {
    res_dir = paste0(res_dir,"_vphi")
}

dir.create(res_dir, showWarnings=FALSE, recursive=TRUE)


registerPlugin("local_include", 
                function() {
                    list(env = list(PKG_CXXFLAGS=paste0("-I",path.expand(file.path(getwd(),"include")))))
                }
              )
sourceCpp(paste0("downscaler_joint_tobit.cpp"),rebuild=TRUE, verbose=FALSE)


load(file=file.path(data_dir,paste0("settings.Rdata")))


for(i in weeks)
{    
    l[[i]]$F = as.matrix(l[[i]]$F, ncol=1)

    cat(i, "/", max(weeks), "\n")
    
    C   = l[[i]]$C
    I   = l[[i]]$I
    F   = l[[i]]$F
    P   = l[[i]]$P
    Q_C = l[[i]]$Q_C
    Q_I = l[[i]]$Q_I
    Q_F = l[[i]]$Q_F
    Q_P = l[[i]]$Q_P
    s_C = l[[i]]$s_C
    s_I = l[[i]]$s_I
    s_F = l[[i]]$s_F
    s_P = l[[i]]$s_P

    # FIXME
    dist = rdist.earth(s_F, miles = TRUE)    


    phi = 3 / c(700,650,725,700,580)
    xi = 3/600

    beta0 = rep(0,6)
    beta1 = rep(0,6)


    # Calc CMAQ other
    Q_C[,6] = Q_C[,6] - apply(Q_C[,1:5],1,sum)
    Q_I[,6] = Q_I[,6] - apply(Q_I[,1:5],1,sum)
    Q_F[,6] = Q_F[,6] - apply(Q_F[,1:5],1,sum)
    Q_P[,6] = Q_P[,6] - apply(Q_P[,1:5],1,sum)

    stopifnot(all(Q_C[,6]>=0))
    stopifnot(all(Q_I[,6]>=0))
    stopifnot(all(Q_F[,6]>=0))
    stopifnot(all(Q_P[,6]>=0))

    if (bivar)
    {
        r = downscaler_joint( C, I, F,
                              Q_C, Q_I, Q_F, Q_P,
                              s_C, s_I, s_F, s_P,
                              phi, xi,
                              beta0, beta1,
                              niter,  nthin, nburn,
                              predict = predict,
                              verbose = verbose,
                              profile = profile,
                              fix_phi = fix_phi,
                              fix_beta1 = fix_beta1,
                              bisp1 = bisp1, bisp2 = bisp2
                            )
    } else {
        r = downscaler_joint( C, I, F,
                              Q_C, Q_I, Q_F, Q_P,
                              s_C, s_I, s_F, s_P,
                              phi, xi,
                              beta0, beta1,
                              niter,  nthin, nburn,
                              predict = predict,
                              verbose = verbose,
                              profile = profile,
                              fix_phi = fix_phi,
                              fix_beta1 = fix_beta1
                            )  
    }
    r$settings$P = P

    if (save)
        save(r, file=file.path(res_dir,paste0(i,".Rdata")))
}

