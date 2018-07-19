library(coda)
library(lubridate)
library(Rcpp)
library(raster)


data_dir = path.expand(file.path(getwd(),"Data"))


registerPlugin("local_include", 
                function() {
                    list(env = list(PKG_CXXFLAGS=paste0("-I",path.expand(file.path(getwd(),"include")))))
                }
              )
sourceCpp("downscaler_joint_predict.cpp",rebuild=TRUE, verbose=FALSE)




sp = c("sulfate", "nitrate", "ammonium", "soil", "carbon", "pm25")


season_months = list(spring = 3:5,
                     summer = 6:8,
                     fall   = 9:11,
                     winter = c(12,1:2),
                     year = 1:12
                     )
week_months = month(mdy("1/1/07")+weeks(1:52-1))

season_weeks = lapply(season_months, function(x) which(week_months %in% x))

load(file=file.path(data_dir,"cmaq_weekly_rasts.Rdata"))


###
#
# Parameters
#
###

n_samples = 100
runs = c("joint_vphi")

for(run in runs)
{
    cat(run,"\n",sep="")
    pb = txtProgressBar(style=3)

    preds = list()
    locs = list()

    for(i in 1:52)
    {
        file = file.path(getwd(), "output", run, paste0(i,".Rdata"))
        if(!file.exists(file))
            next

        load(file=file)

        locs[[i]] = rbind(r$settings$s_C,
                          r$settings$s_I,
                          r$settings$s_F)

        thin = nrow(r$post$beta$beta0)/n_samples

        stopifnot(!is.null(r$post$phi))


        phi = window(mcmc(r$post$phi), thin=thin)
        xi  = window(mcmc(r$post$xi), thin=thin)

        eta2   = window(mcmc(r$post$eta2  ), thin=thin)
        sigma2 = window(mcmc(r$post$sigma2), thin=thin)
        
        beta0 = window(mcmc(r$post$beta$beta0), thin=thin)
        beta1 = window(mcmc(r$post$beta$beta1), thin=thin)

        beta0s_1 = window(mcmc(r$post$beta$beta0s_1), thin=thin)   
        beta0s_2 = window(mcmc(r$post$beta$beta0s_2), thin=thin)   
        beta0s_3 = window(mcmc(r$post$beta$beta0s_3), thin=thin)   
        beta0s_4 = window(mcmc(r$post$beta$beta0s_4), thin=thin)   
        beta0s_5 = window(mcmc(r$post$beta$beta0s_5), thin=thin)   
        beta0s_o = window(mcmc(r$post$beta$beta0s_6), thin=thin)   

        locs = rbind(r$settings$s_C,
                     r$settings$s_I,
                     r$settings$s_F)  
                              
        # Make sure there are no duplicate locations
        stopifnot((nrow(locs) == nrow(unique(locs))))        

        trans = "max"

        rast = cmaq_rast[[1]][[1]]
        rast[] = NA

        pred_locs = xyFromCell(rast,1:length(rast[]))

        cmaq = sapply(cmaq_rast, function(x) x[[i]][])

        sub = apply(cmaq,1,function(x) any(is.na(x)))

        cmaq      = cmaq[!sub,]
        pred_locs = pred_locs[!sub,]

        res = downscaler_joint_predict_phi(cmaq, locs, pred_locs,
                                           phi, xi,
                                           eta2, sigma2,
                                           beta0, beta1,
                                           beta0s_1, beta0s_2, beta0s_3,
                                           beta0s_4, beta0s_5, beta0s_o,
                                           trans)   
    

        res$rast = rast
        res$cells = cellFromXY(rast, pred_locs)
        res$pred_locs = pred_locs


        preds[[i]] = res


        rm(r,res)
        gc()
        
        setTxtProgressBar(pb, i/52)
    } 
    close(pb)

    save(preds,
         file = file.path(getwd(), "output", run, paste0(run,"_rast_predict.Rdata")))


    # Aggregate to seasons

    rasts = list()
    stations = list() 

    for(i in 1:length(season_weeks))
    {
        weeks = season_weeks[[i]]

        stations[[names(season_weeks)[i]]] = unique(do.call("rbind",locs[weeks]))

        data = lapply(weeks, function(x) preds[[x]])

        cells = do.call("cbind", lapply(data, function(x) x[["cells"]]))
        stopifnot(all(apply(cells,2,function(x) all.equal(x, cells[,1]))))
        cells = cells[,1]


        tmp = list()
        for(col in c("Z_1", "Z_2", "Z_3", "Z_4", "Z_5", "Z_o", "Z_tot"))
        {
            tmp[[col]] =  do.call("cbind", lapply(data, function(x) x[[col]]))
        }

        res = lapply(tmp, function(x) 
              {
                  r = data[[1]]$rast
                  r[cells] = apply(x,1,mean)

                  return(r)
              })

        names(res) = c("sulfate", "nitrate", "ammonium", "soil", "carbon", "other", "pm25")

        rasts[[names(season_weeks)[i]]] = res
    }

    save(rasts, stations,
         file = file.path(getwd(), "output", run, paste0(run,"_rasts.Rdata")))
}



