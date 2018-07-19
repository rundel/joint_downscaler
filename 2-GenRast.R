library(raster)
library(lubridate)

load(file="Data/cmaq.Rdata")

sp = c("sulfate", "nitrate", "ammonium", "soil", "carbon", "pm25")

seasons = c("spring","summer","fall","winter")
month_to_season = c(4,4,1,1,1,2,2,2,3,3,3,4)

date_to_season = function(date)
{
    month_to_season[ month(date) ]
}

cmaq$season = date_to_season(cmaq$date)
cmaq$week = week(cmaq$date)

ratio = 5

r=raster(xmn=-130, xmx=-65, ymn=22, ymx=52, ncols=650/ratio, nrows=300/ratio)

cmaq$cell = cellFromXY(r, cmaq[,1:2])

cmaq = cmaq[,c(sp,"season","week","cell")]
gc()


# Seasonal Rasters

cmaq_rast = list()

for(s in sp) # Species
{
    cat(s)
    cmaq_rast[[s]] = list()
 
    for(i in 1:4) # Season
    {
        cat(i)
        cmaq_rast[[s]][[seasons[i]]] = r

        tmp = cmaq[cmaq$season == i, c(s,"cell")]

        z=tapply(tmp[,1],tmp[,2], mean)

        cmaq_rast[[s]][[seasons[i]]][as.integer(names(z))] = z
    }

    cmaq_rast[[s]][["yearly"]] = r

    z = tapply(cmaq[,s],cmaq[,"cell"], mean)
    cmaq_rast[[s]][["yearly"]][as.integer(names(z))] = z

    cat("\n")
}

save(cmaq_rast,file="Data/cmaq_seasonal_rasts.Rdata")


# Weekly Rasters

cmaq_rast = list()

for(s in sp) # Species
{
    cat(s)
    cmaq_rast[[s]] = list()
 
    for(i in 1:52) # Week
    {
        cat(i)
        cmaq_rast[[s]][[i]] = r

        tmp = cmaq[cmaq$week == i, c(s,"cell")]

        z=tapply(tmp[,1],tmp[,2], mean)

        cmaq_rast[[s]][[i]][as.integer(names(z))] = z
    }

    cat("\n")
}

save(cmaq_rast,file="Data/cmaq_weekly_rasts.Rdata")




# load(file="Data/cmaq_weekly_rasts.Rdata")

# season_months = list(spring = 3:5,
#                      summer = 6:8,
#                      fall   = 9:11,
#                      winter = c(12,1:2),
#                      year = 1:12
#                      )
# week_months = month(mdy("1/1/07")+weeks(1:52-1))

# season_weeks = lapply(season_months, function(x) which(week_months %in% x))

# rasts = list()

# for(i in 1:length(season_weeks))
# {
#     weeks = season_weeks[[i]]

#     data = lapply(cmaq_rast, function(sp) 
#             {
#                 sapply(weeks, function(x) sp[[x]][])
#             })

#     res = lapply(data, function(x)
#             {
#                 r = cmaq_rast[[1]][[1]]
#                 r[] = apply(x,1,mean)
#                 return(r)
#             })

#     rasts[[names(season_weeks)[i]]] = res
# }

# save(rasts, file="Data/joint_cmaq_rasts.Rdata")
