library(lubridate)
library(plyr)


CSN_INPUT  = "Data/csn_data_2007.csv"
CSN_OUTPUT = "Data/csn.Rdata"

IMPROVE_INPUT  = "Data/improve_data_2007.csv"
IMPROVE_LOCS   = "Data/improve_sites.csv"
IMPROVE_OUTPUT = "Data/improve.Rdata"

FRM_INPUT  = "Data/AQS_PM25_2007_DS.csv"
FRM_OUTPUT = "Data/frm.Rdata"

CMAQ_INPUT  = "Data/CMAQ_daily_PM_Species_2007.csv.gz"
CMAQ_OUTPUT = "Data/cmaq.Rdata"

source("util/distance_util.R")

na_mean = function(x, na.rm)
{
    m = mean(x, na.rm=na.rm)
    if (is.finite(m)) return(m)
    else return(NA)
}

default_env = c(ls(),"default_env")

# Process CSN data

csn = read.csv(CSN_INPUT,as.is=TRUE,na.strings=".")

site_cols = c("SITE", "STATE_NAME", "COUNTY_NAME", "CITY", 
              "CITY_NAME", "ADDRESS", "AQCR", "AQCR_NAME", "UAR", "UAR_NAME", 
              "LAND_USE", "LOCATION", "LATITUDE", "LONGITUDE")
csn_sites = unique(csn[,site_cols])

csn_full = csn

cols = c("SITE", "LONGITUDE", "LATITUDE", "DATE", "POC", 
         "m_so4", "m_no3", "m_nh4", "oc_adj", "ec_niosh", "PM2.5.Mass")

csn = csn[,cols]

colnames(csn) = c("site", "longitude", "latitude", "date", "POC", 
                  "sulfate", "nitrate", "ammonium", "oc", "ec", "pm25")

csn$carbon = 1.4*csn$oc + csn$ec
csn$soil = with(csn_full, 2.2*aluminum + 2.49*silicon + 1.63*calcium + 2.42*iron + 1.94*titanium)

csn$date = mdy(csn$date)


csn = ddply(csn,
            .(site, longitude, latitude, date),
            summarize,
            sulfate  = na_mean(sulfate , na.rm=TRUE),
            nitrate  = na_mean(nitrate , na.rm=TRUE),
            ammonium = na_mean(ammonium, na.rm=TRUE),
            soil     = na_mean(soil    , na.rm=TRUE),
            oc       = na_mean(oc      , na.rm=TRUE),
            ec       = na_mean(ec      , na.rm=TRUE),
            carbon   = na_mean(carbon  , na.rm=TRUE),
            pm25     = na_mean(pm25    , na.rm=TRUE),  
            .progress="text")


# Check for repeated observations
cols = c("site", "longitude", "latitude", "date")
stopifnot(all(dim(csn[,cols]) == dim(unique(csn[,cols]))))


save(csn, csn_full, file = CSN_OUTPUT)





# Process Improve data

improve = read.csv(IMPROVE_INPUT,as.is=TRUE, na.strings="-999")
imp_site_locs = read.csv(IMPROVE_LOCS,as.is=TRUE,header=FALSE)
colnames(imp_site_locs) = c("site_code","longitude","latitude","unknown")

improve = merge(imp_site_locs[,-4],improve,all.y=TRUE)
rm(imp_site_locs)

improve_full = improve

cols = c("site_code", "longitude", "latitude", "obs_date", "POC", 
         "SO4f_val", "NO3f_val", "NH4f_val", "OCf_val", "ECf_val", "MF_val")


improve = improve[,cols]
colnames(improve) = c("site", "longitude", "latitude", "date", "POC", 
                      "sulfate", "nitrate", "ammonium", "oc", "ec", "pm25")



improve$carbon = 1.4*improve$oc + improve$ec
improve$soil = with(improve_full, 2.20*ALf_val + 2.49*SIf_val + 1.63*CAf_val + 2.42*FEf_val + 1.94*TIf_val)

improve$ammonium = with(improve_full, 18/62*NO3f_val + 36/96*SO4f_val)

improve$date = ymd(improve$date)

improve = ddply(improve,
                .(site, longitude, latitude, date),
                summarize,
                sulfate  = na_mean(sulfate , na.rm=TRUE),
                nitrate  = na_mean(nitrate , na.rm=TRUE),
                ammonium = na_mean(ammonium, na.rm=TRUE),
                soil     = na_mean(soil    , na.rm=TRUE),
                oc       = na_mean(oc      , na.rm=TRUE),
                ec       = na_mean(ec      , na.rm=TRUE),
                carbon   = na_mean(carbon  , na.rm=TRUE),
                pm25     = na_mean(pm25    , na.rm=TRUE),  
                .progress="text")

# Check for repeated observations
cols = c("site", "longitude", "latitude", "date")
stopifnot(all(dim(improve[,cols]) == dim(unique(improve[,cols]))))


save(improve, improve_full, file = IMPROVE_OUTPUT)


# Process FRM data

frm = read.csv(FRM_INPUT, colClasses=c("character","integer","character",rep("numeric",3)))

names(frm) = tolower(names(frm))
frm$date = mdy(frm$date)

frm[frm$pm25 == -999, ] = NA

frm_full = frm

frm = ddply(frm,
            .(site, date, latitude, longitude),
            summarize,
            pm25  = na_mean(pm25 , na.rm=TRUE),       
            .progress="text")

stopifnot(!any(is.nan(as.matrix(frm[,5:ncol(frm)]) )))
stopifnot(!any(duplicated(frm[,1:2])))

save(frm, frm_full, file=FRM_OUTPUT)


# Process CMAQ data

nc = ncol(read.csv(textConnection(readLines(CMAQ_INPUT,5)[5]),head=FALSE))

cmaq = read.csv(CMAQ_INPUT, skip=4, header=FALSE, colClasses=c(rep("integer",2), rep("numeric",2),"character",rep("numeric",nc-5)))
gc()

names(cmaq) = c("column","row","longitude","latitude","date","sulfate","nitrate","ammonium","oc","ec","soil","pm25")

cmaq$column = NULL
cmaq$row = NULL
gc()


cmaq$carbon = 1.4*cmaq$oc + cmaq$ec
cmaq$date = ymd(cmaq$date)


cmaq_locs = unique(cmaq[,c("longitude", "latitude")])
colnames(cmaq_locs) = c("longitude","latitude")

save(cmaq, cmaq_locs, file = CMAQ_OUTPUT)


rm( list = setdiff(ls(),default_env) )



# Weekly Aggregate


data_files = list("csn"     = CSN_OUTPUT,    
                  "improve" = IMPROVE_OUTPUT,
                  "frm"     = FRM_OUTPUT,    
                  "cmaq"    = CMAQ_OUTPUT   
                  )



load(file = data_files[["cmaq"]])
names(cmaq)[1:2] = c("cmaq_long","cmaq_lat")

cmaq$week = week(cmaq$date)


  


for(suff in c("frm","csn","improve")) 
{
    load(file = data_files[[suff]])

    data = get(suff)

    rm(list=c(suff,paste0(suff,"_full")))



    data$week = week(data$date)

    if (suff == "frm")
    {
        data = ddply(data,
                     .(week, site, longitude, latitude),
                     summarize,
                     pm25     = na_mean(pm25    , na.rm=TRUE),    
                     .progress="text"
                    )
    }
    else
    {
        data = ddply(data,
                     .(week, site, longitude, latitude),
                     summarize,
                     sulfate  = na_mean(sulfate , na.rm=TRUE),    
                     nitrate  = na_mean(nitrate , na.rm=TRUE),    
                     ammonium = na_mean(ammonium, na.rm=TRUE),    
                     soil     = na_mean(soil    , na.rm=TRUE),
                     oc       = na_mean(oc      , na.rm=TRUE),
                     ec       = na_mean(ec      , na.rm=TRUE),
                     carbon   = na_mean(carbon  , na.rm=TRUE),    
                     pm25     = na_mean(pm25    , na.rm=TRUE),    
                     .progress="text"
                    )
    }

    data$date = date
    week(data$date) = data$week



    locs = unique(data[,c("site", "longitude", "latitude")])
    rownames(locs) = NULL
    locs$site = as.character(locs$site)
    locs = locs[order(locs$site),]


    d = sp_dist(locs[,c("longitude","latitude")], cmaq_locs[,c("longitude","latitude")])
    z = apply(d,1,which.min)

    locs$cmaq_long = cmaq_locs[z,1]
    locs$cmaq_lat  = cmaq_locs[z,2]
    locs$cmaq_dist = apply(d,1,min)

    # exclude sites not within CMAQ
    locs = locs[locs$cmaq_dist < 1,]

    rm(d); gc()


    d = cmaq[  cmaq$cmaq_long %in% locs$cmaq_long 
             & cmaq$cmaq_lat %in% locs$cmaq_lat
             #& cmaq$date %in% data$date
                    ,  ]

    d = ddply(d,
              .(week, cmaq_long, cmaq_lat),
              summarize,
              sulfate  = na_mean(sulfate , na.rm=TRUE),    
              nitrate  = na_mean(nitrate , na.rm=TRUE),    
              ammonium = na_mean(ammonium, na.rm=TRUE),    
              soil     = na_mean(soil    , na.rm=TRUE),
              oc       = na_mean(oc      , na.rm=TRUE),
              ec       = na_mean(ec      , na.rm=TRUE),
              carbon   = na_mean(carbon  , na.rm=TRUE),    
              pm25     = na_mean(pm25    , na.rm=TRUE),    
              .progress="text"
             )

    d$date = date
    week(d$date) = d$week

    d = d[d$date %in% data$date,] # not all cmaq dates will be in each data source

    d = merge(locs[,c("site", "cmaq_long", "cmaq_lat")],d)

    
    fix_shape = function(data, sps)
    {

        data_sp = list()
        
        for(sp in sps)
        {
            cat("  ",sp,"...\n")

            tmp = reshape(data[,c("site","date",sp)], v.names = sp, idvar = "date", timevar = "site", direction = "wide")
            
            rownames(tmp) = tmp$date
            tmp$date = NULL
            colnames(tmp) = substr(colnames(tmp),nchar(sp)+2,nchar(colnames(tmp)))

            tmp = tmp[order(rownames(tmp)),]
            tmp = tmp[,order(colnames(tmp))]

            data_sp[[sp]] = tmp
        }

        return(data_sp)   
    }

    sps = intersect( colnames(data), c("sulfate", "nitrate", "ammonium", "soil", "oc", "ec", "carbon", "pm25"))
    sps = c("sulfate", "nitrate", "ammonium", "soil", "oc", "ec", "carbon", "pm25")

    cat("CMAQ\n")
    cmaq_sp = fix_shape(d, sps)

    data = data[data$site %in% locs$site,]

    dates = sort(unique(d$date))
    sites = sort(unique(data$site))
    
    all_dates = data.frame(site = rep(sites, each=length(dates)),
                           date = rep(dates, times=length(sites)) 
                          )

    data = merge(all_dates, data, all=TRUE)

    cat(suff,"\n")
    sps = intersect( colnames(data), c("sulfate", "nitrate", "ammonium", "soil", "oc", "ec", "carbon", "pm25"))
    data_sp = fix_shape(data, sps)

    stopifnot(all(sapply(data_sp, function(x) all(colnames(x) == locs$site))))
    stopifnot(all(sapply(cmaq_sp, function(x) all(colnames(x) == locs$site))))

    save(locs, cmaq_sp, data_sp, file=file.path("Data",paste0("combined_",suff,".Rdata")))
}


# Format data for downscaler

data = list()

for(suff in c("csn","improve","frm"))
{
    load(file=file.path("Data",paste0("combined_",suff,".Rdata")))

    if (suff == "frm")
    {
        data_sp = lapply(data_sp, function(x) x[-53,])
        cmaq_sp = lapply(cmaq_sp, function(x) x[-53,])
    }

    data[[suff]] = list()
    data[[suff]]$species = names(data_sp)
    data[[suff]]$data = lapply(data_sp, as.matrix)
    data[[suff]]$cmaq = cmaq_sp
    data[[suff]]$locs = locs[,2:3]

    data[[suff]]$n_dates = unique(c( sapply(data_sp, nrow), sapply(cmaq_sp, nrow)))
    data[[suff]]$n_sites = unique(c( sapply(data_sp, ncol), sapply(cmaq_sp, ncol)))

    rm(data_sp, cmaq_sp, locs)
}


n_dates = unique(sapply(data, function(x) x$n_dates))
sp = c("sulfate", "nitrate", "ammonium", "soil", "carbon", "pm25")

stopifnot(length(n_dates) == 1)

save(data, n_dates, sp, file="Data/data.Rdata")


l = list()
for(i in 1:n_dates)
{
    cat(i,"/",n_dates,"\n")

    C = sapply(sp,     function(s) data$csn$data[[s]][i,])
    I = sapply(sp,     function(s) data$improve$data[[s]][i,])
    F = matrix(sapply("pm25", function(s) data$frm$data[[s]][i,]), ncol=1)

    C[C<=0] = 0
    I[I<=0] = 0
    F[F<=0] = 0

    Q_C = sapply(sp, function(s) as.numeric( data$csn$cmaq[[s]][i,]) )
    Q_I = sapply(sp, function(s) as.numeric( data$improve$cmaq[[s]][i,]) )
    Q_F = sapply(sp, function(s) as.numeric( data$frm$cmaq[[s]][i,]) )

    s_C = as.matrix(data$csn$locs)
    s_I = as.matrix(data$improve$locs)
    s_F = as.matrix(data$frm$locs)

    n_C = nrow(C)
    n_I = nrow(I)
    n_F = nrow(F)

    dups = which(duplicated(rbind(s_C,s_I,s_F)))

    stopifnot(all(dups > n_C+n_I))

    dup_idx = dups - (n_C+n_I)

    F = F[-dup_idx,,drop=FALSE]
    Q_F = Q_F[-dup_idx,,drop=FALSE]
    s_F = s_F[-dup_idx,,drop=FALSE]

    stopifnot(!sum(duplicated(rbind(s_C,s_I,s_F))))

    C_good = apply(C,1, function(x) all(!is.na(x)))
    I_good = apply(I,1, function(x) all(!is.na(x)))
    F_good = apply(F,1, function(x) all(!is.na(x)))

    # Randomly exclude 10% for validation
    C_good[ sample(1:nrow(C), floor(nrow(C)*0.1), replace=FALSE) ] = FALSE
    I_good[ sample(1:nrow(I), floor(nrow(I)*0.1), replace=FALSE) ] = FALSE

    P = rbind(C[!C_good,], I[!I_good,])
    s_P = rbind(s_C[!C_good,], s_I[!I_good,])
    Q_P = rbind(Q_C[!C_good,], Q_I[!I_good,])   

    P = P[!duplicated(s_P),]
    Q_P = Q_P[!duplicated(s_P),]
    s_P = s_P[!duplicated(s_P),]


    C = C[C_good,,drop=FALSE]
    Q_C = Q_C[C_good,,drop=FALSE]
    s_C = s_C[C_good,,drop=FALSE]

    I = I[I_good,,drop=FALSE]
    Q_I = Q_I[I_good,,drop=FALSE]
    s_I = s_I[I_good,,drop=FALSE]

    F = F[F_good,,drop=FALSE]
    Q_F = Q_F[F_good,,drop=FALSE]
    s_F = s_F[F_good,,drop=FALSE]

    l[[i]] = list(C=C,I=I,F=F,P=P,
                  Q_C=Q_C,Q_I=Q_I,Q_F=Q_F,Q_P=Q_P,
                  s_C=s_C,s_I=s_I,s_F=s_F,s_P=s_P)

    rm(C,I,F,P,Q_C,Q_I,Q_F,Q_P,s_C,s_I,s_F,s_P)
}

save(l, file="Data/settings.Rdata")



