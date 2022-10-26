function(d, nsteps, nmix, ...){
    
    install_has_lira <- require("lira")
    if(!install_has_lira) {
        install.packages("lira")
    }
    library(lira)
    
    cmd <- paste(
        "Running: ", 
        "nsteps=", nsteps, ", ",
        "nmix=", nmix, ", ",
        sep=""
    )
    args <- list(...)
    for (arg in names(args)){
        cmd <- paste(cmd, arg, "=", args[[arg]], ", ", sep="")
    }
    print(substr(cmd, 1, nchar(cmd)-2), quote=F)

    if("y_threshold" %in% colnames(d)) {
        print("Threshold detected")
        samples <- lira(
            d$x_obs,
            d$y_obs,
            delta.x=d$x_err,
            delta.y=d$y_err,
            covariance.xy=d$corr * d$x_err * d$y_err,
            y.threshold=d$y_threshold,
            n.mixture=nmix,
            n.iter=nsteps,
            n.adapt=nsteps/4,
            ...
    )
    } else {
        print("No threshold detected")
        samples <- lira(
            d$x_obs,
            d$y_obs,
            delta.x=d$x_err,
            delta.y=d$y_err,
            covariance.xy=d$corr * d$x_err * d$y_err,
            n.mixture=nmix,
            n.iter=nsteps,
            n.adapt=nsteps/4,
            ...
        )
    }
    return(samples[[1]][[1]])
}
