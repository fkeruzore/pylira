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

    # Create a list with the four threshold arguments IF NEEDED.
    # I'll then pass it as an arg to the `lira` func.
    # This is needed because the 4 threshold args (x_threshold, 
    # delta.x_threshold, and same with y) will not function if they
    # are not called properly, and I don't want to write all possible
    # permutations.
    th_args = list()
    cmd_th = ""
    if ("x_threshold" %in% colnames(d)) {
        th_args[["x.threshold"]] <- d$x_threshold
        cmd_th <- paste(cmd_th, "x.threshold,")
    } else th_args[["x.threshold"]] <- rep(NA, length(d$x_obs))
    if ("y_threshold" %in% colnames(d)) {
        th_args[["y.threshold"]] <- d$y_threshold
        cmd_th <- paste(cmd_th, "y.threshold,")
    } else th_args[["y.threshold"]] <- rep(NA, length(d$x_obs))
    if ("delta_x_threshold" %in% colnames(d)) {
        th_args[["delta.x.threshold"]] <- d$delta_x_threshold
        cmd_th <- paste(cmd_th, "delta.x.threshold,")
    } else th_args[["delta.x.threshold"]] <- rep(NA, length(d$x_obs))
    if ("delta_y_threshold" %in% colnames(d)) {
        th_args[["delta.y.threshold"]] <- d$delta_y_threshold
        cmd_th <- paste(cmd_th, "delta.y.threshold,")
    } else th_args[["delta.y.threshold"]] <- rep(NA, length(d$x_obs))
    if (cmd_th != "") {
        cmd_th <- paste("Threshold keywords recognized:", cmd_th)
        print(cmd_th)
    }

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
    return(samples[[1]][[1]])
}
