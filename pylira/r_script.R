function(d, nsteps, nmix, ...){
    
    install_has_lira <- require("lira")
    if(!install_has_lira) {
        install.packages("lira")
    }
    library(lira)
    
    # This piece of code solely aims at calling `lira` with the right set of arguments.
    # The way we have to do it is quite cumbersome due to how R function calls work, and
    #    to how rigid LIRA is in what it expects as inputs.
    # To summarize, in order to not use an argument, it needs to _not be included at all_
    #    in the call to `lira`, which makes stuff messy.
    # So, the best way I've found is to:
    # 1. Define a list of inputs with only named entries (essentially a dict), to which 
    #    I add entries *one by one* *only if they are being used*;
    # 2. Use `do.call(lira, inputs)`, which is essentially the R equivalent to Python's 
    #    `lira(**inputs)` (except it could manage both named and unnamed entries, but we 
    #    shouldn't make it even messier)

    inputs <- list()

    # Mandatory data inputs
    inputs$x             <- d$x_obs
    inputs$y             <- d$y_obs
    inputs$delta.x       <- d$x_err
    inputs$delta.y       <- d$y_err
    inputs$covariance.xy <- d$corr * d$x_err * d$y_err

    # Threshold arguments: only add those included in the dataframe
    cmd_th = "" # For printing purposes later
    if ("y_threshold" %in% colnames(d)) {
        inputs$y.threshold <- d$y_threshold
        cmd_th <- paste(cmd_th, "y.threshold,")
    }

    if ("delta_y_threshold" %in% colnames(d)) {
        inputs$delta.y.threshold <- d$delta_y_threshold
        cmd_th <- paste(cmd_th, "delta.y.threshold,")
    }

    if ("x_threshold" %in% colnames(d)) {
        inputs$x.threshold <- d$x_threshold
        cmd_th <- paste(cmd_th, "x.threshold,")
    }

    if ("delta_x_threshold" %in% colnames(d)) {
        inputs$delta.x.threshold <- d$delta_x_threshold
        cmd_th <- paste(cmd_th, "delta.x.threshold,")
    }

    # Extra arguments (custom priors, etc) passed in the `lira_kwargs` python dict
    extra_args <- list(...)
    for (arg in names(extra_args)) { inputs[[arg]] <- extra_args[[arg]] }

    # <andatory technical inputs
    inputs$n.mixture <- nmix
    inputs$n.iter <- nsteps
    inputs$n.adapt <- nsteps/4  # This is the burn-in

    # Printing technical stuff (number of steps, custom priors, etc)
    cmd <- paste(
        "Running: ", 
        "nsteps=", nsteps, ", ",
        "nmix=", nmix, ", ",
        sep=""
    )
    for (arg in names(extra_args)){
        cmd <- paste(cmd, arg, "=", extra_args[[arg]], ", ", sep="")
    }
    print(substr(cmd, 1, nchar(cmd)-2), quote=F)

    # Printing summoned threshold arguments
    if (cmd_th != "") {
        cmd_th <- paste("Threshold keywords recognized:", cmd_th)
        print(cmd_th)
    }

    samples <- do.call(lira, inputs)

    return(samples[[1]][[1]])
}
