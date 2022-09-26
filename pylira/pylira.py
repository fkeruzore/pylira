import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import pandas as pd
from chainconsumer import ChainConsumer
import os
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

_kwargs = {
    "fmt": "o",
    "mec": "k",
    "ecolor": "k",
    "ms": 4,
    "capsize": 2,
    "elinewidth": 1.0,
    "capthick": 1.0,
    "mew": 1.0,
}


class Dataset:
    def __init__(
        self, x_obs, y_obs, x_err, y_err, corr=None, y_threshold=None,
    ):

        """
        Class for your data points.

        Args:
            x_obs (array-like): your data points along the x-axis,
            y_obs (array-like): your data points along the y-axis,
            x_err (array-like): your observationnal errors on points
                along the x-axis,
            y_err (array-like): your observationnal errors on points
                along the y-axis,
            corr (float or array-like): correlation coefficients
                between x and y,
            y_threshold (float or array-like): the cut in y, can be
                an array for box-like selection
        """

        self.x_obs = x_obs
        self.y_obs = y_obs
        self.x_err = x_err
        self.y_err = y_err

        self.n_pts = self.x_obs.shape[0]

        # ======== Correlated errors
        if corr is None:
            corr = np.zeros(self.n_pts)
        elif isinstance(corr, float) or isinstance(corr, int):
            corr = np.ones(self.n_pts) * corr
        self.corr = corr

        # self.covmats[i, :, :] is the error covariance matrix for the ith data point
        self.covmats = np.array(
            [
                [
                    [
                        self.x_err[i] ** 2,
                        corr[i] * self.x_err[i] * self.y_err[i],
                    ],
                    [
                        corr[i] * self.x_err[i] * self.y_err[i],
                        self.y_err[i] ** 2,
                    ],
                ]
                for i in range(self.n_pts)
            ]
        )

        # ======== y threshold for Malmquist
        if isinstance(y_threshold, float) or isinstance(y_threshold, int):
            y_threshold = np.ones(self.n_pts) * y_threshold
        self.y_threshold = y_threshold

        # ======== Misc
        self.xlabel = "x obs"
        self.ylabel = "y obs"
        self.path_to_results = "./"

    # ------------------------------------------------------- #

    @classmethod
    def from_table(
        cls,
        table,
        x_obs="x_obs",
        y_obs="y_obs",
        x_err="x_err",
        y_err="y_err",
        corr=None,
        y_threshold=None,
    ):
        return cls(
            table[x_obs].data,
            table[y_obs].data,
            table[x_err].data,
            table[y_err].data,
            corr=table[corr].data if corr is not None else None,
            y_threshold=table[y_threshold].data
            if y_threshold is not None
            else None,
        )

    # ------------------------------------------------------- #

    def to_table(self):
        d = {}
        d["x_obs"] = self.x_obs
        d["y_obs"] = self.y_obs
        d["x_err"] = self.x_err
        d["y_err"] = self.y_err
        d["corr"] = self.corr
        if self.y_threshold is not None:
            d["y_threshold"] = self.y_threshold
        return Table(d)

    # ------------------------------------------------------- #

    def set_axesnames(self, x, y):
        self.xlabel = x
        self.ylabel = y

    # ------------------------------------------------------- #

    def set_path_to_results(self, path_to_results):
        self.path_to_results = path_to_results
        os.makedirs(path_to_results, exist_ok=True)

    # ------------------------------------------------------- #

    def fit_bces(self, verbose=True):

        # ======== Slope
        beta_num = np.cov(self.x_obs, self.y_obs, ddof=1)[0, 1] - np.average(
            self.covmats[:, 0, 1]
        )
        beta_den = np.var(self.x_obs, ddof=1) - np.average(self.x_err ** 2)
        beta = beta_num / beta_den

        # ======== Intercept
        alpha = np.average(self.y_obs) - beta * np.average(self.x_obs)

        # ======== Scatter, Pratt+09 eqs 3&4
        vi = self.y_err ** 2 + (beta * self.x_err) ** 2
        wi = (1 / vi) / np.average(1 / vi)
        V_raw = (1 / (self.n_pts - 2)) * np.sum(
            wi * (self.y_obs - beta * self.x_obs - alpha) ** 2
        )
        V_orth = V_raw - (
            (1.0 / (beta ** 2 + 1))
            * (
                (beta * self.x_err) ** 2
                - 2 * beta * self.corr * self.x_err * self.y_err
                + self.y_err ** 2
            )
        )  # intrinsic variance orthogonal to the regression line
        s_orth = np.sqrt(np.mean(V_orth))
        s = s_orth / np.sin(np.arctan(1 / beta))

        if verbose:
            print(f"BCES: beta={beta:.3}, alpha={alpha:.3}, s={s:.3}")

        return alpha, beta, s

    # ------------------------------------------------------- #

    def fit_lira(self, nmix, nsteps, lira_args={}):
        """
        Use the LIRA R package to perform the fit.

        Parameters
        ----------
        nmix : int
            Number of gaussians in the mixture model;
        nsteps : int
            Number of steps to perform in the MCMC;
        lira_args : dict
            Arguments to pass to LIRA. The syntax is the same
            as for the R function.

        Returns
        -------
        chains : DataFrame
            The MCMC chains in the parameter space

        Examples
        --------
        d.fit_lira(3, 10_000, lira_args={"sigma.YIZ.0": "prec.dgamma"} )
        d.fit_lira(3, 10_000, lira_args={
            "sigma.XIZ.0": 0.0, "sigma.YIZ.0": "dunif(0.0, 1.0)"
        })
        """

        # Create R function to be ran
        with open("r_script.R", "r") as f:
            rfunc_str = f.read()

        rfunc = ro.r(rfunc_str)
        lira_py = ro.functions.wrap_r_function(rfunc, "lira_py")

        # Format data to give to LIRA
        d = self.to_table().to_pandas()
        with ro.conversion.localconverter(
            ro.default_converter + pandas2ri.converter
        ):
            d_r = ro.conversion.py2rpy(d)

        # Run the fit
        chains_r = lira_py(d_r, int(nsteps), int(nmix), **lira_args)

        # Fetch chains and format them
        with ro.conversion.localconverter(
            ro.default_converter + pandas2ri.converter
        ):
            chains = ro.conversion.rpy2py(chains_r)
        self.lira_chains = chains

        return chains

    # ------------------------------------------------------- #

    def plot_lira_results(self, nmix=1, truth=None, chains_file=None):
        def rename_params(cc):
            for chain in cc.chains:
                for i, p in enumerate(chain.parameters):
                    new_p = r"$" + p + r"$"
                    new_p = new_p.replace("alpha", r"\alpha")
                    new_p = new_p.replace("beta", r"\beta")
                    new_p = new_p.replace("sigma", r"\sigma")
                    new_p = new_p.replace(".YIZ", r"_{Y|Z}")
                    new_p = new_p.replace(".XIZ", r"_{X|Z}")
                    new_p = new_p.replace(".0", "")
                    chain.parameters[i] = new_p

        if (not hasattr(self, "lira_chains")) and (chains_file is not None):
            self.lira_chains = pd.read_csv(chains_file)
        chains = self.lira_chains
        params_toplot = ["alpha.YIZ", "beta.YIZ", "sigma.YIZ.0"]
        cc = ChainConsumer()
        cc.add_chain(
            chains[params_toplot],
            shade=True,
            shade_alpha=0.3,
            shade_gradient=0.0,
            name="LIRA",
        )
        rename_params(cc)
        cc.configure(cmap="RdBu_r", sigmas=[1], summary=False)

        # Corner plot
        fig_corner = cc.plotter.plot(figsize=(8, 8), truth=truth)
        fig_corner.subplots_adjust(hspace=0, wspace=0)
        fig_corner.align_labels()
        fig_corner.savefig(self.path_to_results + "/corner.pdf")

        # Trace plot
        fig_walk = cc.plotter.plot_walks(truth=truth)
        fig_walk.align_labels()
        fig_walk.savefig(self.path_to_results + "/walk.pdf")
        plt.close(fig_corner)
        plt.close(fig_walk)

        # Latent distribution + model
        x = np.linspace(1.5 * self.x_obs.min(), 1.5 * self.x_obs.max(), 100)
        y = 0.0 * x
        cols = [""]
        for i in range(nmix - 1):
            cols.append(f".mixture.{i+2}.")
        for i, col in enumerate(cols):
            mu = np.median(chains["mu.Z.0" + col])
            sigma = np.median(chains["sigma.Z.0" + col])
            pi = np.median(chains[f"pi.{i+1}."])
            y += pi * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
        y /= np.trapz(y, x)
        fig_hist, ax = plt.subplots()
        ax.hist(
            self.x_obs,
            bins=self.n_pts // 5,
            density=True,
            alpha=0.3,
            color="tab:blue",
            label="Data",
        )
        ax.plot(x, y, color="tab:blue", label="Modeled distribution")
        ax.set_xlabel(self.xlabel)
        ax.legend(frameon=False)
        fig_hist.savefig(self.path_to_results + "/hist.pdf")
        plt.close(fig_hist)

    # ------------------------------------------------------- #

    def plot_data(self, style="errb"):
        """
        Make (x, y) plot

        Args:
            style (str): how to plot data points. Can be
                "errb" (error bars), "ellipse" (error ellipses),
                or "points" (no uncertainty representation)

        Returns:
            (tuple) fig, ax
        """

        kw = _kwargs
        fig, ax = plt.subplots(figsize=(5, 5))

        if style == "ellipse":
            from matplotlib.patches import Ellipse

            def eigsorted(cov):
                vals, vecs = np.linalg.eigh(cov)
                order = vals.argsort()[::-1]
                return vals[order], vecs[:, order]

            for x, y, cov in zip(self.x_obs, self.y_obs, self.covmats):
                eigval, eigvec = eigsorted(cov)
                theta = np.degrees(np.arctan2(*eigvec[:, 0][::-1]))
                w, h = 2 * np.sqrt(eigval)
                ell = Ellipse(
                    xy=(x, y),
                    width=w,
                    height=h,
                    angle=theta,
                    ec="k",
                    fc="#00000000",
                    lw=0.75,
                    ls=":",
                )
                ax.add_patch(ell)

            ax.plot(self.x_obs, self.y_obs, "kD", ms=2)
        elif style == "errb":
            ax.errorbar(
                self.x_obs,
                self.y_obs,
                xerr=self.x_err,
                yerr=self.y_err,
                mfc="w",
                **kw,
            )
        else:
            ax.plot(self.x_obs, self.y_obs, "o", mfc="w", mec="k")

        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        fig.subplots_adjust(left=0.15, bottom=0.15, right=0.85, top=0.85)

        return fig, ax

    # ------------------------------------------------------- #

    def plot_alphabeta(
        self, ax, alpha, beta, label="", addeq=True, setlims=False, **kwargs
    ):
        """
        Plots a line with $y = alpha + beta x$ on an axis.
        You should put the `xlims` before using this function,
            as it uses `ax.get_xlims` to span across the whole
            subplot.
        `label` will be added the line equation if `addeq=True`
        """
        x_span = np.array(ax.get_xlim())
        ax.set_xlim(*x_span)
        if isinstance(alpha, float):
            y = alpha + beta * x_span
            if setlims:
                ax.set_ylim(*y)
            if addeq:
                if label is not None and label != "":
                    label += ": "
                sign = "+" if alpha >= 0.0 else "-"
                label = (
                    label + f"$y = {beta:.2f} x {sign} {np.abs(alpha):.2f}$"
                )
            ax.plot(x_span, y, label=label, **kwargs)
        else:
            x = np.linspace(*x_span, 100)
            alphas, xs = np.meshgrid(alpha, x)
            betas, _ = np.meshgrid(beta, x)
            ys = alphas + betas * xs
            pct = np.percentile(ys, (16, 84), axis=1)
            ax.fill_between(
                x, pct[0], pct[1], alpha=0.3, label=label, **kwargs
            )

