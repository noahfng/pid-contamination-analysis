import numpy as np
import matplotlib.pyplot as plt

# momentum-bin centers [GeV/c]
P_CENTERS = np.array([0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30], dtype=float)

# half-width of each momentum bin
P_ERR = np.array([0.00] * 10, dtype=float)

# nSigma windows
NSIGMA_WINDOWS = [0.5, 1.0, 1.5, 2.0, 3.0]

DATA = {
    0.5: {
        "yield": [65013.072547, 106357.005285, 78386.967561, 70327.284774, 53833.739782, 41120.299542, 15718.021611, 29261.752919, 19544.604272, 12633.957183],
        "yield_err": [299.693662, 618.028111, 829.076639, 299.393663, 308.812943, 309.850640, 504.484679, 380.004971, 296.458288, 187.999573],

        "no_veto": {
            "cont": [0.024512, 0.581808, 0.702450, 0.174662, 0.126509, 0.550888, 0.862739, 0.762053, 0.456986, 0.259749],
            "cont_err": [0.000679, 0.001424, 0.003862, 0.003628, 0.001975, 0.002305, 0.003519, 0.002446, 0.012805, 0.029945],
        },
        "veto": {
            "cont": [0.014914, 0.202886, 0.280829, 0.042942, 0.032283, 0.074491, 0.511166, 0.104141, 0.070094, 0.108948],
            "cont_err": [0.000501, 0.002091, 0.004973, 0.001054, 0.000976, 0.001280, 0.006596, 0.009475, 0.007284, 0.015002],
        },
    },

    1.0: {
        "yield": [115907.163169, 189616.307657, 139750.525294, 125381.492556, 95976.329322, 73310.444837, 28022.537994, 52168.689117, 34844.678888, 22524.179819],
        "yield_err": [534.302424, 1101.838173, 1478.101519, 533.767576, 550.560537, 552.410574, 899.409700, 677.483719, 528.534305, 335.171009],

        "no_veto": {
            "cont": [0.070704, 0.602438, 0.739095, 0.299554, 0.233625, 0.601531, 0.870461, 0.771062, 0.559043, 0.352306],
            "cont_err": [0.001913, 0.001382, 0.003094, 0.004124, 0.003520, 0.001868, 0.003333, 0.002357, 0.009955, 0.029426],
        },
        "veto": {
            "cont": [0.056814, 0.224391, 0.321325, 0.083913, 0.080451, 0.126999, 0.541820, 0.192547, 0.123200, 0.160454],
            "cont_err": [0.001728, 0.002145, 0.004990, 0.001554, 0.002131, 0.001987, 0.006123, 0.007913, 0.009790, 0.017140],
        },
    },

    1.5: {
        "yield": [147095.125923, 240637.712884, 177354.190664, 159118.780338, 121801.361203, 93036.606366, 35562.761112, 66206.088434, 44220.583860, 28584.920692],
        "yield_err": [678.070968, 1398.317587, 1875.824782, 677.392205, 698.703767, 701.051607, 1141.420249, 859.779070, 670.750781, 425.357851],

        "no_veto": {
            "cont": [0.224411, 0.634332, 0.780998, 0.481543, 0.435840, 0.679037, 0.882671, 0.798643, 0.668721, 0.489977],
            "cont_err": [0.005279, 0.001315, 0.002335, 0.003771, 0.004963, 0.001592, 0.003042, 0.002138, 0.007066, 0.025110],
        },
        "veto": {
            "cont": [0.207041, 0.272389, 0.384550, 0.170489, 0.208668, 0.278592, 0.605559, 0.421342, 0.241944, 0.256570],
            "cont_err": [0.004822, 0.002242, 0.005030, 0.002253, 0.004280, 0.004753, 0.005712, 0.005382, 0.012758, 0.018617],
        },
    },

    2.0: {
        "yield": [162055.162606, 265111.324700, 195391.669346, 175301.660469, 134188.942506, 102498.721683, 39179.605704, 72939.455738, 48717.956240, 31492.096979],
        "yield_err": [747.032916, 1540.530881, 2066.602059, 746.285122, 769.764136, 772.350759, 1257.506276, 947.221304, 738.968243, 468.618081],

        "no_veto": {
            "cont": [0.516668, 0.678026, 0.819041, 0.653764, 0.662165, 0.775363, 0.899939, 0.856386, 0.762448, 0.637801],
            "cont_err": [0.007318, 0.001279, 0.001836, 0.002669, 0.004387, 0.002030, 0.002673, 0.001767, 0.005190, 0.017382],
        },
        "veto": {
            "cont": [0.504618, 0.371646, 0.483687, 0.314027, 0.439952, 0.545120, 0.709176, 0.705617, 0.444838, 0.399979],
            "cont_err": [0.006627, 0.002771, 0.005368, 0.002926, 0.006031, 0.006977, 0.006677, 0.003629, 0.014406, 0.016987],
        },
    },

    3.0: {
        "yield": [169321.834886, 276999.110831, 204153.175025, 183162.315424, 140206.073051, 107094.839497, 40936.447944, 76210.114400, 50902.505109, 32904.225692],
        "yield_err": [780.530420, 1609.609415, 2159.270010, 779.749093, 804.280924, 806.983533, 1313.893779, 989.695402, 772.104121, 489.631260],

        "no_veto": {
            "cont": [0.910391, 0.815535, 0.898586, 0.859162, 0.910063, 0.931014, 0.947510, 0.957560, 0.903452, 0.856812],
            "cont_err": [0.002016, 0.001898, 0.001754, 0.001189, 0.001599, 0.001463, 0.002113, 0.000682, 0.003337, 0.005064],
        },
        "veto": {
            "cont": [0.909354, 0.726437, 0.812823, 0.715934, 0.858234, 0.904165, 0.906878, 0.943818, 0.836482, 0.746584],
            "cont_err": [0.001856, 0.003831, 0.005240, 0.003360, 0.003142, 0.002543, 0.004767, 0.000733, 0.007600, 0.007695],
        },
    },
}


def arr(x):
    return np.array(x, dtype=float)

def validate():
    for ns in NSIGMA_WINDOWS:
        if ns not in DATA:
            raise ValueError(f"Missing window {ns}")

        for key in ["yield", "yield_err"]:
            if key not in DATA[ns]:
                raise ValueError(f"Missing '{key}' for nsigma={ns}")
            if len(DATA[ns][key]) != len(P_CENTERS):
                raise ValueError(f"Wrong length for {key} at nsigma={ns}")

        for cfg in ["no_veto", "veto"]:
            if cfg not in DATA[ns]:
                raise ValueError(f"Missing '{cfg}' for nsigma={ns}")
            for key in ["cont", "cont_err"]:
                if key not in DATA[ns][cfg]:
                    raise ValueError(f"Missing '{key}' in {cfg} for nsigma={ns}")
                if len(DATA[ns][cfg][key]) != len(P_CENTERS):
                    raise ValueError(f"Wrong length for {cfg}/{key} at nsigma={ns}")

def weighted_constant_fit(y, yerr):
    y = arr(y)
    yerr = arr(yerr)
    mask = np.isfinite(y) & np.isfinite(yerr) & (yerr > 0)
    if np.sum(mask) == 0:
        return np.nan, np.nan
    w = 1.0 / yerr[mask]**2
    c = np.sum(w * y[mask]) / np.sum(w)
    c_err = np.sqrt(1.0 / np.sum(w))
    return c, c_err

def yield_weighted_average(cont, cont_err, yld, yld_err):
    cont = arr(cont)
    cont_err = arr(cont_err)
    yld = arr(yld)
    yld_err = arr(yld_err)

    mask = np.isfinite(cont) & np.isfinite(cont_err) & np.isfinite(yld) & np.isfinite(yld_err) & (yld >= 0)
    if np.sum(mask) == 0:
        return np.nan, np.nan

    cont = cont[mask]
    cont_err = cont_err[mask]
    yld = yld[mask]
    yld_err = yld_err[mask]

    S = np.sum(yld)
    if S <= 0:
        return np.nan, np.nan

    C = np.sum(cont * yld) / S

    var = 0.0
    for ci, sci, Yi, sYi in zip(cont, cont_err, yld, yld_err):
        dC_dci = Yi / S
        dC_dYi = (ci - C) / S
        var += (dC_dci * sci)**2 + (dC_dYi * sYi)**2

    return C, np.sqrt(var)

def plot_cont_vs_momentum(nsigma=2.0, show_no_veto=True, show_veto=True, title_size=20, label_size=18, tick_size=16, legend_size=14):
    plt.figure(figsize=(8, 5))

    ymax_candidates = []

    if show_no_veto:
        y_no_veto = arr(DATA[nsigma]["no_veto"]["cont"])
        yerr_no_veto = arr(DATA[nsigma]["no_veto"]["cont_err"])
        plt.errorbar(
            P_CENTERS,
            y_no_veto,
            yerr=yerr_no_veto,
            fmt="o-",
            capsize=3,
            label="No Exclusion"
        )
        ymax_candidates.append(np.max(y_no_veto + yerr_no_veto))

    if show_veto:
        y_veto = arr(DATA[nsigma]["veto"]["cont"])
        yerr_veto = arr(DATA[nsigma]["veto"]["cont_err"])
        plt.errorbar(
            P_CENTERS,
            y_veto,
            yerr=yerr_veto,
            fmt="o-",
            capsize=3,
            label="Exclusion"
        )
        ymax_candidates.append(np.max(y_veto + yerr_veto))

    ymax = max(ymax_candidates) if ymax_candidates else 1.0
    y_position = 0.7 * ymax

    plt.xlabel(r"$p$ (GeV/$c$)", fontsize=label_size)
    plt.ylabel("Electron contamination", fontsize=label_size)
    plt.title(rf"Contamination vs momentum ($\sigma_{{\mathrm{{window}}}} = {nsigma}$)",fontsize=title_size)

    plt.axvline(x=0.65, linestyle="--", linewidth=1.5, color="black")
    plt.text(0.66, y_position, "K → K+p exclusion", rotation=90, fontsize=16)

    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)
    plt.ylim(0, 1.05 * ymax)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=legend_size)
    plt.tight_layout()
    plt.show()
    
def plot_cont_vs_momentum_all_windows(cfg="no_Exclusion", title_size=20, label_size=18, tick_size=16, legend_size=14, marker_size=6, line_width=1.5):
    cfg_map = {
        "no_Exclusion": "no_veto",
        "Exclusion": "veto",
        "no_veto": "no_veto",
        "veto": "veto",
    }

    if cfg not in cfg_map:
        raise ValueError(f"Unknown cfg='{cfg}'. Use one of: {list(cfg_map.keys())}")

    data_key = cfg_map[cfg]

    plt.figure(figsize=(9, 6))

    ymax = 0.0

    for ns in NSIGMA_WINDOWS:
        y = arr(DATA[ns][data_key]["cont"])
        yerr = arr(DATA[ns][data_key]["cont_err"])

        plt.errorbar(
            P_CENTERS,
            y,
            yerr=yerr,
            fmt="o-",
            capsize=3,
            markersize=marker_size,
            linewidth=line_width,
            label=rf"$\sigma_{{\mathrm{{window}}}} = {ns}$"
        )

        ymax = max(ymax, np.max(y + yerr))

    y_position = 0.72 * ymax

    plt.xlabel(r"$p$ (GeV/$c$)", fontsize=label_size)
    plt.ylabel("Electron contamination", fontsize=label_size)
    plt.title(f"Contamination vs momentum ({cfg.replace('_', ' ')})", fontsize=title_size)

    plt.axvline(x=0.65, linestyle="--", linewidth=1.5, color="black")
    plt.text(0.66, y_position, "K → K+p Exclusion", rotation=90, fontsize=16)

    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)
    plt.ylim(0, 1.05 * ymax)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=legend_size)
    plt.tight_layout()
    plt.show()

def plot_avg_cont_vs_nsigma(mode="yield_weighted", show_no_veto=True, show_veto=True, title_size=20, label_size=18, tick_size=16, legend_size=14, marker_size=6, line_width=1.5):
    plt.figure(figsize=(8, 5))

    if show_no_veto:
        vals, errs = [], []
        for ns in NSIGMA_WINDOWS:
            if mode == "yield_weighted":
                v, e = yield_weighted_average(
                    DATA[ns]["no_veto"]["cont"],
                    DATA[ns]["no_veto"]["cont_err"],
                    DATA[ns]["yield"],
                    DATA[ns]["yield_err"]
                )
            else:
                v, e = weighted_constant_fit(
                    DATA[ns]["no_veto"]["cont"],
                    DATA[ns]["no_veto"]["cont_err"]
                )
            vals.append(v)
            errs.append(e)

        plt.errorbar(
            NSIGMA_WINDOWS,
            vals,
            yerr=errs,
            fmt="o-",
            capsize=4,
            markersize=marker_size,
            linewidth=line_width,
            label="No Exclusion"
        )

    if show_veto:
        vals, errs = [], []
        for ns in NSIGMA_WINDOWS:
            if mode == "yield_weighted":
                v, e = yield_weighted_average(
                    DATA[ns]["veto"]["cont"],
                    DATA[ns]["veto"]["cont_err"],
                    DATA[ns]["yield"],
                    DATA[ns]["yield_err"]
                )
            else:
                v, e = weighted_constant_fit(
                    DATA[ns]["veto"]["cont"],
                    DATA[ns]["veto"]["cont_err"]
                )
            vals.append(v)
            errs.append(e)

        plt.errorbar(
            NSIGMA_WINDOWS,
            vals,
            yerr=errs,
            fmt="o-",
            capsize=4,
            markersize=marker_size,
            linewidth=line_width,
            label="Exclusion"
        )

    plt.xlabel(r"$\sigma_{\mathrm{window}}$", fontsize=label_size)
    plt.ylabel("Average electron contamination", fontsize=label_size)
    plt.title(f"Yield-weighted average contamination vs selection window", fontsize=title_size)

    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)

    plt.ylim(bottom=0)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=legend_size)
    plt.tight_layout()
    plt.show()

validate()

plot_cont_vs_momentum_all_windows(cfg="no_Exclusion")
plot_cont_vs_momentum_all_windows(cfg="Exclusion")

plot_cont_vs_momentum(nsigma=0.5, show_no_veto=True, show_veto=True)
plot_cont_vs_momentum(nsigma=1.0, show_no_veto=True, show_veto=True)
plot_cont_vs_momentum(nsigma=1.5, show_no_veto=True, show_veto=True)
plot_cont_vs_momentum(nsigma=2.0, show_no_veto=True, show_veto=True)
plot_cont_vs_momentum(nsigma=3.0, show_no_veto=True, show_veto=True)

plot_avg_cont_vs_nsigma(mode="yield_weighted", show_no_veto=True, show_veto=True)
