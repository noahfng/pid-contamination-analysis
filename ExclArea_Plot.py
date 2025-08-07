import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erfc
from matplotlib.backends.backend_pdf import PdfPages

sigKa = np.array([13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 1, 0.75, 0.5, 0.25, 0])
AreaKa =     np.array([0.86, 1.07, 1.52, 2.50, 3.23, 4.12, 5.13, 6.50, 9.23, 14.34, 24.94, 29.88, 35.77, 44.57, 54.62, 67.28, 87.89, 112.10, 141.43, 175.75, 219.02, 263.88, 308.81])
AreaKa_err = np.array([1.31, 2.18, 2.11, 2.32, 2.29, 1.95, 2.21, 2.30, 2.30, 2.37, 2.30, 2.52, 2.63, 2.74, 2.78, 2.97, 3.21, 3.53, 3.75, 4.15, 4.58, 4.90, 5.42])
chi2Ka = np.array([1.66247, 1.66475, 1.65947, 1.65515, 1.65082, 1.6473, 1.63918, 1.63777, 1.63336, 1.62737, 1.6334, 1.63111, 1.61772, 1.61145, 1.60784, 1.61561, 1.6073, 1.64207, 1.70633, 1.80732, 1.90965, 2.00597, 2.12548])
#const Double_t pStart = 0.45, pEnd = 0.55, step = 0.1;
#const Bool_t FitKaonExclComp = true;
#const Bool_t FitProtonExclComp = false;
#manualNGauss = 4;
#manualMeans = {-6, -5, 0, 1.5};
#manualSigmas = {1, 1, 1, 1};
#manualAmps = {1e4, 1e4, 1e2, 1e1};
#par[offM + 2] = -0.271677;
#par[offS + 2] = 1.08472;
#par[offA1 + 2] = 78.0623;
#par[offA2 + 2] = par[offA1 + 2];
#//par[offA2 + 3] = 0;
#std::cout << "Mean: " << par[offM + 2] << ", Sigma: " << par[offS + 2] << ", Amp2: "<< par[offA2 + 2] << std::endl;

sigPr = np.array([13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 1, 0.75, 0.5, 0.25, 0])
AreaPr = np.array([0, 0, 0, 0.09, 0.08, 0.16, 0.65, 0.41, 0.82, 1.51, 4.92, 6.89, 9.30, 12.81, 16.63, 23.85, 34.82, 49.95, 69.80, 93.46, 123.34, 157.39, 189.79])
AreaPr_err = np.array([0, 0, 0.03, 1.49, 1.55, 1.50, 0.52, 1.50, 1.56, 1.66, 1.64, 1.67, 1.67, 1.97, 1.84, 1.97, 2.10, 2.36, 2.69, 3.08, 3.44, 3.85, 3.99])
chi2Pr = np.array([1.43889, 1.43374, 1.43729, 1.4323, 1.43458, 1.43238, 1.43212, 1.43347, 1.43785, 1.42486, 1.41977, 1.43227, 1.41255, 1.40898, 1.4024, 1.39478, 1.39959, 1.38814, 1.38673, 1.36941, 1.35767, 1.36301, 1.37432])
#const Double_t pStart = 0.9, pEnd = 1.0, step = 0.1;
#const Bool_t FitKaonExclComp = false;
#const Bool_t FitProtonExclComp = true;
#manualNGauss = 3;
#manualMeans = {-5, -0.5, 1};
#manualSigmas = {1, 1, 1};
#manualAmps = {1e3, 1e2, 1e2};
#par[offM + 1] = -1;
#par[offS + 1] = 1.46314;
#par[offA1 + 1] = 22.3524;
#par[offA2 + 1] = par[offA1 + 1];
#//par[offA2 + 2] = 0;
#std::cout << "Mean: " << par[offM + 1] << ", Sigma: " << par[offS + 1] << ", Amp2: "<< par[offA2 + 1] << std::endl;


def tail_area(n, A, sigma0):
    return A * erfc(n / (np.sqrt(2) * sigma0))

poptKa, pcovKa = curve_fit(tail_area, sigKa, AreaKa, sigma=AreaKa_err, absolute_sigma=True, p0=[AreaKa.max(), 1.0])
Aka, sigma0_ka = poptKa

poptPr, pcovPr = curve_fit(tail_area, sigPr, AreaPr, sigma=AreaPr_err, absolute_sigma=True, p0=[AreaPr.max(), 1.0])
Apr, sigma0_pr = poptPr


with PdfPages('PeakAreas_Kaon_Proton.pdf') as pdf:
    fig1, ax1 = plt.subplots()
    n_fit_ka = np.linspace(sigKa.min(), sigKa.max(), 200)
    ax1.errorbar(sigKa, AreaKa, yerr=AreaKa_err, fmt='o',markersize=3.5, color="blue", label='Kaon-Area ± error', capsize=3)
    ax1.plot(n_fit_ka, tail_area(n_fit_ka, Aka, sigma0_ka), '-', label=f'Fit: A={Aka:.1f}, σ₀={sigma0_ka:.2f}')
    ax1.set_xlabel('Sigma')
    ax1.set_ylabel('Kaon-Peak-area', rotation=90)
    ax1.tick_params(axis='y')

    ax1b = ax1.twinx()
    ax1b.plot(sigKa, chi2Ka, marker='o', markersize=3.5, color="red", linestyle='', label='χ² / NDF')
    ax1b.set_ylabel('χ² / NDF', rotation=90)
    ax1b.set_ylim(0,5)
    ax1b.tick_params(axis='y')
    plt.title('Kaon-Peak-Area vs Excl-Sigma')
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1b.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1+ labels2, loc= "best")
    pdf.savefig(fig1)
    plt.close(fig1)

    # Seite 2: Proton
    fig2, ax2 = plt.subplots()
    n_fit_pr = np.linspace(sigKa.min(), sigKa.max(), 200)
    ax2.errorbar(sigPr, AreaPr, yerr=AreaPr_err, fmt='o',markersize=3.5, color="blue", label='Proton-Area ± error', capsize=3)
    ax2.plot(n_fit_pr, tail_area(n_fit_pr, Apr, sigma0_pr), '-', label=f'Fit: A={Apr:.1f}, σ₀={sigma0_pr:.2f}')
    ax2.set_xlabel('Sigma')
    ax2.set_ylabel('Proton-Peak-area', rotation=90)
    ax2.tick_params(axis='y')

    ax2b = ax2.twinx()
    ax2b.plot(sigPr, chi2Pr, marker='o', markersize=3.5, color="red", linestyle='', label='χ² / NDF')
    ax2b.set_ylabel('χ² / NDF', rotation=90)
    ax2b.set_ylim(0,5)
    ax2b.tick_params(axis='y')
    plt.title('Proton-Peak-Area vs Excl-Sigma')
    lines3, labels3 = ax2.get_legend_handles_labels()
    lines4, labels4 = ax2b.get_legend_handles_labels()
    ax2.legend(lines3 + lines4, labels3+ labels4, loc= "best")
    pdf.savefig(fig2)
    plt.close(fig2)
