import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erfc
from matplotlib.backends.backend_pdf import PdfPages
from scipy.special import logsumexp

#-------------------------------0.45 < p < 0.55 GeV/c ---------------------------------
sigKa = np.array([13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 1, 0.75, 0.5, 0.25, 0])
AreaKa = np.array([0.86, 1.07, 1.52, 2.50, 3.23, 4.12, 5.13, 6.50, 9.23, 14.34, 24.94, 29.88, 35.77, 44.57, 54.62, 67.28, 87.89, 112.10, 141.43, 175.75, 219.02, 263.88, 308.81])
AreaKa_err = np.array([1.31, 2.18, 2.11, 2.32, 2.29, 1.95, 2.21, 2.30, 2.30, 2.37, 2.30, 2.52, 2.63, 2.74, 2.78, 2.97, 3.21, 3.53, 3.75, 4.15, 4.58, 4.90, 5.42])
chi2Ka = np.array([1.66247, 1.66475, 1.65947, 1.65515, 1.65082, 1.6473, 1.63918, 1.63777, 1.63336, 1.62737, 1.6334, 1.63111, 1.61772, 1.61145, 1.60784, 1.61561, 1.6073, 1.64207, 1.70633, 1.80732, 1.90965, 2.00597, 2.12548])
DNKa = np.array([4.52325, 4.52327, 4.52326, 4.52325, 4.52323, 4.52321, 4.52322, 4.52322, 4.52323, 4.52326, 4.52331, 4.52334, 4.52341, 4.52344, 4.52362, 4.52379, 4.5239, 4.52417, 4.5245, 4.52496, 4.52539, 4.52576, 4.52622])
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

#-------------------------------0.9 < p < 1.0 GeV/c ---------------------------------
sigPr = np.array([13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 1, 0.75, 0.5, 0.25, 0])
AreaPr = np.array([0, 0, 0, 0.09, 0.08, 0.16, 0.65, 0.41, 0.82, 1.51, 4.92, 6.89, 9.30, 12.81, 16.63, 23.85, 34.82, 49.95, 69.80, 93.46, 123.34, 157.39, 189.79])
AreaPr_err = np.array([0, 0, 0.03, 1.49, 1.55, 1.50, 0.52, 1.50, 1.56, 1.66, 1.64, 1.67, 1.67, 1.97, 1.84, 1.97, 2.10, 2.36, 2.69, 3.08, 3.44, 3.85, 3.99])
chi2Pr = np.array([1.43889, 1.43374, 1.43729, 1.4323, 1.43458, 1.43238, 1.43212, 1.43347, 1.43785, 1.42486, 1.41977, 1.43227, 1.41255, 1.40898, 1.4024, 1.39478, 1.39959, 1.38814, 1.38673, 1.36941, 1.35767, 1.36301, 1.37432])
DNPr = np.array([4.53047, 4.53047, 4.53051, 4.53045, 4.53048, 4.53051, 4.53066, 4.53052, 4.53056, 4.53063, 4.53048, 4.53041, 4.53036, 4.53039, 4.53053, 4.5305, 4.53064, 4.53064, 4.53062, 4.53055, 4.53043, 4.53043, 4.5305])
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

#-------------------------------0.45 < p < 0.55 GeV/c ---------------------------------
sigKaFixElArea = np.array([13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 1, 0.75, 0.5, 0.25, 0])
AreaKaElFixElArea = np.array([237.41, 237.61, 237.78, 236.91, 236.90, 236.86, 237.19, 237.65, 238.66, 239.38, 242.39, 242.86, 246.09, 247.90, 255.87, 261.23, 265.20, 276.97, 288.72, 307.53, 321.65, 333.30, 483.79])
AreaKaElFixElArea_err = np.array([5.61, 5.76, 5.66, 6.00, 6.04, 6.09, 6.14, 6.14, 6.26, 6.15, 6.38, 6.42, 6.58, 6.45, 6.94, 6.93, 7.20, 7.74, 8.18, 8.37, 8.93, 9.56, 9.36])
AreaKaFixElArea = np.array([0.44, 0.57, 1.07, 2.70, 3.41, 4.28, 5.03, 6.13, 7.91, 12.09, 20.26, 24.87, 28.25, 35.88, 41.01, 50.51, 68.29, 84.37, 105.54, 126.73, 159.79, 196.27, 94.81])
AreaKaFixElArea_err = np.array([2.40, 3.84, 4.07, 4.25, 4.37, 4.23, 4.37, 4.40, 4.46, 4.72, 5.23, 5.24, 5.37, 5.58, 6.19, 6.40, 6.77, 7.54, 8.81, 10.08, 11.69, 13.32, 6.51])
chi2KaFixElArea = np.array([1.50571, 1.50594, 1.49773, 1.49221, 1.48745, 1.48427, 1.47456, 1.4718, 1.466, 1.45689, 1.45845, 1.4442, 1.41731, 1.39587, 1.34463, 1.30091, 1.23513, 1.1822, 1.15611, 1.12202, 1.10199, 1.07952, 1.05626])
DNKaFixElArea = np.array([4.52145, 4.52145, 4.52143, 4.52141, 4.5214, 4.52139, 4.52138, 4.52138, 4.52137, 4.52138, 4.52139, 4.52136, 4.52132, 4.52127, 4.52116, 4.52107, 4.5209, 4.52077, 4.52071, 4.52063, 4.52058, 4.52052, 4.52047])
#const Double_t pStart = 0.45, pEnd = 0.55, step = 0.1;
#const Bool_t FitKaonExclComp = true;
#const Bool_t FitProtonExclComp = false;
#manualNGauss = 4;
#manualMeans = {-6, -5, 0, 1.5};
#manualSigmas = {1, 1, 1, 1};
#manualAmps = {1e4, 1e4, 1e2, 1e1};
#par[offA2 + 2] = par[offA1 + 2];

#-------------------------------0.9 < p < 1.0 GeV/c ---------------------------------
sigPrFixElArea = np.array([13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 1, 0.75, 0.5, 0.25, 0])
AreaPrElFixElArea = np.array([80.48, 80.84, 80.93, 81.27, 81.35, 81.59, 81.65, 81.96, 82.10, 82.75, 82.71, 82.54, 83.63, 84.00, 83.98, 84.91, 85.46, 88.80, 91.90, 94.40, 102.35, 103.66, 100.75])
AreaPrElFixElArea_err = np.array([3.22, 3.25, 3.28, 3.31, 3.29, 3.43, 3.30, 3.25, 4.21, 4.16, 4.78, 5.05, 0.00, 5.26, 5.16, 5.45, 5.94, 6.93, 8.12, 8.97, 9.43, 0.00, 9.58])
AreaPrFixElArea = np.array([0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 4.19, 6.31, 8.25, 11.61, 15.52, 22.26, 32.97, 46.37, 64.82, 87.68, 114.54, 149.05, 183.34])
AreaPrFixElArea_err = np.array([0.16, 0.16, 0.28, 0.20, 0.14, 0.03, 0.07, 0.39, 0.01, 0.02, 3.29, 3.24, 0.00, 3.41, 3.40, 3.60, 3.86, 3.93, 4.56, 5.17, 5.44, 0.00, 6.37])
chi2PrFixElArea = np.array([1.43852, 1.4335, 1.43703, 1.43201, 1.43436, 1.43222, 1.43182, 1.43171, 1.43322, 1.43751, 1.42459, 1.41911, 1.41129, 1.40775, 1.40187, 1.39383, 1.3986, 1.38589, 1.3829, 1.36369, 1.34542, 1.34465, 1.35345])
DNPrFixElArea = np.array([4.53059, 4.53054, 4.53056, 4.53051, 4.53053, 4.53053, 4.53053, 4.53055, 4.53057, 4.53062, 4.53055, 4.53052, 4.53048, 4.53049, 4.53057, 4.53054, 4.53066, 4.53063, 4.5306, 4.53052, 4.53037, 4.53035, 4.53042])
#const Double_t pStart = 0.9, pEnd = 1.0, step = 0.1;
#const Bool_t FitKaonExclComp = false;
#const Bool_t FitProtonExclComp = true;
#manualNGauss = 3;
#manualMeans = {-5, -0.5, 1};
#manualSigmas = {1, 1, 1};
#manualAmps = {1e3, 1e2, 1e2};
#par[offA2 + 1] = par[offA1 + 1];

##-------------------------------0.55 < p < 0.65 GeV/c ---------------------------------
#sigKa = np.array([13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 1, 0.75, 0.5, 0.25, 0])
#AreaKa = np.array([0.00, 0.00, 0.00, 0.00, 0.32, 1.25, 2.39, 3.78, 5.43, 8.42, 16.69, 21.44, 28.86, 37.45, 52.17, 72.01, 99.01, 137.38, 189.09, 252.84])
#AreaKa_err = np.array([0.00, 0.00, 0.00, 0.02, 4.04, 3.39, 4.60, 4.56, 4.43, 4.36, 4.63, 4.64, 4.65, 5.08, 4.64, 5.78, 6.82, 8.52, 11.54, 16.01])
#chi2Ka = np.array([1.20463, 1.20446, 1.20625, 1.20841, 1.20907, 1.20897, 1.20865, 1.20436, 1.20319, 1.204, 1.20481, 1.19952, 1.19478, 1.19312, 1.18671, 1.18384, 1.18051, 1.17878, 1.18036, 1.19297])
#DNKa = np.array([5.21149, 5.21149, 5.21153, 5.21155, 5.21156, 5.21155, 5.21155, 5.21154, 5.21154, 5.21154, 5.21154, 5.21152, 5.21149, 5.21148, 5.21142, 5.21143, 5.2114, 5.21139, 5.21138, 5.21141])
##const Bool_t FitKaonExclComp = true;
##const Bool_t FitProtonExclComp = false;
##const Double_t pStart = 0.55, pEnd = 0.65, step = 0.1;
##manualNGauss = 4;
##manualMeans = {-6, -5, -2, 0};
##manualSigmas = {1, 1, 1, 1};
##manualAmps = {3e3, 3e3, 2e2, 2e1};
##par[offM + 3] = -0.310652;
##par[offS + 3] = 1.12093;
##par[offA1 + 3] = 43.4094;
##par[offA2 + 3] = par[offA1 + 3];
##//par[offA2 + 2] = 0;
##std::cout << "Mean: " << par[offM + 3] << ", Sigma: " << par[offS + 3] << ", Amp2: "<< par[offA2 + 3] << std::endl;



def tail_area(n, A, sigma0):
    return A * erfc(n / (np.sqrt(2) * sigma0))

poptKa, pcovKa = curve_fit(tail_area, sigKa, AreaKa, sigma=AreaKa_err, absolute_sigma=True, p0=[AreaKa.max(), 1.0])
Aka, sigma0_ka = poptKa

poptPr, pcovPr = curve_fit(tail_area, sigPr, AreaPr, sigma=AreaPr_err, absolute_sigma=True, p0=[AreaPr.max(), 1.0])
Apr, sigma0_pr = poptPr

poptPrFixElArea, pcovPrFixElArea = curve_fit(tail_area, sigPrFixElArea, AreaPrFixElArea, sigma=AreaPrFixElArea_err, absolute_sigma=True, p0=[AreaPrFixElArea.max(), 1.0])
AprfixElArea, sigma0_prfixelarea = poptPrFixElArea

poptKaFixElArea, pcovKaFixElArea = curve_fit(tail_area, sigKaFixElArea, AreaKaFixElArea, sigma=AreaKaFixElArea_err, absolute_sigma=True, p0=[AreaKaFixElArea.max(), 1.0])
AkafixElArea, sigma0_kafixelarea = poptKaFixElArea

def dn_weights(D_over_N, N_ref=None, target_ratio=0.1, quantile=0.9):
    r = np.asarray(D_over_N, float)
    dr = r - r.min()
    if N_ref is None:
        pos = dr[dr > 0]
        if pos.size == 0:
            return np.ones_like(r) / r.size, 0.0
        dq = np.quantile(pos, quantile)
        N_ref = 2.0 * np.log(1.0/target_ratio) / dq
    logw = -0.5 * N_ref * dr
    logw -= logsumexp(logw)
    return np.exp(logw), N_ref


with PdfPages('PeakAreas_Kaon_Proton.pdf') as pdf:
    fig1, ax1 = plt.subplots()
    n_fit_ka = np.linspace(sigKa.min(), sigKa.max(), 200)
    ax1.errorbar(sigKa, AreaKa, yerr=AreaKa_err, fmt='o',markersize=3.5, color="blue", label='Kaon-Area ± error', capsize=3)
    ax1.plot(n_fit_ka, tail_area(n_fit_ka, Aka, sigma0_ka), '-', label=f'Fit: A={Aka:.1f}, $\sigma_{0}$={sigma0_ka:.2f}')
    ax1.set_xlabel('Sigma')
    ax1.set_ylabel('Kaon-Peak-area', rotation=90)
    ax1.tick_params(axis='y')

    ax1b = ax1.twinx()
    ax1b.plot(sigKa, chi2Ka, marker='o', markersize=3.5, color="red", linestyle='', label='$\chi^{2}$/NDF')
    ax1b.set_ylabel('$\chi^{2}$/NDF', rotation=90)
    ax1b.set_ylim(0,5)
    ax1b.tick_params(axis='y')
    plt.title('Kaon-Peak-Area vs Excl-Sigma')
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1b.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1+ labels2, loc= "best")
    pdf.savefig(fig1)
    plt.close(fig1)
    
    fig2, ax2 = plt.subplots()
    w_ka, Nref_ka = dn_weights(DNKa)  
    ax2.plot(sigKa, w_ka, 'o', markersize=3.5)
    ax2.set_xlabel('Sigma-Exclusion')
    ax2.set_ylabel('D/N-weight w($\sigma$)')
    ax2.set_title(f'Kaon: Poisson-Deviance/Number-of-Events-weights vs Excl-Sigma')
    ax2.set_ylim(0, 1)
    pdf.savefig(fig2, bbox_inches='tight', pad_inches=0.1); plt.close(fig2)

    fig3, ax3 = plt.subplots()
    n_fit_pr = np.linspace(sigPr.min(), sigPr.max(), 200)
    ax3.errorbar(sigPr, AreaPr, yerr=AreaPr_err, fmt='o',markersize=3.5, color="blue", label='Proton-Area ± error', capsize=3)
    ax3.plot(n_fit_pr, tail_area(n_fit_pr, Apr, sigma0_pr), '-', label=f'Fit: A={Apr:.1f}, $\sigma_{0}$={sigma0_pr:.2f}')
    ax3.set_xlabel('Sigma')
    ax3.set_ylabel('Proton-Peak-area', rotation=90)
    ax3.tick_params(axis='y')

    ax3b = ax3.twinx()
    ax3b.plot(sigPr, chi2Pr, marker='o', markersize=3.5, color="red", linestyle='', label='$\chi^{2}$/NDF')
    ax3b.set_ylabel('$\chi^{2}$/NDF', rotation=90)
    ax3b.set_ylim(0,5)
    ax3b.tick_params(axis='y')
    plt.title('Proton-Peak-Area vs Excl-Sigma')
    lines3, labels3 = ax3.get_legend_handles_labels()
    lines4, labels4 = ax3b.get_legend_handles_labels()
    ax3.legend(lines3 + lines4, labels3+ labels4, loc= "best")
    pdf.savefig(fig3)
    plt.close(fig3)

    fig4, ax4 = plt.subplots()
    w_pr, Nref_pr = dn_weights(DNPr)
    ax4.plot(sigPr, w_pr, 'o', markersize=3.5)
    ax4.set_xlabel('Sigma-Exclusion')
    ax4.set_ylabel('D/N-weight w($\sigma$)')
    ax4.set_title(f'Proton: Poisson-Deviance/Number-of-Events-weights vs Excl-Sigma')
    ax4.set_ylim(0, 1)
    pdf.savefig(fig4, bbox_inches='tight', pad_inches=0.1); plt.close(fig4)
    
    fig5, ax5 = plt.subplots()
    n_fit_kafixelarea = np.linspace(sigKaFixElArea.min(), sigKaFixElArea.max(), 200)
    ax5.errorbar(sigKaFixElArea, AreaKaFixElArea, yerr=AreaKaFixElArea_err, fmt='o',markersize=3.5, color="blue", label='Kaon-Area ± error', capsize=3)
    ax5.errorbar(sigKaFixElArea, AreaKaElFixElArea, yerr=AreaKaElFixElArea_err, fmt='o', markersize=3.5, color="green", capsize=3, label='Electron-Area ± error')
    ax5.plot(n_fit_kafixelarea, tail_area(n_fit_kafixelarea, AkafixElArea, sigma0_kafixelarea), '-', label=f'Fit: A={AkafixElArea:.1f}, $\sigma_{0}$={sigma0_prfixelarea:.2f}')
    ax5.set_xlabel('Sigma')
    ax5.set_ylabel('KaonFixElArea-Peak-area', rotation=90)
    ax5.tick_params(axis='y')

    ax5b = ax5.twinx()
    ax5b.plot(sigKaFixElArea, chi2KaFixElArea, marker='o', markersize=3.5, color="red", linestyle='', label='$\chi^{2}$/NDF')
    ax5b.set_ylabel('$\chi^{2}$/NDF', rotation=90)
    ax5b.set_ylim(0,5)
    ax5b.tick_params(axis='y')
    plt.title('KaonFixElArea-Area vs Excl-Sigma')
    lines5, labels5 = ax5.get_legend_handles_labels()
    lines6, labels6 = ax5b.get_legend_handles_labels()
    ax5.legend(lines5 + lines6, labels5 + labels6, loc= "best")
    pdf.savefig(fig5)
    plt.close(fig5)
    
    fig6, ax6 = plt.subplots()
    w_pr, Nref_pr = dn_weights(DNKaFixElArea)
    ax6.plot(sigKaFixElArea, w_pr, 'o', markersize=3.5)
    ax6.set_xlabel('Sigma-Exclusion')
    ax6.set_ylabel('D/N-weight w($\sigma$)')
    ax6.set_title(f'KaonFixEl: Poisson-Deviance/Number-of-Events-weights vs Excl-Sigma')
    ax6.set_ylim(0, 1)
    pdf.savefig(fig6, bbox_inches='tight', pad_inches=0.1); plt.close(fig6)
    
    fig7, ax7 = plt.subplots()
    n_fit_prfixelarea = np.linspace(sigPrFixElArea.min(), sigPrFixElArea.max(), 200)
    ax7.errorbar(sigPrFixElArea, AreaPrFixElArea, yerr=AreaPrFixElArea_err, fmt='o',markersize=3.5, color="blue", label='Proton-Area ± error', capsize=3)
    ax7.errorbar(sigPrFixElArea, AreaPrElFixElArea, yerr=AreaPrElFixElArea_err, fmt='o', markersize=3.5, color="green", capsize=3, label='Electron-Area ± error')
    ax7.plot(n_fit_prfixelarea, tail_area(n_fit_prfixelarea, AprfixElArea, sigma0_prfixelarea), '-', label=f'Fit: A={AprfixElArea:.1f}, $\sigma_{0}$={sigma0_prfixelarea:.2f}')
    ax7.set_xlabel('Sigma')
    ax7.set_ylabel('ProtonFixElArea-Peak-area', rotation=90)
    ax7.tick_params(axis='y')

    ax7b = ax7.twinx()
    ax7b.plot(sigPrFixElArea, chi2PrFixElArea, marker='o', markersize=3.5, color="red", linestyle='', label='$\chi^{2}$/NDF')
    ax7b.set_ylabel('$\chi^{2}$/NDF', rotation=90)
    ax7b.set_ylim(0,5)
    ax7b.tick_params(axis='y')
    plt.title('ProtonFixElArea-Area vs Excl-Sigma')
    lines7, labels7 = ax7.get_legend_handles_labels()
    lines8, labels8 = ax7b.get_legend_handles_labels()
    ax7.legend(lines7 + lines8, labels7 + labels8, loc= "best")
    pdf.savefig(fig7)
    plt.close(fig7)
    
    fig8, ax8 = plt.subplots()
    w_pr, Nref_pr = dn_weights(DNPrFixElArea)
    ax8.plot(sigPrFixElArea, w_pr, 'o', markersize=3.5)
    ax8.set_xlabel('Sigma-Exclusion')
    ax8.set_ylabel('D/N-weight w($\sigma$)')
    ax8.set_title(f'ProtonFixEl: Poisson-Deviance/Number-of-Events-weights vs Excl-Sigma')
    ax8.set_ylim(0, 1)
    pdf.savefig(fig8, bbox_inches='tight', pad_inches=0.1); plt.close(fig8)