"""Order-of-magnitude estimates quoted in the text that are not power laws of
the solver: X-ray heating of the belt (Section 4.2), belt dust emission
(Section 5.1), flare variability (Section 5.2), metal-line cooling and H2
formation (Sections 3.2, 5.5), and the zeta propagation (Section 3.3).

    python scripts/lrd_estimates.py
"""
import json, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = json.load(open(os.path.join(ROOT, 'scripts', 'lrd_results.json')))
pc, yr, Msun, kB, mp, c, sigSB = 3.086e18, 3.156e7, 1.989e33, 1.381e-16, 1.673e-24, 2.998e10, 5.670e-5
T3 = 0.9  # window center

def val(key):
    """evaluate a power law from lrd_results.json at m6=1, T3=0.9, x=1, zeta=1"""
    for blk in R.values():
        if isinstance(blk, dict) and key in blk:
            d = blk[key]
            if isinstance(d, dict) and 'prefactor' in d:
                e = d.get('exponents', {})
                return d['prefactor'] * T3 ** e.get('lT3', e.get('T3', 0.0))
    raise KeyError(key)

# --- fiducial belt (Table 1) ---
try:
    RMC = val('R_MC[pc]'); nMC = val('n_MC[cm^-3]'); rb = val('r_b[pc]'); NMC = val('N_MC'); Mdisk = val('M_disk[Msun]')
    NH = val('N_H[cm^-2]'); Rudr = val('R_UDR[pc]'); tc = val('tau_c[yr]'); trad = val('tau_rad[yr]')
except KeyError:
    RMC, nMC, rb, NMC, Mdisk, NH, Rudr, tc, trad = 0.52, 1.9e5, 4.9, 77, 2.8e5, 5.9e22, 0.078, 2.4e5, 31
Geq, Eudr, LEdd, fOm = 1.4e-2, 1.8e50, 1.3e44, 0.22
LH2 = Geq / yr * Eudr

print('--- Section 4.2: X-ray heating ---')
LLTE = 9.5e-22 * T3**3.76
tcool = 1.5 * kB * 1e3 * T3 * 1.16 / LLTE / yr          # per H2, incl. He
print(f'cooling time at T3={T3}: {tcool:.0f} yr')
print(f'strike interval per cloud: tau_c (R_UDR/R_MC)^2 = {tc*(Rudr/RMC)**2:.2e} yr')
Msw = 4*np.pi/3*(Rudr*pc)**3 * nMC*1.32*mp / Msun
print(f'mass swept per remnant at R_UDR: {Msw:.0f} Msun')
print(f'remnant interval 1/Geq = {1/Geq:.0f} yr; radiative time {trad:.0f} yr')
sig1keV = 6.3e-18*(13.6/1000)**3 + 0.079*7.4e-18*(24.6/1000)**3   # H + He, E^-3 approx
print(f'photoabsorption at 1 keV: {sig1keV:.1e} cm^2 per H -> tau=1 at N_H={1/sig1keV:.1e}; ambient N_H={NH:.1e}; cloud column={nMC*2*RMC*pc:.1e}')
print(f'temperature fluctuation for heating x3: {3**(1/3.76)-1:.2f}')
v_sh = 0.4 * Rudr*pc / (trad*yr) / 1e5
print(f'shell velocity at tau_rad: {v_sh:.0f} km/s -> post-shock T = {3/16*0.6*mp*(v_sh*1e5)**2/kB:.1e} K')

print('--- Section 5.1: belt dust ---')
FEdd = 1.9e-2; LAGN = 0.026*LEdd
Lavg = FEdd*LEdd + 0.05*0.3*LEdd + LAGN
print(f'<L> = {Lavg:.1e} erg/s; light-crossing 2 r_b/c = {2*rb*pc/c/yr:.0f} yr')
Td = (Lavg/(16*np.pi*sigSB*(rb*pc)**2))**0.25
print(f'T_dust at r_b = {Td:.0f} K (peak ~{2900/Td*1e0:.0f} um); L_IR <~ fOm <L> = {fOm*Lavg:.1e} erg/s = {fOm*Lavg/LEdd:.3f} L_Edd')
print(f'dust mass in belt: {Mdisk*0.01:.1e} Zp Msun')

print('--- Section 5.2: variability ---')
print(f'F_Edd/F_0.1 = {1.9/7:.2f}; fractional change over dt=0.3-0.5 yr at t=2 yr: {5/3*0.3/2:.2f}-{5/3*0.5/2:.2f}')

print('--- Section 3.2: [OI] cooling ---')
g = np.array([5,3,1]); E = np.array([0,228,326.])
f = g*np.exp(-E/(1e3*T3)); f /= f.sum()
LOI = f[1]*8.9e-5*3.15e-14 + f[2]*1.7e-5*1.37e-14
print(f'[OI] LTE per O: {LOI:.1e}; /2.5 subthermal: {LOI/2.5:.1e}; per H: {LOI/2.5*4.9e-4:.1e} Zp vs H2 {LLTE/2:.1e} -> ratio {LOI/2.5*4.9e-4/(LLTE/2):.2f} Zp')

print('--- Section 5.5: H2 formation ---')
S = 1/(1+0.04*(900+20)**0.5+2e-3*900+8e-6*900**2)
print(f'HM79 sticking factor at 900 K: {S:.2f}')
tform = 1/(3e-17*S*nMC)/yr
print(f't_form = {tform:.1e} Zp^-1 yr; < tau_c for Zp > {tform/tc:.2f}; < t_gas (4e7 yr) for Zp > {tform/4e7:.1e}')
zeta_X = (LLTE/2)/(37*1.602e-12)   # ionizations per H per s, energy deposition = cooling
alpha = 2.6e-13*(900/1e4)**-0.7
xe = np.sqrt(zeta_X/(alpha*nMC))
kHm = 1.4e-18*900**0.928*np.exp(-900/16200)
print(f'X-ray ionization x_e ~ {xe:.1e}; H- channel t_form = {1/(kHm*xe*nMC)/yr:.1e} yr')

print('--- Section 3.3: zeta = 0.4 ---')
print(f'n_MC x{0.4**-0.44:.2f}, r_b x{0.4**0.15:.2f}, A_V x{0.4**-0.29:.2f}')
