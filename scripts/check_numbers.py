"""Evaluate the power laws in scripts/lrd_results.json and scripts/lrd_scan.json at
the points quoted in paper/lrd_engine.tex, and print every derived number the text
uses (Table 1, the temperature window, the mass floor, the (f_*, T) wedge, A_V, the
duty cycle and host abundance), so the manuscript can be checked line by line.

    python scripts/check_numbers.py
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, 'lrd_results.json')))
S = json.load(open(os.path.join(HERE, 'lrd_scan.json')))
b = R['belt']


def ev(key, m6=1., T3=1., x=1., z=1.):
    p, e = b[key]
    return p * m6 ** e['m6'] * T3 ** e['T3'] * x ** e['x'] * z ** e['zeta']


def law(key, unit=''):
    p, e = b[key]
    s = '%.3g' % p
    for n, v in e.items():
        if abs(v) > 1e-6:
            s += ' %s^%+.2f' % (n, v)
    return s + (' ' + unit if unit else '')


print('== collision cap (Sec. 2, eq. Geq) ==')
c = R['C3']
print('sigma_eq = %.0f m6^%.2f km/s;  Gamma_eq = %.2e m6^%.2f /yr' % (c['sig'], c['ss'], c['gp'], c['gs']))
G = 6.67428e-8; msun = 1.9889225e33; pc = 3.0856776e18; yr = 3.1556926e7
ah = 2 * G * 1e6 * msun / (c['sig'] * 1e5) ** 2 / pc
print('a_h = %.3f pc  (exponent %.2f)' % (ah, 1 - 2 * c['ss']))

print('\n== belt power laws (Sec. 4.1) ==')
for k in ['r_b[pc]', 'n_MC[cm-3]', 'tau_c[yr]', 'R_UDR[pc]', 'tau_rad[yr]', 'R_MC[pc]', 'Mach', 'sigma_cl[km/s]',
          'M_MC[Msun]', 'N_MC', 'M_disk[Msun]', 'L_H2[erg/s]', 'n_bg[cm-3]', 'L_AGN/L_Edd', 'N_H[cm-2]',
          'Mdot_TDE/Mdot_amb', 'H/R_MC', 't_gas[yr]', 't_cool[yr]', 'tau_es', 'A_V(Zsun)', 'heat_turb/cool']:
    print('%-20s = %s' % (k, law(k)))
v_udr = 0.4 * ev('R_UDR[pc]') * pc / (ev('tau_rad[yr]') * yr) / 1e5
print('UDR shell speed at the radiative transition, (2/5) R_UDR/tau_rad = %.0f km/s' % v_udr)

print('\n== temperature window (Sec. 4.2) ==')
eH = b['H/R_MC'][1]; eM = b['Mach'][1]
Tmin, Tmax = R['Tmin'] * 1e3, R['Tmax'] * 1e3
print('T_min = %.0f K  m6^%.2f x^%.2f   (H = R_MC)' % (Tmin, -eH['m6'] / eH['T3'], -eH['x'] / eH['T3']))
print('T_max = %.0f K  m6^%.2f x^%.2f   (Mach = 1)' % (Tmax, -eM['m6'] / eM['T3'], -eM['x'] / eM['T3']))
print('T_max/T_min = %.2f m6^%.3f;  geometric center %.0f K m6^%.2f' % (
    Tmax / Tmin, -eM['m6'] / eM['T3'] + eH['m6'] / eH['T3'], np.sqrt(Tmin * Tmax),
    0.5 * (-eH['m6'] / eH['T3'] - eM['m6'] / eM['T3'])))
T200 = 0.2
print('at T = 200 K: R_MC = %.0f pc, Mach = %.0f, M_disk = %.1e Msun' % (
    ev('R_MC[pc]', T3=T200), ev('Mach', T3=T200), ev('M_disk[Msun]', T3=T200)))
p, e = b['heat_turb/cool']
T3tb = p ** (-1 / e['T3'])
print('turbulent thermal balance: T = %.0f K m6^%.2f, R_MC = %.1f pc, M_disk = %.1e Msun' % (
    T3tb * 1e3, -e['m6'] / e['T3'], ev('R_MC[pc]', T3=T3tb), ev('M_disk[Msun]', T3=T3tb)))
for x in [1e-2]:
    print('x_H2 = %g: window %.0f-%.0f K' % (x, Tmin * x ** (-eH['x'] / eH['T3']), Tmax * x ** (-eM['x'] / eM['T3'])))

print('\n== Table 1, metal-poor column at T_MC = 900 K, m6 = 1 ==')
T3c = 0.9
for k in ['R_MC[pc]', 'Mach', 'sigma_cl[km/s]', 'M_MC[Msun]', 'N_MC', 'M_disk[Msun]', 'n_bg[cm-3]', 'N_H[cm-2]',
          'A_V(Zsun)', 'L_AGN/L_Edd', 'Mdot_TDE/Mdot_amb', 't_gas[yr]', 'tau_es']:
    print('%-20s = %.3g' % (k, ev(k, T3=T3c)))
print('M_disk across the window: %.2e - %.2e Msun' % (ev('M_disk[Msun]', T3=Tmax / 1e3), ev('M_disk[Msun]', T3=Tmin / 1e3)))
print('CO cooling per H nucleus 1.3e-27 n_MC = %.2e erg/s;  H2 at 1e3 K L_LTE/2 = %.2e erg/s' % (
    1.3e-27 * ev('n_MC[cm-3]'), 9.5e-22 / 2))

print('\n== compression floor (Sec. 4.3) ==')
pcmp, ecmp = R['comp']
slope = ecmp['m6'] - c['ss']
fl = (c['sig'] / pcmp) ** (1 / slope)
print('sigma_comp = %.0f m6^%.2f T3^%.2f x^%.2f zeta^%.2f km/s' % (pcmp, ecmp['m6'], ecmp['T3'], ecmp['x'], ecmp['zeta']))
print('floor = %.2e T3^%.2f x^%.2f zeta^%.2f Msun; at T3=0.8: %.2e' % (
    fl * 1e6, -ecmp['T3'] / slope, -ecmp['x'] / slope, -ecmp['zeta'] / slope, fl * 1e6 * 0.8 ** (-ecmp['T3'] / slope)))
m6f = 0.7
Tminf = Tmin * m6f ** (-eH['m6'] / eH['T3'])
print('at m6 = 0.7: T_min = %.0f K, floor(T_min) = %.2e Msun' % (Tminf, fl * 1e6 * (Tminf / 1e3) ** (-ecmp['T3'] / slope)))

print('\n== AGN radiation pressure (Sec. 4.4) ==')
fEdd = ev('L_AGN/L_Edd', T3=T3c); LEdd = 1.26e38 * 1e6; tau_es = ev('tau_es', T3=T3c)
pTDE = 1.7e33 * c['gp'] / 1e-2
for Zp in [1.0, 0.1]:
    tau = max(30 * fEdd * Zp, tau_es)
    print('Zp = %.1f: tau = %.3f, p_TDE/p_AGN = %.0f' % (Zp, tau, pTDE / (tau * fEdd * LEdd / 2.99792458e10)))
print('tau_es = %.3f' % tau_es)

print('\n== (f_*, T) wedge (Sec. 4.4, Fig. 2) ==')
for r in S:
    tlo = (1 / r['HR'][0]) ** (1 / r['HR'][1]['T3']) * 1e3
    thi = r['Mach'][0] ** (-1 / r['Mach'][1]['T3']) * 1e3
    print('f_* = %.2f: sigma = %.0f, Gamma = %.2e, rho0 = %.2e, window %.0f-%.0f K, floor %.2e' % (
        r['f'], r['sig'], r['G'], r['rho0'], tlo, thi, r['floor'] * 1e6))
fs = np.array([r['f'] for r in S]); rho0 = np.array([r['rho0'] for r in S])
print('f_* at rho0 = 1e8: %.2f;  at 1e9: %.2f;  d log rho0 / d log f = %.2f' % (
    np.interp(8, np.log10(rho0), fs), np.interp(9, np.log10(rho0), fs),
    np.polyfit(np.log10(fs), np.log10(rho0), 1)[0]))

print('\n== mass range (Sec. 4.5) ==')
for m6 in [0.7, 1, 2, 3, 5]:
    tl = Tmin * m6 ** (-eH['m6'] / eH['T3']); th = Tmax * m6 ** (-eM['m6'] / eM['T3'])
    print('m6 = %g: window %.0f-%.0f K, n_MC = %.1e, Gamma = %.1e, Mach(center) = %.1f' % (
        m6, tl, th, ev('n_MC[cm-3]', m6=m6), c['gp'] * m6 ** c['gs'], ev('Mach', m6=m6, T3=np.sqrt(tl * th) / 1e3)))

print('\n== LRD consequences (Sec. 5) ==')
print('A_V = %.3g Zp at T3=1;  %.3g Zp at T3=0.9' % (ev('A_V(Zsun)'), ev('A_V(Zsun)', T3=T3c)))
print('L_AGN = %.1e erg/s (f_B=1e-3, T3=0.9)' % (fEdd * LEdd))
FEdd = 1.3e-2 * c['gp'] / 1e-2
print('GL26 eq. 13: F_Edd = 1.3e-2 m6^-0.4 Gamma_-2 = %.2e m6^%.2f;  F_0.1 = 10^0.6 F_Edd = %.2e' % (FEdd, -0.4 + c['gs'], FEdd * 10 ** 0.6))
print('n_host = 1e-5 / F_0.1 = %.1e cMpc^-3' % (1e-5 / (FEdd * 10 ** 0.6)))
tgas = ev('t_gas[yr]', T3=T3c)
print('t_gas = %.1e yr; stellar mass consumed 0.25 Msun x Gamma x t_gas = %.1e Msun' % (tgas, 0.25 * c['gp'] * tgas))
print('t_cusp = N_*/Gamma = %.1e yr m6^%.2f' % (4e6 / c['gp'], 1 - c['gs']))
n_form = ev('n_MC[cm-3]')
S900 = 1 / (1 + 0.04 * (900 + 20) ** 0.5 + 2e-3 * 900 + 8e-6 * 900 ** 2)   # HM79 sticking factor at T_gas = 900 K
t_form = 1 / (3e-17 * S900 * n_form) / yr
print('grain H2 formation time 1/(3e-17 S n_MC), S(900 K) = %.2f: %.1e yr / Zp; < tau_c for Zp > %.2f; < t_gas for Zp > %.1e' % (
    S900, t_form, t_form / ev('tau_c[yr]'), t_form / tgas))
