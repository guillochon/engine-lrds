#!/usr/bin/env python3
"""
Metal-poor variant of the Guillochon & Loeb (2026) TDE-engine equilibrium.

C3 (collision cap) is metallicity-independent and is taken from GL2026's
losscone.py (bridged Cohn-Kulsrud flux).  The belt solve (virial, C1, C2)
and the compression criterion C4 are re-done with

  * H2 rovibrational cooling in LTE (Hollenbach & McKee 1979 fit),
    Lambda_vol = (n/2) * L_LTE(T),   L_LTE ~ 9.5e-22 (T/1000K)^3.76 erg/s/molecule
  * a cloud temperature T_MC (fiducial 200 K, the H2 floor) as a parameter
  * primordial composition (mu_H = 1.32, mu = 2.27)
  * a hot-phase cooling suppression zeta for the UDR radiative time,
    tau_rad prop. zeta^(-5/14)

Everything is monomial, so the belt is an exact power law in (m6, T_MC, zeta).
"""
import os, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np
import sympy as sp
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import losscone as L

# ---------------- constants ----------------
G   = sp.Float('6.67428e-8',30); msun=sp.Float('1.9889225e33',30)
yr  = sp.Float('3.1556926e7',30); pc=sp.Float('3.0856776e18',30)
kb  = sp.Float('1.3806503e-16',30); mp=sp.Float('1.67262158e-24',30)
cc  = sp.Float('2.99792458e10',30); sigT=sp.Float('6.6524e-25',30)
km  = sp.Float('1e5',30)
gam = sp.Rational(7,4); mstar = msun/2
epsacc = sp.Rational(1,20)

# symbols
m6, sig, Mach, rc, rb = sp.symbols('m6 sigma Mach R_MC r_b', positive=True)
T2, zeta = sp.symbols('T2 zeta', positive=True)   # T_MC/200K, hot-phase cooling factor
lm, ls, lM, lr, lb, lT, lz = sp.symbols('lm ls lM lr lb lT lz')

mh  = sp.Float('1e6',30)*m6*msun
Nst = 2*mh/mstar
ah  = 2*G*mh/sig**2
n0  = (3-gam)*Nst/(4*sp.pi*ah**3)
trelax_st = sp.Float('0.34',30)*sig**3/(G**2*mstar**2*n0*10)
E   = sp.Float('1.8e50',30)*m6**sp.Rational(1,3)
rhoMC = 9*mh/(4*sp.pi*rb**3)
Vmc = sp.Rational(4,3)*sp.pi*rc**3
Mmc = rhoMC*Vmc
LEdd = 4*sp.pi*G*mh*mp*cc/sigT

def _solve(eqs, unk, params):
    sub = {m6: sp.exp(lm), sig: sp.exp(ls), Mach: sp.exp(lM), rc: sp.exp(lr),
           rb: sp.exp(lb), T2: sp.exp(lT), zeta: sp.exp(lz)}
    n = len(unk); A = sp.zeros(n,n); B = sp.zeros(n,1)
    for i,e in enumerate(eqs):
        Lx = sp.expand_log(sp.log(sp.powsimp(e.subs(sub), force=True)), force=True)
        for j,v in enumerate(unk): A[i,j] = sp.diff(Lx, v)
        rest = Lx.subs({v:0 for v in unk})
        B[i] = -rest
    return A.solve(B)

def make_model(kind):
    """kind = 'CO' (GL2026) or 'H2' (metal-poor)."""
    d = {}
    if kind == 'CO':
        muH, mucl = sp.Float('1.4',30), sp.Float('2.33',30)
        cs = sp.sqrt(kb*10/(mucl*mp))
        nMC = rhoMC/(muH*mp)
        Lvol = sp.Float('1.3e-27',30)*nMC**2
        trad = sp.Float('1.7e4',30)*yr*(E/sp.Float('1e50',30))**sp.Rational(4,17)*nMC**sp.Rational(-9,17)
    else:
        muH, mucl = sp.Float('1.32',30), sp.Float('2.27',30)     # primordial X=0.76
        cs = sp.sqrt(kb*200*T2/(mucl*mp))
        nMC = rhoMC/(muH*mp)
        # HM79 LTE rotational cooling per H2 molecule, power-law form
        L_LTE = sp.Float('9.5e-22',30)*(200*T2/1000)**sp.Float('3.76',30)
        Lvol = (nMC/2)*L_LTE
        trad = sp.Float('1.7e4',30)*yr*(E/sp.Float('1e50',30))**sp.Rational(4,17) \
               *nMC**sp.Rational(-9,17)*zeta**sp.Rational(-5,14)
    Rudr = sp.Float('1.15',30)*(E/rhoMC)**sp.Rational(1,5)*trad**sp.Rational(2,5)
    fUDR = Rudr**2/(4*rb**2)
    sigcl = Mach*cs; tc = 2*rc/sigcl
    d.update(cs=cs, nMC=nMC, Lvol=Lvol, trad=trad, Rudr=Rudr, fUDR=fUDR,
             sigcl=sigcl, tc=tc, muH=muH)
    return d

# ---------- C3 from GL2026 code: sigma_eq(m6), Gamma_eq(m6) at flattening f ----------
def sig_coll(m6v, fst, LamC=3.0):
    mh_ = 1e6*m6v*L.msun
    rt = (mh_/L.MSTAR)**(1/3)*L.RSTAR
    hi = np.sqrt(2*L.G*mh_/(40*rt))
    f = lambda x: np.log(L.Cusp.engine(m6v, np.exp(x)).rate_coll(LamC)/fst/L.Cusp.engine(m6v, np.exp(x)).rate())
    return np.exp(brentq(f, np.log(2e6), np.log(hi)))

def fitpow(fst, xs=(0.03,0.3,3,30,300)):
    s = np.array([sig_coll(m, fst) for m in xs])
    g = np.array([L.Cusp.engine(m, sv).rate()*L.yr for m, sv in zip(xs, s)])
    a = np.log10(np.asarray(xs))
    ps = np.polyfit(a, np.log10(s/L.km), 1); pg = np.polyfit(a, np.log10(g), 1)
    return (10**ps[1], ps[0]), (10**pg[1], pg[0])

def fit_ck_sigma_m6():
    sigs = np.array([150,200,250,300,400,500,650])*L.km; m6s = np.array([0.3,1.0,3.0,10.0])
    rows=[]
    for mv in m6s:
        for sv in sigs:
            try: g = L.Cusp.engine(mv, sv).rate()*L.yr
            except Exception: continue
            if np.isfinite(g) and g>0: rows.append((np.log(mv), np.log(sv/(300*L.km)), np.log(g)))
    rows=np.array(rows); A=np.column_stack([np.ones(len(rows)), rows[:,1], rows[:,0]])
    c,_,_,_=np.linalg.lstsq(A, rows[:,2], rcond=None)
    return (sp.Float(str(np.exp(c[0])),30)/yr*(sig/(300*km))**sp.Float(str(c[1]),30)*m6**sp.Float(str(c[2]),30)), c

def powlaw_str(x, unit=1, params=(lm,)):
    """x is a log-expression in lm (and lT, lz). return prefactor & exponents."""
    x = sp.expand(x)
    pref = float(sp.exp(x.subs({p:0 for p in params}))/unit) if unit!=1 else float(sp.exp(x.subs({p:0 for p in params})))
    exps = {str(p): float(sp.diff(x, p)) for p in params}
    return pref, exps

def belt(kind, Gp, Gs, F_OM, F_ST):
    """Solve virial + C1 + C2 for (Mach, R_MC, r_b) given Gamma_eq = Gp m6^Gs."""
    md = make_model(kind)
    Geq = sp.Float(str(Gp),30)/yr*m6**sp.Float(str(Gs),30)
    beta = 1/F_OM; BEAM = 1
    Nmc = F_OM*4*rb**2/rc**2
    eqs = [Mach**2*md['cs']**2/(sp.Rational(4,5)*sp.pi*G*rhoMC*rc**2),
           Geq*md['fUDR']*md['tc']*beta,
           Geq*E*BEAM/(md['Lvol']*Nmc*Vmc)]
    s = _solve(eqs, [lM, lr, lb], None)
    logsub = {Mach: sp.exp(s[0]), rc: sp.exp(s[1]), rb: sp.exp(s[2]), m6: sp.exp(lm), T2: sp.exp(lT), zeta: sp.exp(lz)}
    def ev(expr, unit=1):
        x = sp.expand_log(sp.log(sp.powsimp((expr/unit).subs(logsub), force=True)), force=True)
        return powlaw_str(x, 1, (lm, lT, lz))
    out = {}
    out['Mach'] = ev(Mach); out['R_MC[pc]'] = ev(rc, pc); out['r_b[pc]'] = ev(rb, pc)
    out['n_MC[cm-3]'] = ev(md['nMC']); out['sigma_cl[km/s]'] = ev(md['sigcl'], km)
    out['tau_c[yr]'] = ev(md['tc'], yr); out['M_MC[Msun]'] = ev(Mmc, msun)
    out['N_MC'] = ev(Nmc); out['M_disk[Msun]'] = ev(Mmc*Nmc, msun)
    out['R_UDR[pc]'] = ev(md['Rudr'], pc); out['tau_rad[yr]'] = ev(md['trad'], yr)
    out['f_UDR'] = ev(md['fUDR']); out['L_UDR[erg/s]'] = ev(Geq*E*BEAM)
    vc2 = G*mh/rb; rhobg = rhoMC*md['sigcl']**2/vc2; nbg = rhobg/(md['muH']*mp)
    mdot = 4*sp.pi*rb**2*rhobg*sp.sqrt(vc2)*sp.Float('1e-3',30)
    out['n_bg[cm-3]'] = ev(nbg)
    out['L_AGN/L_Edd(fB=1e-3)'] = ev(epsacc*mdot*cc**2/LEdd)
    out['N_H,bg[cm-2]'] = ev(nbg*rb)
    out['A_V(solar dust)'] = ev(nbg*rb/sp.Float('2.2e21',30))
    out['Mdot_TDE/Mdot_amb'] = ev((mstar/2)*Geq/mdot)
    out['Jeans check sig^2 R/GM'] = ev(md['sigcl']**2*rc/(G*Mmc))
    out['hits per crossing'] = ev(Geq*md['tc'])
    out['H/R_MC'] = ev(F_ST*rb/rc)
    out['t_gas[yr]'] = ev(Mmc*Nmc/(mstar*Geq), yr)
    # turbulent dissipation heating per volume vs radiative cooling per volume
    heat = rhoMC*md['sigcl']**3/(2*rc)
    out['heat_turb/cool'] = ev(heat/md['Lvol'])
    out['Sigma_gas[Msun/pc2]'] = ev(Mmc*Nmc/(sp.pi*rb**2), msun/pc**2)
    out['R_UDR/R_MC'] = ev(md['Rudr']/rc)
    return out, md, s

def sigma_comp(kind, F_OM, G_ck):
    md = make_model(kind)
    beta = 1/F_OM; Nmc_ = F_OM*4*rb**2/rc**2
    trMC = sp.Float('0.34',30)*sig**3/(G**2*Mmc*rhoMC*10)
    eqs = [Mach**2*md['cs']**2/(sp.Rational(4,5)*sp.pi*G*rhoMC*rc**2),
           G_ck*md['fUDR']*md['tc']*beta,
           G_ck*E/(md['Lvol']*Nmc_*Vmc),
           trMC/trelax_st]
    x = sp.expand(_solve(eqs, [ls, lM, lr, lb], None)[0])
    x = x.subs({m6: sp.exp(lm)})
    return powlaw_str(x, km, (lm, lT, lz))

def show(out, title):
    print('\n=== %s ===' % title)
    for k,(p,e) in out.items():
        s = f"{k:28s} = {p:10.3e}"
        for name,val in e.items():
            if abs(val) > 1e-6:
                lab = {'lm':'m6','lT':'T2','lz':'zeta'}[name]
                s += f" {lab}^{val:+.3f}"
        print(s)

if __name__ == '__main__':
    F_ST = 0.22; F_OM = sp.Rational(22,100)
    (sp_, ss), (gp, gs) = fitpow(F_ST)
    print(f"C3 (GL2026, f_*={F_ST}): sigma_eq = {sp_:.1f} m6^{ss:.3f} km/s ; Gamma_eq = {gp:.3e} m6^{gs:.3f} /yr")
    ah_pc = 2*6.67428e-8*1e6*1.9889225e33/(sp_*1e5)**2/3.0856776e18
    rho0 = 0.5*1.25*(2e6/0.5)/(4*np.pi*ah_pc**3)
    print(f"   a_h = {ah_pc:.3f} pc, rho_0 = {rho0:.2e} Msun/pc^3 at m6=1")
    G_ck, cfit = fit_ck_sigma_m6()
    print(f"   CK fit: Gamma = {np.exp(cfit[0]):.3e} (sig/300)^{cfit[1]:.3f} m6^{cfit[2]:.3f}")

    outCO, _, _ = belt('CO', gp, gs, F_OM, F_ST); show(outCO, 'GL2026 reproduction: CO cooling, T=10 K, f_*=f_Omega=0.22')
    pc_, e_ = sigma_comp('CO', F_OM, G_ck)
    floor = (sp_/pc_)**(1/(e_['lm']-ss))*1e6
    print(f"sigma_comp(CO) = {pc_:.1f} m6^{e_['lm']:.3f} km/s  -> mass floor {floor:.2e} Msun")

    outH2, _, _ = belt('H2', gp, gs, F_OM, F_ST); show(outH2, 'Metal-poor: H2 LTE cooling, T=200 T2 K, zeta')
    pc_, e_ = sigma_comp('H2', F_OM, G_ck)
    print(f"sigma_comp(H2) = {pc_:.1f} m6^{e_['lm']:.3f} T2^{e_['lT']:.3f} zeta^{e_['lz']:.3f} km/s")
    floor = (sp_/pc_)**(1/(e_['lm']-ss))*1e6
    print(f"   -> mass floor {floor:.2e} Msun at T2=1, zeta=1;  floor prop. T2^{-e_['lT']/(e_['lm']-ss):.3f} zeta^{-e_['lz']/(e_['lm']-ss):.3f}")

# ======================================================================
#  Self-consistent temperature: 4 unknowns (Mach, R_MC, r_b, T2),
#  4 equations (virial, C1, C2, local thermal balance heat_turb = cool)
# ======================================================================
def belt_T(Gp, Gs, F_OM, F_ST, Gexpr=None):
    md = make_model('H2')
    Geq = sp.Float(str(Gp),30)/yr*m6**sp.Float(str(Gs),30) if Gexpr is None else Gexpr
    beta = 1/F_OM; BEAM = 1
    Nmc = F_OM*4*rb**2/rc**2
    heat = rhoMC*md['sigcl']**3/(2*rc)
    eqs = [Mach**2*md['cs']**2/(sp.Rational(4,5)*sp.pi*G*rhoMC*rc**2),
           Geq*md['fUDR']*md['tc']*beta,
           Geq*E*BEAM/(md['Lvol']*Nmc*Vmc),
           heat/md['Lvol']]
    s = _solve(eqs, [lM, lr, lb, lT], None)
    logsub = {Mach: sp.exp(s[0]), rc: sp.exp(s[1]), rb: sp.exp(s[2]), T2: sp.exp(s[3]),
              m6: sp.exp(lm), zeta: sp.exp(lz), sig: sp.exp(ls)}
    def ev(expr, unit=1):
        x = sp.expand_log(sp.log(sp.powsimp((expr/unit).subs(logsub), force=True)), force=True)
        return powlaw_str(x, 1, (lm, lz, ls))
    out = {}
    out['T_MC[K]'] = ev(200*T2)
    out['c_s[km/s]'] = ev(md['cs'], km)
    out['Mach'] = ev(Mach); out['R_MC[pc]'] = ev(rc, pc); out['r_b[pc]'] = ev(rb, pc)
    out['n_MC[cm-3]'] = ev(md['nMC']); out['sigma_cl[km/s]'] = ev(md['sigcl'], km)
    out['tau_c[yr]'] = ev(md['tc'], yr); out['M_MC[Msun]'] = ev(Mmc, msun)
    out['N_MC'] = ev(Nmc); out['M_disk[Msun]'] = ev(Mmc*Nmc, msun)
    out['R_UDR[pc]'] = ev(md['Rudr'], pc); out['tau_rad[yr]'] = ev(md['trad'], yr)
    out['f_UDR'] = ev(md['fUDR']); out['L_UDR=L_H2[erg/s]'] = ev(Geq*E*BEAM)
    vc2 = G*mh/rb; rhobg = rhoMC*md['sigcl']**2/vc2; nbg = rhobg/(md['muH']*mp)
    mdot = 4*sp.pi*rb**2*rhobg*sp.sqrt(vc2)*sp.Float('1e-3',30)
    out['n_bg[cm-3]'] = ev(nbg)
    out['L_AGN/L_Edd(fB=1e-3)'] = ev(epsacc*mdot*cc**2/LEdd)
    out['N_H,bg[cm-2]'] = ev(nbg*rb)
    out['A_V(Z=Zsun dust)'] = ev(nbg*rb/sp.Float('2.2e21',30))
    out['Mdot_TDE/Mdot_amb'] = ev((mstar/2)*Geq/mdot)
    out['hits per crossing'] = ev(Geq*md['tc'])
    out['H/R_MC'] = ev(F_ST*rb/rc)
    out['t_gas[yr]'] = ev(Mmc*Nmc/(mstar*Geq), yr)
    out['Sigma_gas[Msun/pc2]'] = ev(Mmc*Nmc/(sp.pi*rb**2), msun/pc**2)
    out['R_UDR/R_MC'] = ev(md['Rudr']/rc)
    out['n_MC/n_cr(1e4)'] = ev(md['nMC']/sp.Float('1e4',30))
    out['L_H2 per gas mass[erg/s/g]'] = ev(Geq*E*BEAM/(Mmc*Nmc))
    out['t_cool = 3/2 kT n/Lvol [yr]'] = ev(sp.Rational(3,2)*kb*200*T2*(md['nMC']/2)/md['Lvol'], yr)
    return out, md, s, logsub

def sigma_comp_T(F_OM, G_ck):
    """C4 with T solved: 5 unknowns (sigma, Mach, R_MC, r_b, T2)."""
    md = make_model('H2')
    beta = 1/F_OM; Nmc_ = F_OM*4*rb**2/rc**2
    trMC = sp.Float('0.34',30)*sig**3/(G**2*Mmc*rhoMC*10)
    heat = rhoMC*md['sigcl']**3/(2*rc)
    eqs = [Mach**2*md['cs']**2/(sp.Rational(4,5)*sp.pi*G*rhoMC*rc**2),
           G_ck*md['fUDR']*md['tc']*beta,
           G_ck*E/(md['Lvol']*Nmc_*Vmc),
           heat/md['Lvol'],
           trMC/trelax_st]
    x = sp.expand(_solve(eqs, [ls, lM, lr, lb, lT], None)[0]).subs({m6: sp.exp(lm)})
    return powlaw_str(x, km, (lm, lz))

if __name__ == '__main__':
    print('\n\n############ SELF-CONSISTENT TEMPERATURE ############')
    outT, mdT, sT, lsub = belt_T(gp, gs, F_OM, F_ST); show(outT, 'Metal-poor H2 engine, T solved from turbulent heating = H2 cooling, f_*=0.22')
    pcT, eT = sigma_comp_T(F_OM, G_ck)
    floorT = (sp_/pcT)**(1/(eT['lm']-ss))*1e6
    print(f"sigma_comp(H2,T solved) = {pcT:.1f} m6^{eT['lm']:.3f} zeta^{eT['lz']:.3f} km/s -> mass floor {floorT:.2e} Msun (zeta=1); floor prop. zeta^{-eT['lz']/(eT['lm']-ss):.3f}")
    for z in [1.0, 0.1, 0.03]:
        fl = (sp_/(pcT*z**eT['lz']))**(1/(eT['lm']-ss))*1e6
        print(f"   zeta={z}: floor = {fl:.2e} Msun")
