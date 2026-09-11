#!/usr/bin/env python3
"""Metal-poor TDE engine: full solve, windows, figures. Uses GL2026 losscone.py for C3."""
import os, sys, warnings, json
warnings.filterwarnings('ignore')
import numpy as np, sympy as sp
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import losscone as L
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family']='serif'; rcParams['font.serif']=['DejaVu Serif']; rcParams['mathtext.fontset']='stix'
rcParams['xtick.direction']='in'; rcParams['ytick.direction']='in'; rcParams['xtick.top']=True; rcParams['ytick.right']=True
rcParams['pdf.fonttype']=42

G=sp.Float('6.67428e-8',30); msun=sp.Float('1.9889225e33',30); yr=sp.Float('3.1556926e7',30)
pc=sp.Float('3.0856776e18',30); kb=sp.Float('1.3806503e-16',30); mp=sp.Float('1.67262158e-24',30)
cc=sp.Float('2.99792458e10',30); sigT=sp.Float('6.6524e-25',30); km=sp.Float('1e5',30)
gam=sp.Rational(7,4); mstar=msun/2; epsacc=sp.Rational(1,20)

m6,sig,Mach,rc,rb,T3,xh2,zeta = sp.symbols('m6 sigma Mach R_MC r_b T3 x zeta', positive=True)
lm,ls,lM,lr,lb,lT,lx,lz = sp.symbols('lm ls lM lr lb lT lx lz')
PAR=[lm,lT,lx,lz]; PNAME={lm:'m6',lT:'T3',lx:'x',lz:'zeta',ls:'sig'}

mh=sp.Float('1e6',30)*m6*msun; Nst=2*mh/mstar; ah=2*G*mh/sig**2
n0=(3-gam)*Nst/(4*sp.pi*ah**3)
trelax_st=sp.Float('0.34',30)*sig**3/(G**2*mstar**2*n0*10)
E=sp.Float('1.8e50',30)*m6**sp.Rational(1,3)
rhoMC=9*mh/(4*sp.pi*rb**3); Vmc=sp.Rational(4,3)*sp.pi*rc**3; Mmc=rhoMC*Vmc
LEdd=4*sp.pi*G*mh*mp*cc/sigT
# primordial composition
muH=sp.Float('1.32',30); mucl=sp.Float('2.27',30)
cs=sp.sqrt(kb*1000*T3/(mucl*mp))
nMC=rhoMC/(muH*mp)
L_LTE=sp.Float('9.5e-22',30)*T3**sp.Float('3.76',30)          # HM79 LTE rotational, per H2
Lvol=xh2*(nMC/2)*L_LTE
trad=sp.Float('1.7e4',30)*yr*(E/sp.Float('1e50',30))**sp.Rational(4,17)*nMC**sp.Rational(-9,17)*zeta**sp.Rational(-5,14)
Rudr=sp.Float('1.15',30)*(E/rhoMC)**sp.Rational(1,5)*trad**sp.Rational(2,5)
fUDR=Rudr**2/(4*rb**2); sigcl=Mach*cs; tc=2*rc/sigcl
SUB={m6:sp.exp(lm),sig:sp.exp(ls),Mach:sp.exp(lM),rc:sp.exp(lr),rb:sp.exp(lb),T3:sp.exp(lT),xh2:sp.exp(lx),zeta:sp.exp(lz)}

def _solve(eqs, unk):
    n=len(unk); A=sp.zeros(n,n); B=sp.zeros(n,1)
    for i,e in enumerate(eqs):
        Lx=sp.expand_log(sp.log(sp.powsimp(e.subs(SUB),force=True)),force=True)
        for j,v in enumerate(unk): A[i,j]=sp.diff(Lx,v)
        B[i]=-Lx.subs({v:0 for v in unk})
    return A.solve(B)

def pl(x, params):
    x=sp.expand(x); pref=float(sp.exp(x.subs({p:0 for p in params})))
    return pref,{PNAME[p]:float(sp.diff(x,p)) for p in params}

# ---------- C3 via GL2026 code ----------
def sig_coll(m6v,fst,LamC=3.0):
    mh_=1e6*m6v*L.msun; rt=(mh_/L.MSTAR)**(1/3)*L.RSTAR; hi=np.sqrt(2*L.G*mh_/(40*rt))
    f=lambda x: np.log(L.Cusp.engine(m6v,np.exp(x)).rate_coll(LamC)/fst/L.Cusp.engine(m6v,np.exp(x)).rate())
    return np.exp(brentq(f,np.log(2e6),np.log(hi)))
def fitpow(fst,xs=(0.03,0.3,3,30,300)):
    s=np.array([sig_coll(m,fst) for m in xs]); g=np.array([L.Cusp.engine(m,sv).rate()*L.yr for m,sv in zip(xs,s)])
    a=np.log10(np.asarray(xs)); ps=np.polyfit(a,np.log10(s/L.km),1); pg=np.polyfit(a,np.log10(g),1)
    return (10**ps[1],ps[0]),(10**pg[1],pg[0])
def fit_ck():
    sigs=np.array([150,200,250,300,400,500,650])*L.km; rows=[]
    for mv in [0.3,1.0,3.0,10.0]:
        for sv in sigs:
            try: g=L.Cusp.engine(mv,sv).rate()*L.yr
            except Exception: continue
            if np.isfinite(g) and g>0: rows.append((np.log(mv),np.log(sv/(300*L.km)),np.log(g)))
    rows=np.array(rows); A=np.column_stack([np.ones(len(rows)),rows[:,1],rows[:,0]])
    c=np.linalg.lstsq(A,rows[:,2],rcond=None)[0]
    return sp.Float(str(np.exp(c[0])),30)/yr*(sig/(300*km))**sp.Float(str(c[1]),30)*m6**sp.Float(str(c[2]),30)

def belt(Gp,Gs,F_OM,F_ST,Zp_sym=None):
    Geq=sp.Float(str(Gp),30)/yr*m6**sp.Float(str(Gs),30); beta=1/F_OM
    Nmc=F_OM*4*rb**2/rc**2
    eqs=[Mach**2*cs**2/(sp.Rational(4,5)*sp.pi*G*rhoMC*rc**2), Geq*fUDR*tc*beta, Geq*E/(Lvol*Nmc*Vmc)]
    s=_solve(eqs,[lM,lr,lb])
    sub={Mach:sp.exp(s[0]),rc:sp.exp(s[1]),rb:sp.exp(s[2])}; sub.update({k:v for k,v in SUB.items() if k in (m6,T3,xh2,zeta)})
    ev=lambda expr,unit=1: pl(sp.expand_log(sp.log(sp.powsimp((expr/unit).subs(sub),force=True)),force=True),PAR)
    vc2=G*mh/rb; rhobg=rhoMC*sigcl**2/vc2; nbg=rhobg/(muH*mp)
    mdot=4*sp.pi*rb**2*rhobg*sp.sqrt(vc2)*sp.Float('1e-3',30); fEdd=epsacc*mdot*cc**2/LEdd
    out={'Mach':ev(Mach),'R_MC[pc]':ev(rc,pc),'r_b[pc]':ev(rb,pc),'n_MC[cm-3]':ev(nMC),'c_s[km/s]':ev(cs,km),
         'sigma_cl[km/s]':ev(sigcl,km),'tau_c[yr]':ev(tc,yr),'M_MC[Msun]':ev(Mmc,msun),'N_MC':ev(Nmc),
         'M_disk[Msun]':ev(Mmc*Nmc,msun),'R_UDR[pc]':ev(Rudr,pc),'tau_rad[yr]':ev(trad,yr),'f_UDR':ev(fUDR),
         'L_H2[erg/s]':ev(Geq*E),'n_bg[cm-3]':ev(nbg),'L_AGN/L_Edd':ev(fEdd),'N_H[cm-2]':ev(nbg*rb),
         'A_V(Zsun)':ev(nbg*rb/sp.Float('2.2e21',30)),'Mdot_TDE/Mdot_amb':ev((mstar/2)*Geq/mdot),
         'H/R_MC':ev(F_ST*rb/rc),'t_gas[yr]':ev(Mmc*Nmc/(mstar*Geq),yr),'Sigma_gas[Msun/pc2]':ev(Mmc*Nmc/(sp.pi*rb**2),msun/pc**2),
         'heat_turb/cool':ev(rhoMC*sigcl**3/(2*rc)/Lvol),'v_c(r_b)[km/s]':ev(sp.sqrt(vc2),km),
         't_cool[yr]':ev(sp.Rational(3,2)*kb*1000*T3*nMC/(mucl/muH)/Lvol,yr), 'tau_es':ev(nbg*rb*sigT),
         'M_MC/M_h':ev(Mmc/mh), 'q_tidal Mmc/(rhoMC*rb^3)':ev(Mmc*4*sp.pi/(3*Mmc)*0+Mmc/(mh)*(rb/rc)**3)}
    return out,s

def sigma_comp(F_OM,G_ck):
    beta=1/F_OM; Nmc_=F_OM*4*rb**2/rc**2; trMC=sp.Float('0.34',30)*sig**3/(G**2*Mmc*rhoMC*10)
    eqs=[Mach**2*cs**2/(sp.Rational(4,5)*sp.pi*G*rhoMC*rc**2), G_ck*fUDR*tc*beta, G_ck*E/(Lvol*Nmc_*Vmc), trMC/trelax_st]
    x=sp.expand(_solve(eqs,[ls,lM,lr,lb])[0])
    p,e=pl(x,PAR); return p/1e5,e

def show(out,title):
    print('\n=== %s ==='%title)
    for k,(p,e) in out.items():
        s=f"{k:24s} = {p:10.3e}"
        for n,v in e.items():
            if abs(v)>1e-6: s+=f" {n}^{v:+.3f}"
        print(s)

if __name__=='__main__':
    F_ST=0.22; F_OM=sp.Rational(22,100)
    (sp_,ss),(gp,gs)=fitpow(F_ST); G_ck=fit_ck()
    print(f"C3: sigma={sp_:.1f} m6^{ss:.3f}; Gamma={gp:.3e} m6^{gs:.3f}")
    out,s=belt(gp,gs,F_OM,F_ST); show(out,'metal-poor belt, f_*=f_Om=0.22, T3 x zeta family')
    pcomp,ecomp=sigma_comp(F_OM,G_ck)
    print(f"sigma_comp = {pcomp:.1f} m6^{ecomp['m6']:.3f} T3^{ecomp['T3']:.3f} x^{ecomp['x']:.3f} zeta^{ecomp['zeta']:.3f}")
    slope=ecomp['m6']-ss
    fl=(sp_/pcomp)**(1/slope)
    print(f"mass floor m6 = {fl:.3f} T3^{-ecomp['T3']/slope:.3f} x^{-ecomp['x']/slope:.3f} zeta^{-ecomp['zeta']/slope:.3f}")
    # windows at m6=1, x=1, zeta=1
    HR=out['H/R_MC']; Ma=out['Mach']
    Tmin=(1/HR[0])**(1/HR[1]['T3']); Tmax=(Ma[0])**(-1/Ma[1]['T3'])
    print(f"T window (m6=1): H=R_MC at T3={Tmin:.3f}; Mach=1 at T3={Tmax:.3f}; ratio {Tmin/Tmax:.3f}")
    print(f"   T_min prop. m6^{-HR[1]['m6']/HR[1]['T3']:.3f}, T_max prop. m6^{-Ma[1]['m6']/Ma[1]['T3']:.3f}")
    print(f"   T_min prop. x^{-HR[1]['x']/HR[1]['T3']:.3f}, T_max prop. x^{-Ma[1]['x']/Ma[1]['T3']:.3f}")
    res={'C3':dict(sig=sp_,ss=ss,gp=gp,gs=gs),'belt':{k:[p,e] for k,(p,e) in out.items()},'comp':[pcomp,ecomp],'Tmin':Tmin,'Tmax':Tmax}
    json.dump(res,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'lrd_results.json'),'w'),indent=1)

    # ---------- scan over f_* : constraints in (f_*, T3) plane at m6=1 ----------
    fs=np.round(np.arange(0.10,0.601,0.05),2)
    rows=[]
    for f in fs:
        (sf,ssf),(gf,gsf)=fitpow(f); fr=sp.Rational(int(round(f*100)),100)
        o,_=belt(gf,gsf,fr,f)
        ah_pc=2*6.67428e-8*1e6*1.9889225e33/(sf*1e5)**2/3.0856776e18
        rho0=0.5*1.25*(2e6/0.5)/(4*np.pi*ah_pc**3)
        pc_,ec_=sigma_comp(fr,G_ck); fl_=(sf/pc_)**(1/(ec_['m6']-ssf))
        rows.append(dict(f=f,sig=sf,G=gf,Gs=gsf,rho0=rho0,HR=o['H/R_MC'],Mach=o['Mach'],fEdd=o['L_AGN/L_Edd'],
                         nH=o['N_H[cm-2]'],tau_es=o['tau_es'],floor=fl_,floorT=-ec_['T3']/(ec_['m6']-ssf),Mdisk=o['M_disk[Msun]']))
        print(f"f={f:.2f} sig={sf:.0f} G={gf:.2e} rho0={rho0:.2e} HR0={o['H/R_MC'][0]:.3e} Mach0={o['Mach'][0]:.3e} fEdd0={o['L_AGN/L_Edd'][0]:.3e} floor={fl_:.2e} Mdisk0={o['M_disk[Msun]'][0]:.2e}")
    json.dump(rows,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'lrd_scan.json'),'w'),indent=1,default=float)
