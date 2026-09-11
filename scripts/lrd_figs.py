"""Draw paper/twindow.pdf and paper/fplane.pdf from scripts/lrd_results.json and scripts/lrd_scan.json
(written by lrd_final.py).  python scripts/lrd_figs.py"""
import os, json, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family']='serif'; rcParams['mathtext.fontset']='stix'
rcParams['xtick.direction']='in'; rcParams['ytick.direction']='in'; rcParams['xtick.top']=True; rcParams['ytick.right']=True
rcParams['pdf.fonttype']=42; rcParams['axes.linewidth']=0.9
HERE=os.path.dirname(os.path.abspath(__file__)); PAPER=os.path.join(os.path.dirname(HERE),'paper')
R=json.load(open(os.path.join(HERE,'lrd_results.json'))); S=json.load(open(os.path.join(HERE,'lrd_scan.json')))
LEG=dict(frameon=True,facecolor='white',framealpha=0.8,edgecolor='none',labelspacing=0.32)
b=R['belt']
def ev(key, m6=1., T3=1., x=1., z=1.):
    p,e=b[key]; return p*m6**e['m6']*T3**e['T3']*x**e['x']*z**e['zeta']
BLUE=(0.15,0.40,0.80); RED=(0.80,0.15,0.15); GOLD=(0.85,0.62,0.10); GREEN=(0.55,0.82,0.55); PURP=(0.5,0.3,0.7)

# ---------------- Figure 1: belt vs T at m6=1, with GL2026 CO values ----------------
T=np.logspace(np.log10(300),np.log10(3000),400); T3=T/1e3
fig,ax=plt.subplots(1,2,figsize=(9.6,3.7),gridspec_kw={'wspace':0.32})
Tmin,Tmax=R['Tmin']*1e3,R['Tmax']*1e3
for a in ax: a.axvspan(Tmin,Tmax,color=GREEN,alpha=0.45,lw=0,zorder=0)
ax[0].plot(T,ev('R_MC[pc]',T3=T3),color=BLUE,lw=1.8,label=r'$R_{\rm MC}$ (pc)')
ax[0].plot(T,ev('r_b[pc]',T3=T3)*0.22,color=BLUE,lw=1.2,ls='--',label=r'$H=f_\Omega r_{\rm b}$ (pc)')
ax[0].plot(T,ev('Mach',T3=T3),color=RED,lw=1.8,label=r'$\mathcal{M}$')
ax[0].axhline(1.0,color=RED,lw=1.0,ls=':')
ax[0].plot(T,ev('M_disk[Msun]',T3=T3)/1e5,color=GOLD,lw=1.8,label=r'$M_{\rm disk}/10^{5}\,M_\odot$')
# GL2026 CO-cooled values (T=10 K) as markers at the left edge for reference
ax[0].scatter([310],[0.71],marker='o',color=BLUE,zorder=5,s=28)
ax[0].scatter([310],[3.9],marker='o',color=GOLD,zorder=5,s=28)
ax[0].text(325,0.71,r'GL26 $R_{\rm MC}$',fontsize=7.5,color=BLUE,va='center'); ax[0].text(325,3.9,r'GL26 $M_{\rm disk}$',fontsize=7.5,color=GOLD,va='center')
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xlim(300,3000); ax[0].set_ylim(0.05,60)
ax[0].set_xlabel(r'$T_{\rm MC}$ (K)',fontsize=12); ax[0].set_ylabel('belt properties',fontsize=12)
ax[0].legend(fontsize=8.5,loc='lower left',**LEG)
ax[0].text(np.sqrt(Tmin*Tmax),30,'allowed',color=(0.15,0.5,0.15),fontsize=10,fontweight='bold',ha='center')
ax[0].set_title(r'$M_{\rm h}=10^{6}\,M_\odot$, $f_\ast=f_\Omega=0.22$',fontsize=11,pad=6)
# right: constraint ratios
ax[1].plot(T,ev('H/R_MC',T3=T3),color=BLUE,lw=1.8,label=r'$H/R_{\rm MC}$')
ax[1].plot(T,ev('Mach',T3=T3),color=RED,lw=1.8,label=r'$\mathcal{M}$')
ax[1].plot(T,ev('n_MC[cm-3]',T3=T3)/1e4,color=PURP,lw=1.4,ls='-.',label=r'$n_{\rm MC}/n_{\rm cr}$, $n_{\rm cr}=10^{4}$ cm$^{-3}$')
ax[1].axhline(1.0,color='k',lw=1.0)
ax[1].set_xscale('log'); ax[1].set_yscale('log'); ax[1].set_xlim(300,3000); ax[1].set_ylim(0.05,60)
ax[1].set_xlabel(r'$T_{\rm MC}$ (K)',fontsize=12); ax[1].set_ylabel('constraint ratio (need $>1$)',fontsize=12)
ax[1].legend(fontsize=8.5,loc='lower left',**LEG)
ax[1].set_title(r'window $%.0f\,{\rm K} \leq T_{\rm MC} \leq %.0f$ K'%(Tmin,Tmax),fontsize=11,pad=6)
fig.tight_layout(pad=0.4); fig.savefig(os.path.join(PAPER,'twindow.pdf')); plt.close(fig)

# ---------------- Figure 2: (f_*, T) plane at m6=1 ----------------
fs=np.array([r['f'] for r in S]); 
Tlo=np.array([ (1/r['HR'][0])**(1/r['HR'][1]['T3']) for r in S])*1e3
Thi=np.array([ (r['Mach'][0])**(-1/r['Mach'][1]['T3']) for r in S])*1e3
rho0=np.array([r['rho0'] for r in S]); G=np.array([r['G'] for r in S]); sig=np.array([r['sig'] for r in S])
f_rho8=np.interp(np.log10(1e8),np.log10(rho0),fs); f_rho9=np.interp(np.log10(1e9),np.log10(rho0),fs)
fig,(axL,axR)=plt.subplots(1,2,figsize=(9.6,3.7),gridspec_kw={'wspace':0.32})
ff=np.linspace(fs[0],fs[-1],400); Tlo_i=np.interp(ff,fs,Tlo); Thi_i=np.interp(ff,fs,Thi)
m8=ff<=f_rho8; m9=(ff>=f_rho8)&(ff<=f_rho9)
axL.fill_between(ff[m8],Tlo_i[m8],Thi_i[m8],color=GREEN,alpha=0.55,lw=0,label=r'allowed ($\rho_0<10^{8}\,M_\odot\,{\rm pc}^{-3}$)')
axL.fill_between(ff[m9],Tlo_i[m9],Thi_i[m9],color=GREEN,alpha=0.22,lw=0,label=r'allowed if $\rho_0<10^{9}\,M_\odot\,{\rm pc}^{-3}$')
axL.plot(fs,Tlo,color=BLUE,lw=1.8,label=r'$H=R_{\rm MC}$ (lower bound)')
axL.plot(fs,Thi,color=RED,lw=1.8,label=r'$\mathcal{M}=1$ (upper bound)')
axL.axvline(f_rho8,color=GOLD,lw=1.6,ls='--')
axL.axvline(f_rho9,color=GOLD,lw=1.6,ls='--')
axL.set_xlim(0.1,0.6); axL.set_ylim(500,1800)
axL.set_xlabel(r'$f_\ast$',fontsize=13); axL.set_ylabel(r'$T_{\rm MC}$ (K)',fontsize=12)
axL.legend(fontsize=8,loc='upper left',**LEG)
axL.set_title(r'allowed region, $M_{\rm h}=10^{6}\,M_\odot$',fontsize=11,pad=6)
axL.scatter([0.22],[900],marker='*',s=90,color='k',zorder=6); axL.text(0.235,905,'fiducial',fontsize=8.5,va='center')
axR.plot(fs,sig,'o-',color=BLUE,ms=5,lw=1.6); axR.set_ylabel(r'$\sigma_{\rm eq}$ (km s$^{-1}$)',fontsize=12,color=BLUE)
axR.tick_params(axis='y',colors=BLUE); axR.set_xlabel(r'$f_\ast$',fontsize=13); axR.set_xlim(0.1,0.6)
axR.set_title('equilibrium across the range (metallicity-independent)',fontsize=10.5,pad=6)
axG=axR.twinx(); axG.plot(fs,G,'s--',color=RED,ms=5,lw=1.6); axG.set_yscale('log')
axG.set_ylabel(r'$\Gamma_{\rm eq}$ (yr$^{-1}$)',fontsize=12,color=RED); axG.tick_params(colors=RED)
fig.tight_layout(pad=0.4); fig.savefig(os.path.join(PAPER,'fplane.pdf')); plt.close(fig)
print('f at rho0=1e8:',f_rho8,' 1e9:',f_rho9)
print('T windows (f, Tlo, Thi):'); [print(f, round(a), round(c)) for f,a,c in zip(fs,Tlo,Thi)]
# mass dependence of window
e=b['H/R_MC'][1]; eM=b['Mach'][1]
for m6 in [0.7,1,2,3,5]:
    tl=(1/ev('H/R_MC',m6=m6))**(1/e['T3'])*1e3; th=(ev('Mach',m6=m6))**(-1/eM['T3'])*1e3
    print(f"m6={m6}: T window {tl:.0f}-{th:.0f} K; n_MC={ev('n_MC[cm-3]',m6=m6):.2e}; Gamma={R['C3']['gp']*m6**R['C3']['gs']:.2e}; Mdisk(Tmid)={ev('M_disk[Msun]',m6=m6,T3=np.sqrt(tl*th)/1e3):.2e}; AV(Z=0.1,Tmid)={0.1*ev('A_V(Zsun)',m6=m6,T3=np.sqrt(tl*th)/1e3):.2f}; fEdd={ev('L_AGN/L_Edd',m6=m6,T3=np.sqrt(tl*th)/1e3):.3f}")
