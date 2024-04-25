import numpy as np
import matplotlib.pyplot as plt
import os

# load the .npy files with the saved numpy arrays
# Bstar=0,1500,2000,2500
Bstarlist=[0,1500,2000,2500]
# 'SEMF','FRDM2012','WS4RBF','BSKG3','HFB31','DDMED','DDPC1','SVMIN','UNEDF0'
modellist=['SEMF','FRDM2012','WS4RBF','BSKG3','HFB31','DDMED','DDPC1','SVMIN','UNEDF0']

def keyname(prefix, model, field):
    
    '''
    return a key of the dictionary data 
    
    '''
    return f'{prefix}_m{model}_f{field}'

data={}

# read the .npy files and assign the value to variables having the same name as
# the filename without the extension .npy
for filename in os.listdir():
    if filename.endswith('.npy'):
        array=np.load(filename)
        key=filename[:-4]
        data[key]=array

del array, filename, key

#%% plot1
plt.figure()
for m in modellist:
    plt.plot(Bstarlist, [data[keyname('presslist',m,f)][-1] for f in Bstarlist],
               'o-', label=m)
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.grid()
plt.legend()
plt.xlabel(r'Magnetic field strength $B_{\star}$')
plt.ylabel(r'Pressure at neutron-drip transition [$MeVfm^{-3}$]')
plt.title(r'$P_{drip}$ vs $B_{\star}$')
plt.show()   
plt.savefig('Pdrip_vs_Bstar.png')

#%% plot2, now for the same model, plot for all Bstar
for m in modellist:
    plt.figure()
    for f in Bstarlist:
        plt.loglog(data[keyname('nbminlist',m,f)],
                   data[keyname('presslist',m,f)], '-',
                   linewidth=1, label=rf'$B_*=${f}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Baryon number density [$fm^{-3}$]')
    plt.ylabel(r'Pressure [$MeVfm^{-3}$]')
    plt.title(rf'P vs $n_b$ (model: {m})')
    plt.show()   
    plt.savefig(f'Ptot_vs_nb_m{m}.png')

# now for the same Bstar, plot for all models
for f in Bstarlist:
    plt.figure()
    for m in modellist:
        plt.loglog(data[keyname('nbminlist',m,f)],
                   data[keyname('presslist',m,f)], '-',
                   linewidth=1, label=rf'model: {m}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Baryon number density [$fm^{-3}$]')
    plt.ylabel(r'Pressure [$MeVfm^{-3}$]')
    plt.title(rf'P vs $n_b$ ($B_*=${f})')
    plt.show()   
    plt.savefig(f'Ptot_vs_nb_f{f}.png')

#%% plot3, now for the same model, plot for all Bstar
for m in modellist:
    plt.figure()
    for f in Bstarlist:
        plt.loglog(data[keyname('etotminlist',m,f)],
                   data[keyname('presslist',m,f)], '-',
                   linewidth=1, label=rf'$B_*=${f}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
    plt.ylabel(r'Pressure [$MeVfm^{-3}$]')
    plt.title(rf'P vs $\epsilon$ (model: {m})')
    plt.show()   
    plt.savefig(f'Ptot_vs_etot_m{m}.png')

# now for the same Bstar, plot for all models
for f in Bstarlist:
    plt.figure()
    for m in modellist:
        plt.loglog(data[keyname('etotminlist',m,f)],
                   data[keyname('presslist',m,f)], '-',
                   linewidth=1, label=rf'model: {m}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
    plt.ylabel(r'Pressure [$MeVfm^{-3}$]')
    plt.title(rf'P vs $\epsilon$ ($B_*=${f})')
    plt.show()   
    plt.savefig(f'Ptot_vs_etot_f{f}.png')

#%% plot4, now for the same model, plot for all Bstar
for m in modellist:
    plt.figure()
    for f in Bstarlist:
        plt.loglog(data[keyname('etotminlist',m,f)],
                   data[keyname('shearlist',m,f)], '-',
                   linewidth=1, label=rf'$B_*=${f}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
    plt.ylabel(r'Effective shear modulus [$MeVfm^{-3}$]')
    plt.title(rf'S (shear modulus) vs $\epsilon$ (model: {m})')
    plt.show()   
    plt.savefig(f'S_vs_etot_m{m}.png')

# now for the same Bstar, plot for all models
for f in Bstarlist:
    plt.figure()
    for m in modellist:
        plt.loglog(data[keyname('etotminlist',m,f)],
                   data[keyname('shearlist',m,f)], '-',
                   linewidth=1, label=rf'model: {m}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
    plt.ylabel(r'Effective shear modulus [$MeVfm^{-3}$]')
    plt.title(rf'S (shear modulus) vs $\epsilon$ ($B_*=${f})')
    plt.show()   
    plt.savefig(f'S_vs_etot_f{f}.png')

#%% plot5, now for the same model, plot for all Bstar
for m in modellist:
    plt.figure()
    for f in Bstarlist:
        plt.loglog(data[keyname('nbminlist',m,f)],
                   data[keyname('shearlist',m,f)], '-',
                   linewidth=1, label=rf'$B_*=${f}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Baryon number density [$fm^{-3}$]')
    plt.ylabel(r'Effective shear modulus [$MeVfm^{-3}$]')
    plt.title(rf'S (shear modulus) vs $n_b$ (model: {m})')
    plt.show()   
    plt.savefig(f'S_vs_nb_m{m}.png')

# now for the same Bstar, plot for all models
for f in Bstarlist:
    plt.figure()
    for m in modellist:
        plt.loglog(data[keyname('nbminlist',m,f)],
                   data[keyname('shearlist',m,f)], '-',
                   linewidth=1, label=rf'model: {m}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Baryon number density [$fm^{-3}$]')
    plt.ylabel(r'Effective shear modulus [$MeVfm^{-3}$]')
    plt.title(rf'S (shear modulus) vs $n_b$ ($B_*=${f})')
    plt.show()   
    plt.savefig(f'S_vs_nb_f{f}.png')

#%% plot6, now for the same model, plot for all Bstar
for m in modellist:
    plt.figure()
    for f in Bstarlist:
        plt.semilogx(data[keyname('nbminlist',m,f)],
                   data[keyname('gamma',m,f)], '-',
                   linewidth=1, label=rf'$B_*=${f}') 
        if f>1304:
            plt.ylim(1, 3.5)    
    plt.grid()
    plt.axhline (y=4/3, color='black', linestyle='--', linewidth=1)
    plt.legend()
    plt.xlabel(r'Baryon number density [$fm^{-3}$]')
    plt.ylabel(r'Adiabatic index')
    plt.title(rf'$\Gamma$ vs $n_b$ (model: {m})')
    plt.show()   
    plt.savefig(f'Gamma_vs_nb_m{m}.png')

# now for the same Bstar, plot for all models
for f in Bstarlist:
    plt.figure()
    for m in modellist:
        plt.semilogx(data[keyname('nbminlist',m,f)],
                   data[keyname('gamma',m,f)], '-',
                   linewidth=1, label=rf'model: {m}') 
        if f>1304:
            plt.ylim(1.5, 3.5)    
    plt.grid()
    plt.axhline (y=4/3, color='black', linestyle='--', linewidth=1)
    plt.legend()
    plt.xlabel(r'Baryon number density [$fm^{-3}$]')
    plt.ylabel(r'Adiabatic index')
    plt.title(rf'$\Gamma$ vs $n_b$ ($B_*=${f})')
    plt.show()   
    plt.savefig(f'Gamma_vs_nb_f{f}.png')

#%% plot7, now for the same model, plot for all Bstar
for m in modellist:
    plt.figure()
    for f in Bstarlist:
        plt.semilogx(data[keyname('etotminlist',m,f)],
                   data[keyname('gamma',m,f)], '-',
                   linewidth=1, label=rf'$B_*=${f}') 
        if f>1304:
            plt.ylim(1, 3.5)        
    plt.grid()
    plt.axhline (y=4/3, color='black', linestyle='--', linewidth=1)
    plt.legend()
    plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
    plt.ylabel(r'Adiabatic index')
    plt.title(rf'$\Gamma$ vs $\epsilon$ (model: {m})')
    plt.show()   
    plt.savefig(f'Gamma_vs_etot_m{m}.png')

# now for the same Bstar, plot for all models
for f in Bstarlist:
    plt.figure()
    for m in modellist:
        plt.semilogx(data[keyname('etotminlist',m,f)],
                   data[keyname('gamma',m,f)], '-',
                   linewidth=1, label=rf'model: {m}') 
        if f>1304:
            plt.ylim(1.5, 3.5)      
    plt.grid()
    plt.axhline (y=4/3, color='black', linestyle='--', linewidth=1)
    plt.legend()
    plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
    plt.ylabel(r'Adiabatic index')
    plt.title(rf'$\Gamma$ vs $\epsilon$ ($B_*=${f})')
    plt.show()   
    plt.savefig(f'Gamma_vs_etot_f{f}.png')

#%% plot8, now for the same model, plot for all Bstar
for m in modellist:
    plt.figure()
    for f in Bstarlist:
        plt.semilogx(data[keyname('nbminlist',m,f)],
                   data[keyname('soundspeed',m,f)], '-',
                   linewidth=1, label=rf'$B_*=${f}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Baryon number density [$fm^{-3}$]')
    plt.ylabel(r'Speed of sound [c]')
    plt.title(rf'$c_s/c$ vs $n_b$ (model: {m})')
    plt.show()   
    plt.savefig(f'speed_vs_nb_m{m}.png')

# now for the same Bstar, plot for all models
for f in Bstarlist:
    plt.figure()
    for m in modellist:
        plt.semilogx(data[keyname('nbminlist',m,f)],
                   data[keyname('soundspeed',m,f)], '-',
                   linewidth=1, label=rf'model: {m}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Baryon number density [$fm^{-3}$]')
    plt.ylabel(r'Speed of sound [c]')
    plt.title(rf'$c_s/c$ vs $n_b$ ($B_*=${f})')
    plt.show()   
    plt.savefig(f'speed_vs_nb_f{f}.png')

#%% plot9, now for the same model, plot for all Bstar
for m in modellist:
    plt.figure()
    for f in Bstarlist:
        plt.semilogx(data[keyname('etotminlist',m,f)],
                   data[keyname('soundspeed',m,f)], '-',
                   linewidth=1, label=rf'$B_*=${f}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
    plt.ylabel(r'Speed of sound [c]')
    plt.title(rf'$c_s/c$ vs $\epsilon$ (model: {m})')
    plt.show()   
    plt.savefig(f'speed_vs_etot_m{m}.png')

# now for the same Bstar, plot for all models
for f in Bstarlist:
    plt.figure()
    for m in modellist:
        plt.semilogx(data[keyname('etotminlist',m,f)],
                   data[keyname('soundspeed',m,f)], '-',
                   linewidth=1, label=rf'model: {m}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
    plt.ylabel(r'Speed of sound [c]')
    plt.title(rf'$c_s/c$ vs $\epsilon$ ($B_*=${f})')
    plt.show()   
    plt.savefig(f'speed_vs_etot_f{f}.png')

#%% plot10, now for the same model, plot for all Bstar
for m in modellist:
    plt.figure()
    for f in Bstarlist:
        plt.loglog(data[keyname('etotminlist',m,f)],
                   data[keyname('compression',m,f)], '-',
                   linewidth=1, label=rf'$B_*=${f}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
    plt.ylabel(r'Compression modulus [$MeVfm^{-3}$]')
    plt.title(rf'$K$ (compression modulus) vs $\epsilon$ (model: {m})')
    plt.show()   
    plt.savefig(f'compression_vs_etot_m{m}.png')

# now for the same Bstar, plot for all models
for f in Bstarlist:
    plt.figure()
    for m in modellist:
        plt.loglog(data[keyname('etotminlist',m,f)],
                   data[keyname('compression',m,f)], '-',
                   linewidth=1, label=rf'model: {m}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
    plt.ylabel(r'Compression modulus [$MeVfm^{-3}$]')
    plt.title(rf'$K$ (compression modulus) vs $\epsilon$ ($B_*=${f})')
    plt.show()   
    plt.savefig(f'compression_vs_etot_f{f}.png')

#%% plot11, now for the same model, plot for all Bstar
for m in modellist:
    plt.figure()
    for f in Bstarlist:
        plt.loglog(data[keyname('nbminlist',m,f)],
                   data[keyname('compression',m,f)], '-',
                   linewidth=1, label=rf'$B_*=${f}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Baryon number density [$fm^{-3}$]')
    plt.ylabel(r'Compression modulus [$MeVfm^{-3}$]')
    plt.title(rf'$K$ (compression modulus) vs $n_b$ (model: {m})')
    plt.show()   
    plt.savefig(f'compression_vs_nb_m{m}.png')

# now for the same Bstar, plot for all models
for f in Bstarlist:
    plt.figure()
    for m in modellist:
        plt.loglog(data[keyname('nbminlist',m,f)],
                   data[keyname('compression',m,f)], '-',
                   linewidth=1, label=rf'model: {m}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Baryon number density [$fm^{-3}$]')
    plt.ylabel(r'Compression modulus [$MeVfm^{-3}$]')
    plt.title(rf'$K$ (compression modulus) vs $n_b$ ($B_*=${f})')
    plt.show()   
    plt.savefig(f'compression_vs_nb_f{f}.png')

#%% plot12, now for the same model, plot for all Bstar
for m in modellist:
    plt.figure()
    for f in Bstarlist:
        plt.semilogx(data[keyname('nbminlist',m,f)],
                   data[keyname('shearlist',m,f)]/data[keyname('compression',m,f)], '-',
                   linewidth=1, label=rf'$B_*=${f}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Baryon number density [$fm^{-3}$]')
    plt.ylabel(r'$\frac{shear}{compression}$')
    plt.title(rf'S/K vs $n_b$ (model: {m})')
    plt.show()   
    plt.savefig(f'SK_vs_nb_m{m}.png')

# now for the same Bstar, plot for all models
for f in Bstarlist:
    plt.figure()
    for m in modellist:
        plt.semilogx(data[keyname('nbminlist',m,f)],
                   data[keyname('shearlist',m,f)]/data[keyname('compression',m,f)], '-',
                   linewidth=1, label=rf'model: {m}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Baryon number density [$fm^{-3}$]')
    plt.ylabel(r'$\frac{shear}{compression}$')
    plt.title(rf'S/K vs $n_b$ ($B_*=${f})')
    plt.show()   
    plt.savefig(f'SK_vs_nb_f{f}.png')

#%% plot13, now for the same model, plot for all Bstar
for m in modellist:
    plt.figure()
    for f in Bstarlist:
        plt.semilogx(data[keyname('etotminlist',m,f)],
                   data[keyname('shearlist',m,f)]/data[keyname('compression',m,f)], '-',
                   linewidth=1, label=rf'$B_*=${f}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
    plt.ylabel(r'$\frac{shear}{compression}$')
    plt.title(rf'S/K vs $\epsilon$ (model: {m})')
    plt.show()   
    plt.savefig(f'SK_vs_etot_m{m}.png')

# now for the same Bstar, plot for all models
for f in Bstarlist:
    plt.figure()
    for m in modellist:
        plt.semilogx(data[keyname('etotminlist',m,f)],
                   data[keyname('shearlist',m,f)]/data[keyname('compression',m,f)], '-',
                   linewidth=1, label=rf'model: {m}') 
    plt.grid()
    plt.legend()
    plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
    plt.ylabel(r'$\frac{shear}{compression}$')
    plt.title(rf'S/K vs $\epsilon$ ($B_*=${f})')
    plt.show()   
    plt.savefig(f'SK_vs_etot_f{f}.png')

#%% plot14, now for the same model, plot for all Bstar
for m in modellist:
    plt.figure(figsize=(8,6))
    plt.suptitle(rf'model: {m}')
    plt.subplot(3,1,1)
    for f in Bstarlist:
        plt.semilogx(data[keyname('nbminlist',m,f)],
                     data[keyname('Zminlist',m,f)], '-',
                     linewidth=1, label=rf'$B_*=${f}')
    plt.ylabel(r'Z')
    plt.grid()
    plt.legend()
    plt.subplot(3,1,2)
    for f in Bstarlist:
        plt.semilogx(data[keyname('nbminlist',m,f)],
                     data[keyname('Nminlist',m,f)], '-',
                     linewidth=1, label=rf'$B_*=${f}')
    plt.ylabel(r'N')
    plt.grid()
    plt.legend()
    plt.subplot(3,1,3)
    for f in Bstarlist:
        plt.semilogx(data[keyname('nbminlist',m,f)],
                     data[keyname('ZoverA',m,f)], '-',
                     linewidth=1, label=rf'$B_*=${f}')
    plt.xlabel(r'Baryon number density [$fm^{-3}$]')
    plt.ylabel(r'Z/A')
    plt.grid()
    plt.legend()
    plt.show()   
    plt.savefig(f'Z_N_ZoverA_vs_nb_m{m}.png')

# now for the same Bstar, plot for all models
for f in Bstarlist:
    plt.figure(figsize=(8,6))
    plt.suptitle(rf'$B_*=${f}')
    plt.subplot(3,1,1)
    for m in modellist:
        plt.semilogx(data[keyname('nbminlist',m,f)],
                     data[keyname('Zminlist',m,f)], '-',
                     linewidth=1, label=rf'model: {m}')
    plt.ylabel(r'Z')
    plt.grid()
    plt.legend()
    plt.subplot(3,1,2)
    for m in modellist:
        plt.semilogx(data[keyname('nbminlist',m,f)],
                     data[keyname('Nminlist',m,f)], '-',
                     linewidth=1, label=rf'model: {m}')
    plt.ylabel(r'N')
    plt.grid()
    plt.legend()
    plt.subplot(3,1,3)
    for m in modellist:
        plt.semilogx(data[keyname('nbminlist',m,f)],
                     data[keyname('ZoverA',m,f)], '-',
                     linewidth=1, label=rf'model: {m}')
    plt.xlabel(r'Baryon number density [$fm^{-3}$]')
    plt.ylabel(r'Z/A')
    plt.grid()
    plt.legend()
    plt.show()   
    plt.savefig(f'Z_N_ZoverA_vs_nb_f{f}.png')

#%% plot15, now for the same model, plot for all Bstar
for m in modellist:
    plt.figure(figsize=(8,6))
    plt.suptitle(rf'model: {m}')
    plt.subplot(3,1,1)
    for f in Bstarlist:
        plt.semilogx(data[keyname('etotminlist',m,f)],
                     data[keyname('Zminlist',m,f)], '-',
                     linewidth=1, label=rf'$B_*=${f}')
    plt.ylabel(r'Z')
    plt.grid()
    plt.legend()
    plt.subplot(3,1,2)
    for f in Bstarlist:
        plt.semilogx(data[keyname('etotminlist',m,f)],
                     data[keyname('Nminlist',m,f)], '-',
                     linewidth=1, label=rf'$B_*=${f}')
    plt.ylabel(r'N')
    plt.grid()
    plt.legend()
    plt.subplot(3,1,3)
    for f in Bstarlist:
        plt.semilogx(data[keyname('etotminlist',m,f)],
                     data[keyname('ZoverA',m,f)], '-',
                     linewidth=1, label=rf'$B_*=${f}')
    plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
    plt.ylabel(r'Z/A')
    plt.grid()
    plt.legend()
    plt.show()   
    plt.savefig(f'Z_N_ZoverA_vs_etot_m{m}.png')

# now for the same Bstar, plot for all models
for f in Bstarlist:
    plt.figure(figsize=(8,6))
    plt.suptitle(rf'$B_*=${f}')
    plt.subplot(3,1,1)
    for m in modellist:
        plt.semilogx(data[keyname('etotminlist',m,f)],
                     data[keyname('Zminlist',m,f)], '-',
                     linewidth=1, label=rf'model: {m}')
    plt.ylabel(r'Z')
    plt.grid()
    plt.legend()
    plt.subplot(3,1,2)
    for m in modellist:
        plt.semilogx(data[keyname('etotminlist',m,f)],
                     data[keyname('Nminlist',m,f)], '-',
                     linewidth=1, label=rf'model: {m}')
    plt.ylabel(r'N')
    plt.grid()
    plt.legend()
    plt.subplot(3,1,3)
    for m in modellist:
        plt.semilogx(data[keyname('etotminlist',m,f)],
                     data[keyname('ZoverA',m,f)], '-',
                     linewidth=1, label=rf'model: {m}')
    plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
    plt.ylabel(r'Z/A')
    plt.grid()
    plt.legend()
    plt.show()   
    plt.savefig(f'Z_N_ZoverA_vs_etot_f{f}.png')

#%% plot16, now for the same model, plot for all Bstar
for m in modellist:
    plt.figure(figsize=(8,6))
    for f in Bstarlist:
        plt.plot(data[keyname('Rlist',m,f)], data[keyname('Mcr1_2',m,f)],
                     'o-', label=rf'M=1.2, $B_*=${f}')
        plt.plot(data[keyname('Rlist',m,f)], data[keyname('Mcr1_4',m,f)],
                 'o-', label=rf'M=1.4, $B_*=${f}')
        plt.plot(data[keyname('Rlist',m,f)], data[keyname('Mcr1_6',m,f)],
                 'o-', label=rf'M=1.6, $B_*=${f}')
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    plt.grid()
    plt.legend()
    plt.xlabel(r'Total radius of NS [km]')
    plt.ylabel(r'Grav. mass of outer crust [$M_{\odot}$]')
    plt.title(rf'Mcr vs R (model: {m})')
    plt.show()   
    plt.savefig(f'Mcr_m{m}.png')

# now for the same Bstar, plot for all models
for f in Bstarlist:
    plt.figure(figsize=(8,6))
    for m in modellist:
        plt.plot(data[keyname('Rlist',m,f)], data[keyname('Mcr1_2',m,f)],
                     'o-', label=rf'M=1.2, model: {m}')
        plt.plot(data[keyname('Rlist',m,f)], data[keyname('Mcr1_4',m,f)],
                 'o-', label=rf'M=1.4, model: {m}')
        plt.plot(data[keyname('Rlist',m,f)], data[keyname('Mcr1_6',m,f)],
                 'o-', label=rf'M=1.6, model: {m}')
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    plt.grid()
    plt.legend()
    plt.xlabel(r'Total radius of NS [km]')
    plt.ylabel(r'Grav. mass of outer crust [$M_{\odot}$]')
    plt.title(rf'Mcr vs R ($B_*=${f})')
    plt.show()   
    plt.savefig(f'Mcr_f{f}.png')

#%% plot17, now for the same model, plot for all Bstar
for m in modellist:
    plt.figure(figsize=(8,6))
    for f in Bstarlist:
        plt.plot(data[keyname('Rlist',m,f)], data[keyname('zdrip1_2',m,f)],
                     'o-', label=rf'M=1.2, $B_*=${f}')
        plt.plot(data[keyname('Rlist',m,f)], data[keyname('zdrip1_4',m,f)],
                 'o-', label=rf'M=1.4, $B_*=${f}')
        plt.plot(data[keyname('Rlist',m,f)], data[keyname('zdrip1_6',m,f)],
                 'o-', label=rf'M=1.6, $B_*=${f}')
    plt.grid()
    plt.legend()
    plt.xlabel(r'Total radius of NS [km]')
    plt.ylabel(r'Proper thickness of outer crust [km]')
    plt.title(rf'$z_d$ vs R (model: {m})')
    plt.show()   
    plt.savefig(f'zdrip_m{m}.png')

# now for the same Bstar, plot for all models
for f in Bstarlist:
    plt.figure(figsize=(8,6))
    for m in modellist:
        plt.plot(data[keyname('Rlist',m,f)], data[keyname('zdrip1_2',m,f)],
                     'o-', label=rf'M=1.2, model: {m}')
        plt.plot(data[keyname('Rlist',m,f)], data[keyname('zdrip1_4',m,f)],
                 'o-', label=rf'M=1.4, model: {m}')
        plt.plot(data[keyname('Rlist',m,f)], data[keyname('zdrip1_6',m,f)],
                 'o-', label=rf'M=1.6, model: {m}')
    plt.grid()
    plt.legend()
    plt.xlabel(r'Total radius of NS [km]')
    plt.ylabel(r'Proper thickness of outer crust [km]')
    plt.title(rf'$z_d$ vs R ($B_*=${f})')
    plt.show()   
    plt.savefig(f'zdrip_f{f}.png')
