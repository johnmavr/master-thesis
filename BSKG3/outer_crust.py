import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
import fortranformat 
import time

# parameters of beSEMF from https://www.actaphys.uj.edu.pl/R/37/6/1833/pdf
AV, AS, AC, AA, AP=15.777, 18.341, 0.710, 23.211, 11.996
# definition of useful constants 
Mn=939.565421 # neutron mass MeV/c**2
Mp=938.272088 # proton mass MeV/c**2
Me=0.510999 # electron mass MeV/c**2
AMU=931.494102 # atomic mass unit MeV/c**2
HBARC=197.326980 # MeVfm
ALPHA=7.297353e-3 # fine structure constant
LAMDAe=386.159268 # reduced electron Compton wavelength fm
K1=Me*LAMDAe**(-3) # MeVfm**(-3)
CM=0.895929255682 # Madelung constant for bcc lattice
# conversion of units
MeV_fm3to_km2=1.3234e-6 # 1MeVfm**(-3)=1.3234e-6 km**(-2)
MeV_c2toMsolar=8.964992783e-61 # 1MeV/c**2=8.964992783e-61 Msolar

def beSEMF(A,Z):
    
    '''
    return the nuclear binding energy for SEMF model, given mass number A and 
    atomic number Z, in MeV
    
    '''
    N=A-Z
    
    be=AV*A-AS*A**(2/3)-AC*Z*(Z-1)*A**(-1/3)-AA*((A-2*Z)**2)/A
    
    if (Z%2==0) and (N%2==0):
        be=be+AP*A**(-1/2)
    elif (Z%2!=0) and (N%2!=0):
        be=be-AP*A**(-1/2)
    
    return be

def Sp(A,Z):
    
    '''
    return the one-proton separation energy
    
    '''
    return beSEMF(A,Z)-beSEMF(A-1,Z-1)

def S2p(A,Z):
    
    '''
    return the two-proton separation energy
    
    '''
    return beSEMF(A,Z)-beSEMF(A-2,Z-2)

def Sn(A,Z):
    
    '''
    return the one-neutron separation energy
    
    '''
    return beSEMF(A,Z)-beSEMF(A-1,Z)

def S2n(A,Z):
    
    '''
    return the two-neutron separation energy
    
    '''
    return beSEMF(A,Z)-beSEMF(A-2,Z)

def Bel(Z):
    
    '''
    return the binding energy of atomic electrons in MeV
    
    '''
    return 1.44381e-5*Z**(2.39)+1.55468e-12*Z**(5.35)

def WN(A,Z,tabular_form,external=None):
    
    '''
    return the mass of nucleus in MeV
    
    '''
    if tabular_form=='semf_be':
        return Z*Mp+(A-Z)*Mn-beSEMF(A,Z)
    elif tabular_form=='be':
        return Z*Mp+(A-Z)*Mn+external # be is tabulated as a negative number, otherwise change here the sign to -
    elif tabular_form=='mass_excess':
        return external+A*AMU-Z*Me+Bel(Z)    

def ge(xe):
    
    ''' 
    return the nondimensional Fermi energy (chemical potential)
    of electrons ge=mue/(Me*c**2), (xe is the nondimensional Fermi momentum xe=peF/(Me*c))
    
    '''
    return (xe**2+1)**(1/2)

def mue(ge):
    
    '''
    return the chemical potential of electrons in MeV
    
    '''
    return Me*ge
    
def ne(xe,B):
    
    ''' 
    return the number density of electrons in fm**(-3)
    
    '''
    if B==0:
        return (xe**3)/(3*np.pi**2*LAMDAe**3)
    else:
        return (B*xe)/(2*np.pi**2*LAMDAe**3)

def nb(A,Z,xe,B):
    
    ''' 
    return the number density of baryons in fm**(-3)
    
    '''
    return (A/Z)*ne(xe,B)

def nN(Z,xe,B):
    
    ''' 
    return the number density of nuclei in fm**(-3)
    
    '''
    return ne(xe,B)/Z

def rcell(Z,xe,B):
    
    ''' 
    return the Rcell in fm

    '''
    return (4*np.pi*nN(Z,xe,B)/3)**(-1/3)

def shear(Z,xe,B):
    
    ''' 
    return the effective shear modulus in MeVfm**(-3)

    '''
    return 0.1194*nN(Z,xe,B)*Z**2*ALPHA*HBARC/rcell(Z,xe,B)

def f(x,B):
    
    '''
    function useful for electron energy density
    
    '''
    if B==0:
        return (1/(8*np.pi**2))*(x*(1+x**2)**(1/2)*(1+2*x**2)-np.log(x+(1+x**2)**(1/2)))  
    else:
        return x*(1+x**2)**(1/2)+np.log(x+(1+x**2)**(1/2))
    
def g(x,B): 
                   
    '''
    function useful for electron pressure
    
    '''
    if B==0:
        return (1/(8*np.pi**2))*(x*(1+x**2)**(1/2)*((2/3)*x**2-1)+np.log(x+(1+x**2)**(1/2)))  
    else:
        return x*(1+x**2)**(1/2)-np.log(x+(1+x**2)**(1/2))
    
def ee(xe,B):
    
    ''' 
    return the energy density of electrons in MeVfm**(-3)

    '''
    if B==0:
        return K1*f(xe,B)
    else:
        return (K1*B*f(xe,B))/((2*np.pi)**2)

def Pe(xe,B):
    
    ''' 
    return the pressure of electrons in MeVfm**(-3)

    '''
    if B==0:
        return K1*g(xe,B)
    else:
        return (K1*B*g(xe,B))/((2*np.pi)**2)

def el(Z,xe,B):
    
    ''' 
    return the energy density of lattice in MeVfm**(-3)

    '''
    return -CM*(4*np.pi/3)**(1/3)*Z**(2/3)*ALPHA*HBARC*(ne(xe,B))**(4/3)

def Pl(Z,xe,B):
    
    ''' 
    return the pressure of lattice in MeVfm**(-3)

    '''
    return el(Z,xe,B)/3

def Ptot(Z,xe,B):
    
    ''' 
    return the total pressure in MeVfm**(-3)

    '''
    return Pe(xe,B)+Pl(Z,xe,B)

def etot(A,Z,xe,B,tabular_form,external=None):
    
    ''' 
    return the total energy density in MeVfm**(-3)

    '''
    return nN(Z,xe,B)*WN(A,Z,tabular_form,external)+ee(xe,B)+el(Z,xe,B)

def gibbs(A,Z,xe,B,tabular_form,external=None):
    
    ''' 
    return Gcell/A (chemical potential of baryons) in MeV

    '''
    return (etot(A,Z,xe,B,tabular_form,external)+Ptot(Z,xe,B))/nb(A,Z,xe,B)

def funtosolve(xe,Z,Pfix,B):
    
    '''   
    for fixed Z and Pfix, find numerically the xroot of this function
    
    '''    
    return Ptot(Z,xe,B)-Pfix

#%%--------------------------------USER GIVES values----------------------------------------------------------------

model='BSKG3' # GIVE the name of the nuclear model: 'SEMF' or ...
tabulation='mass_excess' # GIVE the tabular_form: 'semf_be' or 'be' or 'mass_excess'
Bstar=0 # GIVE a magnetic field in Bcr units (Bcr=4.414e13G): Bstar>1304 or Bstar=0
filename='BSKG3_table.dat' # GIVE the filename of the nuclear mass table: 'FRDM2012_table.dat' or ..., or None for SEMF
extra_tabulation=['Z','A'] # GIVE what the header of file contains: (['Z','A'], if contains only Z and A in any order) or (['Z','N'], if contains only Z and N in any order) or (['Z','N','A'], if contains only Z, N and A in any order), or None for SEMF
sepform='\\s+' # GIVE how the columns are saparated in the file: ' ' or '\\s+' or ',' or ..., or None for SEMF

#%% set-up of variables for minimization process
field=(f'{Bstar}').replace(".", "_")

# the endpoint must be bigger than Pdrip~5e-4 MeVfm**(-3) (without magnetic field)
pressinterval=(6.082e-15, 20.0e-4) 
numberpoints=18000 
presslist=np.logspace(np.log10(pressinterval[0]), np.log10(pressinterval[1]),
                      numberpoints, dtype=np.float64)

Aminlist=np.zeros((numberpoints), dtype=np.int16)
Zminlist=np.zeros((numberpoints), dtype=np.int16)
Nminlist=np.zeros((numberpoints), dtype=np.int16)
xeminlist=np.zeros((numberpoints), dtype=np.float64)
gminlist=np.zeros((numberpoints), dtype=np.float64)
etotminlist=np.zeros((numberpoints), dtype=np.float64)

Zmin_interest=24
Zmax_interest=56
Nmin_interest=24
Nmax_interest=132

reader=fortranformat.FortranRecordReader('(a1,i3,i5,i5,i5,1x,a3,a4,1x,f14.6)')

with open('AME2020_table.dat','r') as fileAME:
    
    lines=fileAME.readlines()
    
dfexp=pd.DataFrame(columns=['cc', 'N-Z', 'N', 'Z', 'A', 'el', 'origin', 'mass_excess[keV]'])
        
for index in range(1,len(lines)):
        
    try:
        dfexp.loc[len(dfexp)]=reader.read(lines[index]) 
    except ValueError as message:
        print(message,'\nIt is not experimentally measured mass\n')
        continue

del lines
dfexp=dfexp.drop(['cc','N-Z','el','origin'], axis=1)
dfexp=dfexp.sort_values(by=['Z','N'])
dfexp_interest=dfexp.query(f'{Zmin_interest}<=Z<={Zmax_interest} and {Nmin_interest}<=N<={Nmax_interest}')
Zexp_array=dfexp_interest['Z'].to_numpy()
Aexp_array=dfexp_interest['A'].to_numpy()
externalexp_array=dfexp_interest['mass_excess[keV]'].to_numpy()    
externalexp_array=externalexp_array*10**(-3) # exp. mass excess in MeV

#%% start of gibbs minimization process
i=0

if (tabulation=='be') or (tabulation=='mass_excess'):
    
    if extra_tabulation==['Z', 'A']:
        
        df=pd.read_csv(filename, sep=sepform, header=0, usecols=['Z','A',tabulation])
        df=df.sort_values(by=['Z','A'])
        df['N']=df['A']-df['Z']
        df_interest=df.query(f'{Zmin_interest}<=Z<={Zmax_interest} and {Nmin_interest}<=N<={Nmax_interest}')
        Z_array=df_interest['Z'].to_numpy()
        A_array=df_interest['A'].to_numpy()
        external_array=df_interest[tabulation].to_numpy()
        
    elif extra_tabulation==['Z', 'N']:
        
        df=pd.read_csv(filename, sep=sepform, header=0, usecols=['Z','N',tabulation])
        df=df.sort_values(by=['Z','N'])
        df['A']=df['Z']+df['N']
        df_interest=df.query(f'{Zmin_interest}<=Z<={Zmax_interest} and {Nmin_interest}<=N<={Nmax_interest}')
        Z_array=df_interest['Z'].to_numpy()
        A_array=df_interest['A'].to_numpy()
        external_array=df_interest[tabulation].to_numpy()
        
    elif extra_tabulation==['Z', 'N', 'A']:
        
        df=pd.read_csv(filename, sep=sepform, header=0, usecols=['Z','N','A',tabulation])
        df=df.sort_values(by=['Z','N'])
        df_interest=df.query(f'{Zmin_interest}<=Z<={Zmax_interest} and {Nmin_interest}<=N<={Nmax_interest}')
        Z_array=df_interest['Z'].to_numpy()
        A_array=df_interest['A'].to_numpy()
        external_array=df_interest[tabulation].to_numpy()
      
    delete_list=[]   
    
    for index1 in range(len(external_array)):
        
        for index2 in range(len(externalexp_array)):
            
            if (Z_array[index1]==Zexp_array[index2]) and (A_array[index1]==Aexp_array[index2]):
                
                delete_list.append(index1)
                
    Z_array=np.delete(Z_array,delete_list)
    A_array=np.delete(A_array,delete_list)
    external_array=np.delete(external_array,delete_list)                                      
    
    print('Start of gibbs minimization process...\n')
    start_time=time.time()    
    while True:
        
        Pfix=presslist[i]
        mini=100000.0
        flag=0.0
        ref=True
        
        for index in range(len(external_array)): 
            
            if Z_array[index]!=Z_array[index-1]:
            
                xroot=sp.optimize.fsolve(funtosolve, x0=100,
                                         args=(Z_array[index],Pfix,Bstar), xtol=1e-10) 
                
            flag=gibbs(A_array[index], Z_array[index], xroot, Bstar,
                       tabulation, external_array[index])
            
            if flag<mini:
                
                mini=flag
                Amin=A_array[index]
                Zmin=Z_array[index]
                Nmin=A_array[index]-Z_array[index]
                xemin=xroot
                externalmin=external_array[index]
                ref=False
                
        for index in range(len(externalexp_array)): 
            
            if Zexp_array[index]!=Zexp_array[index-1]:
            
                xroot=sp.optimize.fsolve(funtosolve, x0=100,
                                         args=(Zexp_array[index],Pfix,Bstar), xtol=1e-10) 
                
            flag=gibbs(Aexp_array[index], Zexp_array[index], xroot, Bstar,
                       'mass_excess', externalexp_array[index])
            
            if flag<mini:
                
                mini=flag
                Amin=Aexp_array[index]
                Zmin=Zexp_array[index]
                Nmin=Aexp_array[index]-Zexp_array[index]
                xemin=xroot
                externalmin=externalexp_array[index]
                ref=True
                                      
        if mini > Mn:
            break
        
        Aminlist[i]=Amin
        Zminlist[i]=Zmin
        Nminlist[i]=Nmin
        xeminlist[i]=xemin
        gminlist[i]=mini
        if ref==False:
            etotminlist[i]=etot(Amin, Zmin, xemin, Bstar, tabulation, externalmin)
        elif ref==True:
            etotminlist[i]=etot(Amin, Zmin, xemin, Bstar, 'mass_excess', externalmin)

        i=i+1
      
    end_time=time.time()
    print(f'Execution time for gibbs minimization: {(end_time-start_time):.2f} sec\n')
        
elif tabulation=='semf_be' :
    
    Z_array=np.zeros(((Zmax_interest-Zmin_interest+1)*(Nmax_interest-Nmin_interest+1)), dtype=np.int16)
    A_array=np.zeros(((Zmax_interest-Zmin_interest+1)*(Nmax_interest-Nmin_interest+1)), dtype=np.int16)
    external_array=np.zeros(((Zmax_interest-Zmin_interest+1)*(Nmax_interest-Nmin_interest+1)), dtype=np.float64)
    index=0
    
    for Z in range(Zmin_interest,Zmax_interest+1):
        
        for N in range(Nmin_interest,Nmax_interest+1):
            
            A=Z+N
            
            if (beSEMF(A, Z)<=0.0) or (Sp(A, Z)<=0.0) or (S2p(A, Z)<=0.0) or (Sn(A, Z)<=0.0) or (S2n(A, Z)<=0.0):
                continue
            
            Z_array[index]=Z
            A_array[index]=A
            external_array[index]=beSEMF(A, Z)
            index=index+1
    
    Z_array=Z_array[:index]
    A_array=A_array[:index]
    external_array=external_array[:index]
    
    delete_list=[]   
    
    for index1 in range(len(external_array)):
        
        for index2 in range(len(externalexp_array)):
            
            if (Z_array[index1]==Zexp_array[index2]) and (A_array[index1]==Aexp_array[index2]):
                
                delete_list.append(index1)
                
    Z_array=np.delete(Z_array,delete_list)
    A_array=np.delete(A_array,delete_list)
    external_array=np.delete(external_array,delete_list)  
     
    print('Start of gibbs minimization process...\n')
    start_time=time.time()     
    while True:
        
        Pfix=presslist[i]
        mini=100000.0
        flag=0.0
        ref=True
        
        for index in range(len(external_array)): 
            
            if Z_array[index]!=Z_array[index-1]:
            
                xroot=sp.optimize.fsolve(funtosolve, x0=100,
                                         args=(Z_array[index],Pfix,Bstar), xtol=1e-10) 
                
            flag=gibbs(A_array[index], Z_array[index], xroot, Bstar, tabulation)
            
            if flag<mini:
                
                mini=flag
                Amin=A_array[index]
                Zmin=Z_array[index]
                Nmin=A_array[index]-Z_array[index]
                xemin=xroot
                ref=False
                
        for index in range(len(externalexp_array)): 
            
            if Zexp_array[index]!=Zexp_array[index-1]:
            
                xroot=sp.optimize.fsolve(funtosolve, x0=100,
                                         args=(Zexp_array[index],Pfix,Bstar), xtol=1e-10) 
                
            flag=gibbs(Aexp_array[index], Zexp_array[index], xroot, Bstar,
                       'mass_excess', externalexp_array[index])
            
            if flag<mini:
                
                mini=flag
                Amin=Aexp_array[index]
                Zmin=Zexp_array[index]
                Nmin=Aexp_array[index]-Zexp_array[index]
                xemin=xroot
                externalmin=externalexp_array[index]
                ref=True
                      
        if mini > Mn:
            break
        
        Aminlist[i]=Amin
        Zminlist[i]=Zmin
        Nminlist[i]=Nmin
        xeminlist[i]=xemin
        gminlist[i]=mini
        if ref==False:
            etotminlist[i]=etot(Amin, Zmin, xemin, Bstar, tabulation)
        elif ref==True:
            etotminlist[i]=etot(Amin, Zmin, xemin, Bstar, 'mass_excess', externalmin)        

        i=i+1
      
    end_time=time.time()
    print(f'Execution time for gibbs minimization: {(end_time-start_time):.2f} sec\n')
                 
#%% end of gibbs minimization process    
presslist=presslist[:i]
Aminlist=Aminlist[:i]
Zminlist=Zminlist[:i]
Nminlist=Nminlist[:i]
xeminlist=xeminlist[:i]
gminlist=gminlist[:i]
etotminlist=etotminlist[:i]

mueminlist=mue(ge(xeminlist))
geminlist=ge(xeminlist)
nbminlist=nb(Aminlist, Zminlist, xeminlist, Bstar)
ZoverA=Zminlist/Aminlist
shearlist=shear(Zminlist, xeminlist, Bstar)
gamma=np.gradient(np.log(presslist), np.log(nbminlist)) 
compression=gamma*presslist
soundspeed=(np.gradient(presslist,etotminlist))**(1/2) 

changelist=[]

# find phase transition
for index in range(i-1):
    
    if (Zminlist[index]!=Zminlist[index+1]) or (Nminlist[index]!=Nminlist[index+1]):
        
        changelist.append(index)
        
# min - max values and A, Z, N, Z/A in each layer
pressminmax=np.zeros((len(changelist)+1, 2), dtype=np.float64)
nbminmax=np.zeros((len(changelist)+1, 2), dtype=np.float64)
etotminmax=np.zeros((len(changelist)+1, 2), dtype=np.float64)
mubminmax=np.zeros((len(changelist)+1, 2), dtype=np.float64)
mueminmax=np.zeros((len(changelist)+1, 2), dtype=np.float64)
Alayer=np.zeros((len(changelist)+1, 1), dtype=np.int16)
Zlayer=np.zeros((len(changelist)+1, 1), dtype=np.int16)
Nlayer=np.zeros((len(changelist)+1, 1), dtype=np.int16)
ZoverAlayer=np.zeros((len(changelist)+1, 1), dtype=np.float64)

for index in range(len(changelist)+1):
    if len(changelist)==0:
        pressminmax[index]=[presslist[index],presslist[-1]]
        nbminmax[index]=[nbminlist[index],nbminlist[-1]]
        etotminmax[index]=[etotminlist[index],etotminlist[-1]]
        mubminmax[index]=[gminlist[index],gminlist[-1]]
        mueminmax[index]=[mueminlist[index],mueminlist[-1]]
        Alayer[index]=Aminlist[index]
        Zlayer[index]=Zminlist[index]
        Nlayer[index]=Nminlist[index]
        ZoverAlayer[index]=ZoverA[index]
                
    elif index==0:
        pressminmax[index]=[presslist[index],presslist[changelist[index]]]
        nbminmax[index]=[nbminlist[index],nbminlist[changelist[index]]]
        etotminmax[index]=[etotminlist[index],etotminlist[changelist[index]]]
        mubminmax[index]=[gminlist[index],gminlist[changelist[index]]]
        mueminmax[index]=[mueminlist[index],mueminlist[changelist[index]]]
        Alayer[index]=Aminlist[index]
        Zlayer[index]=Zminlist[index]
        Nlayer[index]=Nminlist[index]
        ZoverAlayer[index]=ZoverA[index]
                
    elif index==len(changelist):
        pressminmax[index]=[presslist[changelist[index-1]+1],presslist[-1]]
        nbminmax[index]=[nbminlist[changelist[index-1]+1],nbminlist[-1]]
        etotminmax[index]=[etotminlist[changelist[index-1]+1],etotminlist[-1]]
        mubminmax[index]=[gminlist[changelist[index-1]+1],gminlist[-1]]
        mueminmax[index]=[mueminlist[changelist[index-1]+1],mueminlist[-1]]
        Alayer[index]=Aminlist[changelist[index-1]+1]
        Zlayer[index]=Zminlist[changelist[index-1]+1]
        Nlayer[index]=Nminlist[changelist[index-1]+1]
        ZoverAlayer[index]=ZoverA[changelist[index-1]+1]
        
    else:
        pressminmax[index]=[presslist[changelist[index-1]+1],
                         presslist[changelist[index]]]
        nbminmax[index]=[nbminlist[changelist[index-1]+1],
                         nbminlist[changelist[index]]]
        etotminmax[index]=[etotminlist[changelist[index-1]+1],
                         etotminlist[changelist[index]]]
        mubminmax[index]=[gminlist[changelist[index-1]+1],
                         gminlist[changelist[index]]]
        mueminmax[index]=[mueminlist[changelist[index-1]+1],
                         mueminlist[changelist[index]]]
        Alayer[index]=Aminlist[changelist[index-1]+1]
        Zlayer[index]=Zminlist[changelist[index-1]+1]
        Nlayer[index]=Nminlist[changelist[index-1]+1]
        ZoverAlayer[index]=ZoverA[changelist[index-1]+1]

#%% plots and interpolations of the above data
plt.figure()
plt.loglog(nbminlist, presslist, 'b-', linewidth=1) 
plt.grid()
plt.xlabel(r'Baryon number density [$fm^{-3}$]')
plt.ylabel(r'Pressure [$MeVfm^{-3}$]')
plt.title(rf'P vs $n_b$ (model: {model}, $B_\star={Bstar:.1f}$)')
plt.show()   
plt.savefig(f'Ptot_vs_nb_m{model}_f{field}.png')

plt.figure()
plt.loglog(etotminlist, presslist, 'b-', linewidth=1)
plt.grid()
plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
plt.ylabel(r'Pressure [$MeVfm^{-3}$]')
plt.title(rf'P vs $\epsilon$ (model: {model}, $B_\star={Bstar:.1f}$)')
plt.show()    
plt.savefig(f'Ptot_vs_etot_m{model}_f{field}.png')
   
plt.figure()
plt.loglog(etotminlist, shearlist, 'b-', linewidth=1)
plt.grid()
plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
plt.ylabel(r'Effective shear modulus [$MeVfm^{-3}$]')
plt.title(rf'S (shear modulus) vs $\epsilon$ (model: {model}, $B_\star={Bstar:.1f}$)')
plt.show()
plt.savefig(f'S_vs_etot_m{model}_f{field}.png')

plt.figure()
plt.loglog(nbminlist, shearlist, 'b-', linewidth=1)
plt.grid()
plt.xlabel(r'Baryon number density [$fm^{-3}$]')
plt.ylabel(r'Effective shear modulus [$MeVfm^{-3}$]')
plt.title(rf'S (shear modulus) vs $n_b$ (model: {model}, $B_\star={Bstar:.1f}$)')
plt.show()
plt.savefig(f'S_vs_nb_m{model}_f{field}.png')

plt.figure()
plt.semilogx(nbminlist, gamma, 'b-', linewidth=1)
if Bstar>1304:
    plt.ylim(1.5, 3.5)
plt.grid()
plt.xlabel(r'Baryon number density [$fm^{-3}$]')
plt.ylabel(r'Adiabatic index')
plt.title(rf'$\Gamma$ vs $n_b$ (model: {model}, $B_\star={Bstar:.1f}$)')
plt.axhline (y=4/3, color='red', linestyle='--', linewidth=1)
plt.show()  
plt.savefig(f'Gamma_vs_nb_m{model}_f{field}.png')

plt.figure()
plt.semilogx(etotminlist, gamma, 'b-', linewidth=1)
if Bstar>1304:
    plt.ylim(1.5, 3.5)
plt.grid()
plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
plt.ylabel(r'Adiabatic index')
plt.title(rf'$\Gamma$ vs $\epsilon$ (model: {model}, $B_\star={Bstar:.1f}$)')
plt.axhline (y=4/3, color='red', linestyle='--', linewidth=1)
plt.show()  
plt.savefig(f'Gamma_vs_etot_m{model}_f{field}.png')

plt.figure()
plt.semilogx(nbminlist, soundspeed, 'b-', linewidth=1)
plt.grid()
plt.xlabel(r'Baryon number density [$fm^{-3}$]')
plt.ylabel(r'Speed of sound [c]')
plt.title(rf'$c_s/c$ vs $n_b$ (model: {model}, $B_\star={Bstar:.1f}$)')
plt.show()  
plt.savefig(f'speed_vs_nb_m{model}_f{field}.png')

plt.figure()
plt.semilogx(etotminlist, soundspeed, 'b-', linewidth=1)
plt.grid()
plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
plt.ylabel(r'Speed of sound [c]')
plt.title(rf'$c_s/c$ vs $\epsilon$ (model: {model}, $B_\star={Bstar:.1f}$)')
plt.show()  
plt.savefig(f'speed_vs_etot_m{model}_f{field}.png')

plt.figure()
plt.loglog(etotminlist, compression, 'b-', linewidth=1)
plt.grid()
plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
plt.ylabel(r'Compression modulus [$MeVfm^{-3}$]')
plt.title(rf'$K$ (compression modulus) vs $\epsilon$ (model: {model}, $B_\star={Bstar:.1f}$)')
plt.show()
plt.savefig(f'compression_vs_etot_m{model}_f{field}.png')

plt.figure()
plt.loglog(nbminlist, compression, 'b-', linewidth=1)
plt.grid()
plt.xlabel(r'Baryon number density [$fm^{-3}$]')
plt.ylabel(r'Compression modulus [$MeVfm^{-3}$]')
plt.title(rf'$K$ (compression modulus) vs $n_b$ (model: {model}, $B_\star={Bstar:.1f}$)')
plt.show()
plt.savefig(f'compression_vs_nb_m{model}_f{field}.png')

plt.figure()
plt.semilogx(nbminlist, shearlist/compression, 'b-', linewidth=1)
plt.grid()
plt.xlabel(r'Baryon number density [$fm^{-3}$]')
plt.ylabel(r'$\frac{shear}{compression}$')
plt.title(rf'S/K vs $n_b$ (model: {model}, $B_\star={Bstar:.1f}$)')
plt.show()
plt.savefig(f'SK_vs_nb_m{model}_f{field}.png')

plt.figure()
plt.semilogx(etotminlist, shearlist/compression, 'b-', linewidth=1)
plt.grid()
plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
plt.ylabel(r'$\frac{shear}{compression}$')
plt.title(rf'S/K vs $\epsilon$ (model: {model}, $B_\star={Bstar:.1f}$)')
plt.show()
plt.savefig(f'SK_vs_etot_m{model}_f{field}.png')

plt.figure(figsize=(8,6))
plt.suptitle(rf'model: {model}, $B_\star={Bstar:.1f}$')
plt.subplot(3,1,1)
plt.semilogx(nbminlist, Zminlist, 'b-', linewidth=1)
plt.grid()
plt.ylabel(r'Z')
plt.subplot(3,1,2)
plt.semilogx(nbminlist, Nminlist, 'r-', linewidth=1)
plt.grid()
plt.ylabel(r'N')
plt.subplot(3,1,3)
plt.semilogx(nbminlist, ZoverA, 'g-', linewidth=1)
plt.grid()
plt.xlabel(r'Baryon number density [$fm^{-3}$]')
plt.ylabel(r'Z/A')
plt.show()  
plt.savefig(f'Z_N_ZoverA_vs_nb_m{model}_f{field}.png')

plt.figure(figsize=(8,6))
plt.suptitle(rf'model: {model}, $B_\star={Bstar:.1f}$')
plt.subplot(3,1,1)
plt.semilogx(etotminlist, Zminlist, 'b-', linewidth=1)
plt.grid()
plt.ylabel(r'Z')
plt.subplot(3,1,2)
plt.semilogx(etotminlist, Nminlist, 'r-', linewidth=1)
plt.grid()
plt.ylabel(r'N')
plt.subplot(3,1,3)
plt.semilogx(etotminlist, ZoverA, 'g-', linewidth=1)
plt.grid()
plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
plt.ylabel(r'Z/A')
plt.show()  
plt.savefig(f'Z_N_ZoverA_vs_etot_m{model}_f{field}.png')

# fitting Ptot=Ptot(etot)
def fitfun(x,K,G):
    
    ''' 
    a simple polytropic EoS for fitting
    
    '''
    return K*x**G

parameters=sp.optimize.curve_fit(fitfun, etotminlist, presslist, method='lm')[0]

plt.figure()
plt.loglog(etotminlist, presslist, 'bo', markersize=2, label='Data')
plt.loglog(etotminlist, fitfun(etotminlist, *parameters),
           'r-', linewidth=2, label=f'fitted curve\n K={parameters[0]:.4f}\n G={parameters[1]:.4f}')
plt.grid()
plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
plt.ylabel(r'Pressure [$MeVfm^{-3}$]')
plt.legend()
plt.title(rf'fitting curve: $P=K\epsilon^G$ (model: {model}, $B_\star={Bstar:.1f}$)')
plt.show()
plt.savefig(f'polytropicEOSfit_m{model}_f{field}.png')

# interpolation Ptot=Ptot(etot)
Ptot_interp=sp.interpolate.interp1d(etotminlist, presslist, kind='linear',
                                    bounds_error=False, fill_value=(presslist[0], presslist[-1]))
etotlistnew=np.logspace(np.log10(etotminlist[0]), np.log10(etotminlist[-1]),
                        25000, dtype=np.float64)
plt.figure()
plt.loglog(etotminlist, presslist, 'bo', markersize=2, label='Data')
plt.loglog(etotlistnew, Ptot_interp(etotlistnew), 'r-', linewidth=2, label='Interpolated function')
plt.grid()
plt.legend()
plt.xlabel(r'Energy density [$MeVfm^{-3}$]')
plt.ylabel(r'Pressure [$MeVfm^{-3}$]')
plt.title(rf'P vs $\epsilon$ (model: {model}, $B_\star={Bstar:.1f}$)')
plt.show()   
plt.savefig(f'Interpolated_function_Ptot_interp_m{model}_f{field}.png')

# interpolation etot=etot(Ptot)        
etot_interp=sp.interpolate.interp1d(presslist, etotminlist, kind='linear',
                                    bounds_error=False, fill_value=(etotminlist[0], etotminlist[-1]))
presslistnew=np.logspace(np.log10(presslist[0]), np.log10(presslist[-1]),
                         25000, dtype=np.float64)
plt.figure()
plt.loglog(presslist, etotminlist, 'bo', markersize=2, label='Data')
plt.loglog(presslistnew, etot_interp(presslistnew), 'r-', linewidth=2, label='Interpolated function')
plt.grid()
plt.legend()
plt.xlabel(r'Pressure [$MeVfm^{-3}$]')
plt.ylabel(r'Energy density [$MeVfm^{-3}$]')
plt.title(rf'$\epsilon$ vs P (model: {model}, $B_\star={Bstar:.1f}$)')
plt.show()   
plt.savefig(f'Interpolated_function_etot_interp_m{model}_f{field}.png') 

#%% general relativity - TOV equations
# convert units in powers of km (Geometrized unit system)
MSOLARgeom=1.4766 # solar mass in km
presslistgeom=MeV_fm3to_km2*presslist # pressure in km**(-2)
etotminlistgeom=MeV_fm3to_km2*etotminlist # energy density in km**(-2)
etotminmaxgeom=MeV_fm3to_km2*etotminmax

# interpolation etot=etot(Ptot), now in Geometrized unit system
etotgeom_interp=sp.interpolate.interp1d(presslistgeom, etotminlistgeom, kind='linear',
                                        bounds_error=False, fill_value=(etotminlistgeom[0], etotminlistgeom[-1]))

# TOV equations
def dm_dr(r,p):
    
    return 4*np.pi*r**2*etotgeom_interp(p)

def dp_dr(r,m,p):
    
    return -(((p+etotgeom_interp(p))*(m+4*np.pi*r**3*p))/(r*(r-2*m)))

# function for baryonic mass
def DMb(r, farg1, farg2):
    
  return 4*np.pi*r**2*farg1(r)*(1-2*farg2(r)/r)**(-1/2)

# function for proper depth
def Dz(r, farg):
    
    return (1-2*farg(r)/r)**(-1/2)

# total grav. mass M in solar mass unit, total radius R in km , h step in km
def tovresults(M, R, choice=False, h=-2e-6, tolerance=1e-15):
    
    '''
    RK4 method for the solution of TOV equations
    
    '''    
    rlist=np.zeros((int(R/(2*abs(h)))), dtype=np.float64)
    mlist=np.zeros((int(R/(2*abs(h)))), dtype=np.float64)
    plist=np.zeros((int(R/(2*abs(h)))), dtype=np.float64)
    
    rlist[0]=R
    mlist[0]=M*MSOLARgeom
    plist[0]=presslistgeom[0]
    
    for j in range(1,len(rlist)):
        
        k1=h*dm_dr(rlist[j-1], plist[j-1])
        k2=h*dp_dr(rlist[j-1], mlist[j-1], plist[j-1])
        s1=h*dm_dr(rlist[j-1]+h/2, plist[j-1]+k2/2)
        s2=h*dp_dr(rlist[j-1]+h/2, mlist[j-1]+k1/2, plist[j-1]+k2/2)
        l1=h*dm_dr(rlist[j-1]+h/2, plist[j-1]+s2/2)
        l2=h*dp_dr(rlist[j-1]+h/2, mlist[j-1]+s1/2, plist[j-1]+s2/2)
        p1=h*dm_dr(rlist[j-1]+h, plist[j-1]+l2)
        p2=h*dp_dr(rlist[j-1]+h, mlist[j-1]+l1, plist[j-1]+l2)
        
        mlist[j]=mlist[j-1]+(k1+2*s1+2*l1+p1)/6
        plist[j]=plist[j-1]+(k2+2*s2+2*l2+p2)/6
        rlist[j]=rlist[j-1]+h
        
        if presslistgeom[-1]-plist[j] < tolerance:
            break

    mlist=mlist[:j+1]
    plist=plist[:j+1]
    rlist=rlist[:j+1]
    
    etotlist=etotgeom_interp(plist)
    
    etotgeomfromr_interp=sp.interpolate.interp1d(rlist, etotlist, kind='linear',
                                                 bounds_error=False, fill_value=(etotlist[-1], etotlist[0]))
    mgeomfromr_interp=sp.interpolate.interp1d(rlist, mlist, kind='linear',
                                              bounds_error=False, fill_value=(mlist[-1], mlist[0]))
    rgeomfrometot_interp=sp.interpolate.interp1d(etotlist, rlist, kind='linear',
                                                   bounds_error=False, fill_value=(rlist[0], rlist[-1])) 
    
    Mcr=(mlist[0]-mlist[-1])/MSOLARgeom # mass of outer crust in solar mass unit
    Mcr_M=Mcr/M # fractional mass of outer crust
    DRcr=rlist[0]-rlist[-1] # coordinate thickness of outer crust in km
    DRcr_R=DRcr/R # fractional coordinate thickness of outer crust
      
    # baryonic mass of outer crust in solar mass unit
    Mbcr=(sp.integrate.quad(DMb, rlist[-1], rlist[0], limit=8000,
                             args=(etotgeomfromr_interp, mgeomfromr_interp))[0])/MSOLARgeom
    
    # relative fraction of outer crust's mass due to gravitational binding energy
    gravbe=(Mbcr-Mcr)/(Mcr)
    
    # proper depth of drip point in km
    zdrip=sp.integrate.quad(Dz, rlist[-1], rlist[0], limit=8000, args=(mgeomfromr_interp))[0]
    
    # total number of baryons in outer crust
    DN=(Mbcr)/(Mn*MeV_c2toMsolar)
    
    # now for each layer, defined by etotmin-etotmax 
    if choice==True:
        
        rminmax=np.zeros((len(changelist)+1, 2), dtype=np.float64) # in km
        DRlayer=np.zeros((len(changelist)+1, 1), dtype=np.float64) # in km
        DRlayer_DRcr=np.zeros((len(changelist)+1, 1), dtype=np.float64)
        zlayer=np.zeros((len(changelist)+1, 1), dtype=np.float64) # in km
        zlayer_zdrip=np.zeros((len(changelist)+1, 1), dtype=np.float64)
        Mblayer_Mbcr=np.zeros((len(changelist)+1, 1), dtype=np.float64)
        Mlayer_Mcr=np.zeros((len(changelist)+1, 1), dtype=np.float64)
        Mblayer=np.zeros((len(changelist)+1, 1), dtype=np.float64) # in solar mass unit
        Mlayer=np.zeros((len(changelist)+1, 1), dtype=np.float64) # in solar mass unit
        DNlayer=np.zeros((len(changelist)+1, 1), dtype=np.float64)
        DNlayer_DN=np.zeros((len(changelist)+1, 1), dtype=np.float64)
        
        for index in range(len(etotminmaxgeom)):
            
            rminmax[index]=[rgeomfrometot_interp(etotminmaxgeom[index][1]),
                            rgeomfrometot_interp(etotminmaxgeom[index][0])]
            DRlayer[index]=rminmax[index][1]-rminmax[index][0]
            DRlayer_DRcr[index]=DRlayer[index]/DRcr
            zlayer[index]=sp.integrate.quad(Dz, rminmax[index][0],
                                            rlist[0], limit=8000, args=(mgeomfromr_interp))[0]
            zlayer_zdrip[index]=zlayer[index]/zdrip
            Mblayer[index]=(sp.integrate.quad(DMb, rminmax[index][0],
                                              rminmax[index][1], limit=8000, args=(etotgeomfromr_interp, mgeomfromr_interp))[0])/MSOLARgeom
            Mblayer_Mbcr[index]=Mblayer[index]/Mbcr
            Mlayer[index]=(mgeomfromr_interp(rminmax[index][1])-mgeomfromr_interp(rminmax[index][0]))/MSOLARgeom
            Mlayer_Mcr[index]=Mlayer[index]/Mcr
            DNlayer[index]=(Mblayer[index])/(Mn*MeV_c2toMsolar)
            DNlayer_DN[index]=DNlayer[index]/DN
               
        return rminmax, DRlayer, DRlayer_DRcr, zlayer, zlayer_zdrip, Mblayer, Mblayer_Mbcr, Mlayer, Mlayer_Mcr, DNlayer, DNlayer_DN
      
    return Mcr, Mcr_M, DRcr, DRcr_R, Mbcr, gravbe, zdrip, DN
    
# unpack for choice==False, for a typical NS 
Mtyp, Rtyp = 1.4, 12.5
Mtyps, Rtyps =  (f'{Mtyp:.2f}').replace(".", "_"),  (f'{Rtyp:.2f}').replace(".", "_")

Mcr, Mcr_M, DRcr, DRcr_R, Mbcr, gravbe, zdrip, DN = tovresults(Mtyp, Rtyp)
print(f'nuclear model: {model}\nmagnetic field: {Bstar:.1f} (in critical field unit)\nfor M={Mtyp:.2f} (in solar mass unit) and R={Rtyp:.2f}km:\n\tMcr={Mcr:.3e}\n\tMcr_M={Mcr_M:.3e}\n\tDRcr={DRcr:.3e}\n\tDRcr_R={DRcr_R:.3e}\n\tMbcr={Mbcr:.3e}\n\tgravbe={gravbe:.3e}\n\tzdrip={zdrip:.3e}\n\tDN={DN:.3e}')
# unpack for choice==True, for a typical NS 
rminmax, DRlayer, DRlayer_DRcr, zlayer, zlayer_zdrip, Mblayer, Mblayer_Mbcr, Mlayer, Mlayer_Mcr, DNlayer, DNlayer_DN = tovresults(Mtyp, Rtyp, choice=True)

#%% plots from TOV equations
Mlist=np.array([1.2, 1.4, 1.6])
Rlist=np.array([10.0, 12.0, 14.0])

Mcr1_2=np.zeros(len(Rlist), dtype=np.float64)
Mcr1_4=np.zeros(len(Rlist), dtype=np.float64)
Mcr1_6=np.zeros(len(Rlist), dtype=np.float64)

zdrip1_2=np.zeros(len(Rlist), dtype=np.float64)
zdrip1_4=np.zeros(len(Rlist), dtype=np.float64)
zdrip1_6=np.zeros(len(Rlist), dtype=np.float64)

for index in range(len(Rlist)):
    
    result=tovresults(Mlist[0], Rlist[index])
    Mcr1_2[index]=result[0]
    zdrip1_2[index]=result[6]
    
    result=tovresults(Mlist[1], Rlist[index])
    Mcr1_4[index]=result[0]
    zdrip1_4[index]=result[6]
    
    result=tovresults(Mlist[2], Rlist[index])
    Mcr1_6[index]=result[0]
    zdrip1_6[index]=result[6]
    
del result

plt.figure()
plt.plot(Rlist, Mcr1_2, 'o-b', label=r'M=1.2$M_{\odot}$')
plt.plot(Rlist, Mcr1_4, 'o-r', label=r'M=1.4$M_{\odot}$')
plt.plot(Rlist, Mcr1_6, 'o-g', label=r'M=1.6$M_{\odot}$')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.grid()
plt.legend()
plt.xlabel(r'Total radius of NS [km]')
plt.ylabel(r'Grav. mass of outer crust [$M_{\odot}$]')
plt.title(rf'Mcr vs R (model: {model}, $B_\star={Bstar:.1f}$)')
plt.show()
plt.savefig(f'Mcr_m{model}_f{field}.png') 

plt.figure()
plt.plot(Rlist, zdrip1_2, 'o-b', label=r'M=1.2$M_{\odot}$')
plt.plot(Rlist, zdrip1_4, 'o-r', label=r'M=1.4$M_{\odot}$')
plt.plot(Rlist, zdrip1_6, 'o-g', label=r'M=1.6$M_{\odot}$')
plt.grid()
plt.legend()
plt.xlabel(r'Total radius of NS [km]')
plt.ylabel(r'Proper thickness of outer crust [km]')
plt.title(rf'$z_d$ vs R (model: {model}, $B_\star={Bstar:.1f}$)')
plt.show()
plt.savefig(f'zdrip_m{model}_f{field}.png') 

#%% save useful data arrays as txt files or as .npy files

def format_column(dataframe, column, fmt):
    
    '''
    return the desired float format for each column
    
    '''

    dataframe[column]=dataframe[column].astype(str).apply(lambda x: format(float(x), fmt))
    
    return dataframe

df1=pd.DataFrame(np.c_[presslist, etotminlist, nbminlist, Aminlist, Zminlist,
                       Nminlist, ZoverA, gminlist, mueminlist, xeminlist, geminlist],
                columns=['Ptot[MeVfm**(-3)]', 'etot[MeVfm**(-3)]',
                         'nb[fm**(-3)]', 'A', 'Z', 'N', 'Z/A', 'mu_b[MeV]',
                         'mu_e[MeV]', 'xe', 'ge'])

df1=format_column(df1,'Ptot[MeVfm**(-3)]','.3e')
df1=format_column(df1,'etot[MeVfm**(-3)]','.3e')
df1=format_column(df1,'nb[fm**(-3)]','.3e')
df1=format_column(df1,'A','.0f')
df1=format_column(df1,'Z','.0f')
df1=format_column(df1,'N','.0f')
df1=format_column(df1,'Z/A','.3f')
df1=format_column(df1,'mu_b[MeV]','.3f')
df1=format_column(df1,'mu_e[MeV]','.3f')
df1=format_column(df1,'xe','.3f')
df1=format_column(df1,'ge','.3f')

df1.to_csv(f'eos_m{model}_f{field}.dat', index=False,
           header=True, sep='\t', na_rep='N/A', mode='w')

df2=pd.DataFrame(np.hstack((pressminmax, etotminmax, nbminmax, Alayer, Zlayer,
                            Nlayer, ZoverAlayer, mubminmax, mueminmax)),
                 columns=['Pmin[MeVfm**(-3)]', 'Pmax[MeVfm**(-3)]',
                          'etotmin[MeVfm**(-3)]', 'etotmax[MeVfm**(-3)]',
                          'nbmin[fm**(-3)]', 'nbmax[fm**(-3)]',
                          'A', 'Z', 'N', 'Z/A', 'mu_bmin[MeV]', 'mu_bmax[MeV]',
                          'mu_emin[MeV]', 'mu_emax[MeV]'])

df2=format_column(df2,'Pmin[MeVfm**(-3)]','.3e')
df2=format_column(df2,'Pmax[MeVfm**(-3)]','.3e')
df2=format_column(df2,'etotmin[MeVfm**(-3)]','.3e')
df2=format_column(df2,'etotmax[MeVfm**(-3)]','.3e')
df2=format_column(df2,'nbmin[fm**(-3)]','.3e')
df2=format_column(df2,'nbmax[fm**(-3)]','.3e')
df2=format_column(df2,'A','.0f')
df2=format_column(df2,'Z','.0f')
df2=format_column(df2,'N','.0f')
df2=format_column(df2,'Z/A','.3f')
df2=format_column(df2,'mu_bmin[MeV]','.3f')
df2=format_column(df2,'mu_bmax[MeV]','.3f')
df2=format_column(df2,'mu_emin[MeV]','.3f')
df2=format_column(df2,'mu_emax[MeV]','.3f')

df2.to_csv(f'layers_m{model}_f{field}.dat', index=False,
           header=True, sep='\t', na_rep='N/A', mode='w')


df3=pd.DataFrame(np.hstack((rminmax, DRlayer, DRlayer_DRcr, zlayer, zlayer_zdrip, Mblayer, Mblayer_Mbcr, Mlayer, Mlayer_Mcr, DNlayer, DNlayer_DN)),
                 columns=['rmin[km]', 'rmax[km]', 'DRlayer[km]', 'DRlayer_DRcr', 
                          'zlayer[km]', 'zlayer_zdrip', 'Mblayer[solarmass]', 
                          'Mblayer_Mbcr', 'Mlayer[solarmass]', 'Mlayer_Mcr', 
                          'DNlayer', 'DNlayer_DN'])

df3=format_column(df3,'rmin[km]','.3f')
df3=format_column(df3,'rmax[km]','.3f')
df3=format_column(df3,'DRlayer[km]','.3e')
df3=format_column(df3,'DRlayer_DRcr','.3e')
df3=format_column(df3,'zlayer[km]','.3f')
df3=format_column(df3,'zlayer_zdrip','.3e')
df3=format_column(df3,'Mblayer[solarmass]','.3e')
df3=format_column(df3,'Mblayer_Mbcr','.3e')
df3=format_column(df3,'Mlayer[solarmass]','.3e')
df3=format_column(df3,'Mlayer_Mcr','.3e')
df3=format_column(df3,'DNlayer','.3e')
df3=format_column(df3,'DNlayer_DN','.3e')
                          
df3.to_csv(f'layers_m{model}_f{field}_M{Mtyps}_R{Rtyps}.dat', index=False,
           header=True, sep='\t', na_rep='N/A', mode='w')

np.save(f'nbminlist_m{model}_f{field}.npy', nbminlist)
np.save(f'presslist_m{model}_f{field}.npy', presslist)
np.save(f'etotminlist_m{model}_f{field}.npy', etotminlist)
np.save(f'shearlist_m{model}_f{field}.npy', shearlist)
np.save(f'gamma_m{model}_f{field}.npy', gamma)
np.save(f'soundspeed_m{model}_f{field}.npy', soundspeed)
np.save(f'compression_m{model}_f{field}.npy', compression)
np.save(f'Zminlist_m{model}_f{field}.npy', Zminlist)
np.save(f'Nminlist_m{model}_f{field}.npy', Nminlist)
np.save(f'ZoverA_m{model}_f{field}.npy', ZoverA)
np.save(f'Rlist_m{model}_f{field}.npy', Rlist)
np.save(f'Mcr1_2_m{model}_f{field}.npy', Mcr1_2)
np.save(f'Mcr1_4_m{model}_f{field}.npy', Mcr1_4)
np.save(f'Mcr1_6_m{model}_f{field}.npy', Mcr1_6)
np.save(f'zdrip1_2_m{model}_f{field}.npy', zdrip1_2)
np.save(f'zdrip1_4_m{model}_f{field}.npy', zdrip1_4)
np.save(f'zdrip1_6_m{model}_f{field}.npy', zdrip1_6)
