from matplotlib import pyplot as plt
import numpy as np
import camb
from camb.camb import csf_background
from camb.constants import Mpc, const_pi, G, m_p, eV, hbar



#### 宇宙学参数
H0 = 67.88                 # Hubble 率
h = H0/100
omk = 0                   # 空间曲率
ombh2 = 0.022             # 重子物质
omdmh2 = 0.122            # 暗物质

##### 暗物质
######## CSF 参数
m_phi = 1e-22            # 粒子质量，单位 eV
Num_density = 1e7        # 粒子数密度， 单位 cm^(-3)
omphih2 = 0.99*omdmh2      # 密度分数
dotR =  1e30                  # 径向速度，eV^2
a_ini = 1e-20


######## CDM
omch2 = omdmh2 - omphih2  # CDM

##### 有质量中微子
mnu = 0.06               # 中微子质量， eV


##### 再电离
tau = 0.06               # 再电离光深


##### 宇宙初始条件
As = 2e-9                # 标量功率谱振幅
ns = 0.965               # 标量谱指数

##### 非线性

halofit_version = 'mead'

##### CMB
lmax = 3000



#### CAMB 参数设置
pars = camb.set_params(H0 = H0, ombh2 = ombh2, omch2 = omch2, mnu = mnu, omk = omk, tau = tau,
                       As = As, ns = ns, halofit_version = halofit_version, lmax = lmax,
                       m_phi = m_phi, Num_density = Num_density, omphih2 = omphih2, dotR = dotR)

# print(pars)


csf_background(pars, get_phi=True)
CSFData = pars.DarkMatter




a, rho,rho_n, P, w, cs2 = CSFData.get_fulid()
a, phiR, phiI, dotphiR, dotphiI = CSFData.get_CSF()

H = CSFData.H_arr



from camb.constants import c
from scipy.integrate import quad, cumulative_trapezoid

def dtauda(ai):
    '''
    Lambda CDM 宇宙
    '''
    
    omb = ombh2/h**2 
    omc = omdmh2/h**2
    omr = 9e-5
    omde = 1- omb - omc - omr -omk
    H = H0*1e3/c*np.sqrt( omr/ai**4 + omb/ai**3 + omc/ai**3 + omk/ai**2 + omde )
    dtauda = 1/(ai**2*H)
    return dtauda


def tau(ai):
    tau, error = quad(dtauda,0,ai)
    return tau 



dtaudascf = 1/(a**2*H)


tauscf = cumulative_trapezoid(dtaudascf,a,initial=2.18638e-16)


N = len(a)
taus = np.zeros_like(a)

for ix in range(N):
    taus[ix] = tau(a[ix])


plt.plot(a, w)
plt.xscale('log')
plt.xlabel('a')
plt.ylabel('w_phi')
plt.show()



plt.plot(a,taus,label='lcdm')
plt.plot(a,tauscf,label='csf')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('a')
plt.ylabel(r'$\tau$')
plt.legend()
plt.show()



    




