from matplotlib import pyplot as plt
import numpy as np
import camb
from camb.camb import csf_background
from camb.constants import Mpc, const_pi, G, m_p, eV, hbar
import sympy as sp


#### 宇宙学参数
H0 = 67.88                 # Hubble 率
h = H0/100
omk = 0                   # 空间曲率
ombh2 = 0.022             # 重子物质
omdmh2 = 0.122            # 暗物质

##### 暗物质
######## CSF 参数
m_phi = 1e-26            # 粒子质量，单位 eV
Num_density = 1e-20        # 粒子数密度， 单位 cm^(-3)
omphih2 = 0.99*omdmh2      # 密度分数
dotR =  1e-20                  # 径向速度，eV^2
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

csf_background(pars, get_phi=True)

pars.WantScalars = True
pars.WantTensors = False

CSFData = pars.DarkMatter

a_osc = CSFData.a_osc
a_wkb = CSFData.a_wkb

data= camb.get_background(pars)

eta = 10**(np.linspace(-2, 4, 1000))


k = 5
if (k<0.01):
    ktxt = '$k = 10^{' + sp.latex(int(np.log10(k)))  +'}$ h/Mpc'
else:
    ktxt = ' k = ' + sp.latex(k) + 'h/Mpc'

    
ks =  k/h
ev1 = np.array(data.get_time_evolution(ks, eta))

a = ev1[:,17]
clxc = ev1[:,1]
clx_phi = ev1[:,13]

v_phi = ev1[:,14]

clx_n = ev1[:,15]

q = ev1[:,16]

print(clx_phi)

plt.plot(a,abs(clx_phi),ls='-',color='C0',label = ktxt)
plt.plot(a,abs(clxc),ls='--',color='C1',label= r'CDM')
plt.axvline(x=a_osc, color='k', linestyle='--', label=r'$a = a_{osc}$')
plt.axvline(x=a_wkb, color='k', linestyle='--', label=r'$a = a_{wkb}$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('a')
plt.ylabel('clx')
# plt.xlim(1e-7,1e0)
# plt.ylim(1e-18,1e3)
plt.legend()
plt.show()


plt.plot(a,abs(v_phi),ls='-',color='C0',label = ktxt)
plt.axvline(x=a_osc, color='k', linestyle='--', label=r'$a = a_{osc}$')
plt.axvline(x=a_wkb, color='k', linestyle='--', label=r'$a = a_{wkb}$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('a')
plt.ylabel('v')
# plt.xlim(1e-7,1e0)
# plt.ylim(1e-18,1e3)
plt.legend()
plt.show()

plt.plot(a,abs(clx_n),ls='-',color='C0',label = ktxt)
plt.axvline(x=a_osc, color='k', linestyle='--', label=r'$a = a_{osc}$')
plt.axvline(x=a_wkb, color='k', linestyle='--', label=r'$a = a_{wkb}$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('a')
plt.ylabel('clx_n')
# plt.xlim(1e-7,1e0)
# plt.ylim(1e-18,1e3)
plt.legend()
plt.show()


plt.plot(a,(q),ls='-',color='C0',label = ktxt)
plt.axvline(x=a_osc, color='k', linestyle='--', label=r'$a = a_{osc}$')
plt.axvline(x=a_wkb, color='k', linestyle='--', label=r'$a = a_{wkb}$')
plt.xscale('log')
# plt.yscale('log')
plt.xlabel('a')
plt.ylabel('q')
# plt.xlim(1e-7,1e0)
# plt.ylim(1e-18,1e3)
plt.legend()
plt.show()
