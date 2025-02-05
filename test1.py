from matplotlib import pyplot as plt
import numpy as np
import camb
from camb.camb import csf_background
from camb.constants import Mpc, const_pi, G, m_p, eV, hbar


#### 宇宙学参数
H0 = 67.5                 # Hubble 率
omk = 0                   # 空间曲率
ombh2 = 0.022             # 重子物质
omdmh2 = 0.122            # 暗物质

##### 暗物质
######## CSF 参数
m_phi = 1e-22             # 粒子质量，单位 eV
Num_density = 10**(-21)          # 粒子数密度， 单位 m^(-3)
omphih2 = 0.056             # 密度分数
dotR = 10**(-22.20033)*(eV)**(1/2)               # CSF 的径向速度，单位 eV^(1/2)
a_ini = 1e-10


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
                       m_phi = m_phi, Num_density = Num_density, omphih2 = omphih2, dotR = dotR, a_ini = a_ini)

csf_background(pars, get_phi=True)

pars.WantScalars = True
pars.WantTensors = False

CSFData = pars.DarkMatter

a_osc = CSFData.a_osc
a_wkb = CSFData.a_wkb

data= camb.get_background(pars)

eta = 10**(np.linspace(0, 4, 500))
y = 10**(-2)*eta**2
ks = 0.2
ev = data.get_time_evolution(ks, eta)
ev = np.array(ev)

x2 = np.array([a_wkb,a_wkb])
y2 = np.array([1e-40,-3e-40])

plt.plot(eta,abs(ev[:,1]))
plt.plot(eta,abs(ev[:,4]))
# plt.plot(eta,y)
plt.xscale('log')
plt.yscale('log')
plt.show()


plt.plot(eta,abs(ev[:,3]))
plt.plot(eta,abs(ev[:,7]))
# plt.plot(eta,y)
plt.xscale('log')
plt.yscale('log')
plt.show()

np.savetxt('ev.data', ev)