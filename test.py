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
m_phi = 1e-27             # 粒子质量，单位 eV
Num_density = 10**(21)          # 粒子数密度， 单位 m^(-3)
omphih2 = 0.056             # 密度分数
dotR = 10**(22.20033)*(eV)**(1/2)               # CSF 的径向速度，单位 eV^(1/2)
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

# print(pars)


csf_background(pars, get_phi=True)

# print(pars)

CSFData = pars.DarkMatter

# results = camb.get_results(pars)



# print(results)



# age = camb.get_age(pars)

# print(age)


# # print(CSFData)

a_osc = CSFData.a_osc
# a_wkb = CSFData.a_wkb
# print('a_wkb =',a_wkb)
# print('a_osc =',a_osc)

a_arr = np.array(CSFData.a_arr)
H_arr = np.array(CSFData.H_arr)
rho_phi = np.array(CSFData.rho_phi)
p_phi = np.array(CSFData.p_phi)
w_phi = np.array(CSFData.w_phi)
cs2 = np.array(CSFData.cs2)
phi1 = np.array(CSFData.phi1_arr)
phi2 = np.array(CSFData.phi2_arr)
dotphi1 = np.array(CSFData.dotphi1_arr)
dotphi2 = np.array(CSFData.dotphi2_arr)
R = np.sqrt((phi1)**2 +(phi2)**2)


###### 绘图
# N = len(H_arr)
# print(a_arr[N-1])

# plt.plot(a_arr, H_arr)
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('a')
# plt.ylabel('H')
# plt.show()

plt.plot(a_arr, rho_phi)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('a')
plt.ylabel('rho_phi')
plt.show()

plt.plot(a_arr, p_phi)
plt.xscale('log')
plt.xlabel('a')
plt.ylabel('p_phi')
plt.show()

plt.plot(a_arr, w_phi)
plt.xscale('log')
plt.xlabel('a')
plt.ylabel('w_phi')
plt.show()

# plt.plot(a_arr, cs2)
# plt.xscale('log')
# plt.xlabel('a')
# plt.ylabel('cs2')
# plt.show()

# plt.plot(a_arr, phi1, linewidth = 2.5)
# plt.xscale("log")
# # plt.yticks([])
# plt.tick_params(axis="both",labelsize = 12)
# plt.gca().get_xticklabels(prop)
# plt.xlabel(r'$a$', fontproperties= prop ,fontsize = 15)
# plt.ylabel(r'$\phi_R$', fontproperties= prop,fontsize = 15)
# plt.show()

# plt.plot(a_arr, phi2)
# plt.xscale('log')
# plt.xlabel('a')
# plt.ylabel('phi2')
# plt.show()

# plt.plot(a_arr, R)
# plt.xscale('log')
# plt.xlabel('a')
# plt.ylabel('R')
# plt.show()


# plt.plot(a_arr, dotphi1)
# plt.xscale('log')
# plt.xlabel('a')
# plt.ylabel('dotphi1')
# plt.show()

# plt.plot(a_arr, dotphi2)
# plt.xscale('log')
# plt.xlabel('a')
# plt.ylabel('dotphi2')
# plt.show()