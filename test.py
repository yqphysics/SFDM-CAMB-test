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
m_phi = 1e-26            # 粒子质量，单位 eV
Num_density = 1e2        # 粒子数密度， 单位 cm^(-3)
omphih2 = 0.99*omdmh2      # 密度分数
dotR =  1e2                  # 径向速度，eV^2
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

print(pars)


csf_background(pars, get_phi=True)

# print(pars.DarkMatter)

CSFData = pars.DarkMatter

# print(CSFData)

a_osc = CSFData.a_osc
a_wkb = CSFData.a_wkb


a = np.array(CSFData.a_arr)
grho_phi = np.array(CSFData.grho_phi)
dotgrho_phi = np.array(CSFData.dotgrho_phi)


# plt.plot(a, grho_phi*a**2)
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('a')
# plt.ylabel('grho_phi')
# plt.show()




age = camb.get_age(pars)

print(age)



a, rho,rho_n, P, w, cs2 = CSFData.get_fulid()
# a, phiR, phiI, dotphiR, dotphiI = CSFData.get_CSF()

###### 绘图

# plt.plot(a, H)
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('a')
# plt.ylabel('H')
# plt.show()

# plt.plot(a, rho*a**3)
# plt.plot(a, rho_n*a**3)
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('a')
# plt.ylabel('rho')
# plt.show()

# plt.plot(a, P)
# plt.xscale('log')
# plt.xlabel('a')
# plt.ylabel('P')
# plt.show()

plt.plot(a, w)
plt.axvline(x=a_osc, color='k', linestyle='--', label=r'$a = a_{osc}$')
plt.axvline(x=a_wkb, color='k', linestyle='--', label=r'$a = a_{wkb}$')
plt.xscale('log')
plt.xlabel('a')
plt.ylabel('w_phi')
plt.show()



# plt.plot(a, cs2)
# plt.xscale('log')
# plt.xlabel('a')
# plt.ylabel('cs2')
# plt.show()




# plt.plot(a, phiR)
# plt.xscale('log')
# plt.xlabel('a')
# plt.ylabel('phiR')
# plt.show()




# plt.plot(a, phiI)
# plt.xscale('log')
# plt.xlabel('a')
# plt.ylabel('phiI')
# plt.show()

# plt.plot( phiR,  phiI)
# # plt.xscale('log')
# plt.xlabel('phiR')
# plt.ylabel('phiI')
# plt.show()

# dotgrho = CSFData.dotgpres_phi



# plt.plot( a,  dotgrho)
# plt.xscale('log')
# plt.xlabel('a')
# plt.ylabel('dotgrho')
# plt.show()



# data = np.loadtxt('data.dat')


# ain = data[:, 0]  # 第一列
# u1 = data[:, 1]  # 第二列
# u2 = data[:, 2]  # 第三列
# u3 = data[:, 3]  # 第四列
# u4 = data[:, 4]  # 第五列

# y1 = 1e-5*ain
# y2 = 2.65e3*(ain/1e-20)**(0.45)
# # y3 = ain**(-2)

# plt.plot( ain,  abs(u1),label='fit')
# # plt.plot( ain, y1 )
# # plt.plot( ain, y2 )
# plt.axvline(x=a_osc, color='k', linestyle='--', label=r'$a = a_{osc}$')
# plt.axvline(x=a_wkb, color='k', linestyle='--', label=r'$a = a_{wkb}$')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('a')
# plt.ylabel('rho')
# plt.show()

# plt.plot( ain,  abs(u2),label='fit')
# # plt.plot( ain, y1 )
# # plt.plot( ain, y2 )
# plt.xscale('log')
# # plt.yscale('log')
# plt.xlabel('a')
# plt.ylabel('u2')
# plt.show()


# plt.plot( ain,  abs(u3),label='fit')
# # plt.plot( ain, y2 )
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('a')
# plt.ylabel('u3')
# plt.show()

# plt.plot( ain, abs(u4),label='fit')
# plt.plot( ain, y2 )
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('a')
# plt.ylabel('u4')
# plt.show()