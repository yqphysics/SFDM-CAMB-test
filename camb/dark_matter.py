from .baseconfig import F2003Class, fortran_class, numpy_1d, CAMBError, np, \
    AllocatableArrayDouble, f_pointer, camblib
from ctypes import c_int, c_double, byref, POINTER, c_bool
from .constants import eV, hbar, c, Mpc, kappa
TeV = 1e12*eV   # TeV
GeV = 1e9*eV
J_TeV = 1/TeV
J_GeV = 1/GeV
meter_TeV = 5.067730716156394e+18
meter_GeV = 5067730716156394.0

int_arg = POINTER(c_int)
d_arg = POINTER(c_double)
 
@fortran_class
class DarkMatterCSF(F2003Class):
    '''
    Abstract base class for dark matter model implementations.
    '''
    _fortran_class_module_ = 'DarkMatterInterface'
    _fortran_class_name_ = 'TDarkMatterCSF'


    _fields_ = [
        ("calculated_CSF",c_bool,"判断是否完成背景complex SFDM的计算"),
        ("m_phi", c_double,"标量场暗物质粒子的质量"),
        ("Num_density",c_double,"共动粒子数密度 （正粒子数密度 - 反粒子数密度）"),
        ("omphih2",c_double,"标量场的质量密度分数"),
        ("dotR",c_double,"CSF 的径向速度"),
        ("a_ini",c_double,"开始演化 CSF 的初始尺度因子"),
        ("mtoH_osc",c_double,"开始振荡时的 m/H "),
        ("mtoH_wkb",c_double,"切换WKB近似时的 m/H"),
        ("N_ex",c_int,"精确解数组的大小"),
        ("N_tot",c_int,"数组大小"),
        ("mnorm",c_double,"无量纲质量"),
        ("Nnorm",c_double,"无量纲共动粒子数密度"),
        ("dotR_norm",c_double,"无量纲径向速度"),
        ("a_osc",c_double,"开始振荡的尺度因子"),
        ("a_wkb",c_double,"转换计算 WKB 近似的尺度因子"),
        ("a_arr",AllocatableArrayDouble,"完整的宇宙尺度因子演化"),
        ("H_arr",AllocatableArrayDouble,"完整的宇宙学 Hubble 参数演化"),
        ("grho_phi",AllocatableArrayDouble,"无量纲 CSF 能量密度"),
        ("dotgrho_phi",AllocatableArrayDouble,"无量纲 CSF 能量密度时间导数"),
        ("gpres_phi",AllocatableArrayDouble,"无量纲 CSF 压强"),
        ("dotgpres_phi",AllocatableArrayDouble,"无量纲 CSF 压强时间导数"),
        ("grhon_phi",AllocatableArrayDouble,"无量纲 CSF 粒子数密度"),
        ("w_phi",AllocatableArrayDouble,"CSF 状态方程"),
        ("cs2",AllocatableArrayDouble,"CSF 绝热声速"),
        ("phi1",AllocatableArrayDouble,"phi1_arr"),
        ("phi2",AllocatableArrayDouble,"phi2_arr"),
        ("dotphi1",AllocatableArrayDouble,"dotphi1_arr"),
        ("dotphi2",AllocatableArrayDouble,"dotphi2_arr"),
    ]

    # _methods_ = [("grho_phia2",[d_arg],c_double)
    #              ]


    def set_params(self,m_phi, Num_density, omphih2, dotR,*keys):
        self.m_phi = m_phi
        self.Num_density = Num_density
        self.omphih2 = omphih2
        self.dotR = dotR
        if not keys:
            self.mtoH_osc= 3
            self.mtoH_wkb = 100
            self.N_tot= 1000
            self.a_ini = 1e-20



    def get_fulid(self,*keys):
        '''
        将计算结果转换为国际单位制
        '''
        a = np.array(self.a_arr)
        rho = np.array(self.grho_phi)/a**2*(c**4/(kappa*Mpc**2))*J_GeV/(meter_GeV**3)
        rho_n = np.array(self.grhon_phi)/a**2*(c**4/(kappa*Mpc**2))*J_GeV/(meter_GeV**3)
        P = np.array(self.gpres_phi)/a**2*(c**4/(kappa*Mpc**2))
        w = np.array(self.w_phi)
        cs2 = np.array(self.cs2)

        ### 列表
        txt = ['rho', 'rho_n','P','w','cs2']
        
        all = [a, rho, rho_n, P, w, cs2]

        if not keys:
            results = all
        else:
            results = [a]
            for ix in range(len(keys)):
                index = txt.index(keys[ix])
                results.append(all[index+1])
        return results


    def get_CSF(self,*keys):
        '''
        获得场的变化
        '''
        a = np.array(self.a_arr)
        phiR = np.array(self.phi1)*(c**4/kappa)**(0.5)
        phiI = np.array(self.phi2)*(c**4/kappa)**(0.5)
        dotphiR = np.array(self.dotphi1)/a*(c**4/kappa)**(0.5)
        dotphiI = np.array(self.dotphi2)/a*(c**4/kappa)**(0.5)

        txt = [ 'phiR', 'phiI', 'dotphiR', 'dotphiI']
        all = [a, phiR, phiI, dotphiR, dotphiI]

        if not keys:
            results = all
        else:
            results = [a]
            for ix in range(len(keys)):
                index = txt.index(keys[ix])
                results.append(all[index+1])
        return results
        

