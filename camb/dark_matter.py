from .baseconfig import F2003Class, fortran_class, numpy_1d, CAMBError, np, \
    AllocatableArrayDouble, f_pointer, camblib
from ctypes import c_int, c_double, byref, POINTER, c_bool
from .constants import eV, hbar


 
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
        ("Num_density",c_double,"粒子数密度 （正粒子数密度 - 反粒子数密度）"),
        ("omphih2",c_double,"标量场的质量密度分数"),
        ("dotR",c_double,"CSF 的径向速度"),
        ("a_ini",c_double,"开始演化 CSF 的初始尺度因子"),
        ("a_osc",c_double,"开始振荡的尺度因子"),
        ("a_wkb",c_double,"转换计算 WKB 近似的尺度因子"),
        ("a_arr",AllocatableArrayDouble,"完整的宇宙尺度因子演化"),
        ("H_arr",AllocatableArrayDouble,"完整的宇宙学 Hubble 参数演化"),
        ("rho_phi",AllocatableArrayDouble,"CSF 能量密度"),
        ("p_phi",AllocatableArrayDouble,"CSF 压强"),
        ("w_phi",AllocatableArrayDouble,"CSF 状态方程"),
        ("cs2",AllocatableArrayDouble,"CSF 绝热声速"),
        ("phi1_arr",AllocatableArrayDouble,"phi1_arr"),
        ("phi2_arr",AllocatableArrayDouble,"phi2_arr"),
        ("dotphi1_arr",AllocatableArrayDouble,"dotphi1_arr"),
        ("dotphi2_arr",AllocatableArrayDouble,"dotphi2_arr"),
        ("R0",c_double,"初始 R0"),
        ("dotR0",c_double,"初始 dotR0"),
        ("theta0",c_double,"初始 theta0"),
        ("dottheta0",c_double,"初始 dottheta0"),
    ]


    def set_params(self,m_phi, Num_density, omphih2, dotR, a_ini):
        self.m_phi = m_phi
        self.Num_density = Num_density
        self.omphih2 = omphih2
        self.dotR = dotR
        self.a_ini = a_ini





    

