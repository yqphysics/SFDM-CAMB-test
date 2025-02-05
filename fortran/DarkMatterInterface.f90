module DarkMatterInterface
    use constants
    use precision
    use Interpolation
    use classes
    implicit none

    private

    type, extends(TCambComponent) :: TDarkMatterCSF
        !!!! 初始化参数
        logical  :: calculated_CSF = .false. 
        real(dl) :: m_phi = 0.0_dl             ! 标量场暗物质粒子的质量， 单位 eV
        real(dl) :: Num_density = 0.0_dl       ! 粒子数密度 （正粒子数密度 - 反粒子数密度， 单位 m^(-3)
        real(dl) :: omphih2 = 0.0_dl           ! 标量场的质量密度分数
        real(dl) :: dotR = 0.0_dl              ! CSF 的径向速度，单位 单位 kg^(1/2)·m^(-1/2)
        real(dl) :: a_ini = 10._dl**(-5._dl)   ! 开始演化 CSF 的初始尺度因子
        
        !!! 计算
        real(dl) :: a_osc = 1.0_dl             ! 开始振荡的尺度因子
        real(dl) :: a_wkb = 1.0_dl             ! 转换到 WKB 近似的尺度因子

        !!!!!! 精确部分和 WKB 近似部分
        real(dl), allocatable :: a_arr(:)       ! 完整的宇宙尺度因子演化
        real(dl), allocatable :: H_arr(:)       ! 完整的宇宙学 Hubble 参数演化
        real(dl), allocatable :: rho_phi(:)     ! CSF 能量密度
        real(dl), allocatable :: p_phi(:)       ! CSF 压强
        real(dl), allocatable :: w_phi(:)       ! CSF 状态方程
        real(dl), allocatable :: cs2(:)         ! CSF 绝热声速
        real(dl), allocatable :: phi1_arr(:)
        real(dl), allocatable :: phi2_arr(:)
        real(dl), allocatable :: dotphi1_arr(:)
        real(dl), allocatable :: dotphi2_arr(:)

        !!!!!!!! a = a_ini 的 CSF 的场量值
        real(dl) :: R0 , dotR0, theta0, dottheta0

    contains
    procedure, nopass :: PythonClass => TDarkMatterCSF_PythonClass
    procedure, nopass :: SelfPointer => TDarkMatterCSF_SelfPointer
    procedure :: grho_phia2
    procedure :: State_equation
    
    end type TDarkMatterCSF

    public TDarkMatterCSF

    
contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! TDarkMatterCSF

    function TDarkMatterCSF_PythonClass()
        character(LEN=:), allocatable :: TDarkMatterCSF_PythonClass
        
        TDarkMatterCSF_PythonClass = 'DarkMatterCSF'

    end function TDarkMatterCSF_PythonClass

    subroutine TDarkMatterCSF_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TDarkMatterCSF), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TDarkMatterCSF_SelfPointer

    function grho_phia2(this, a)
        class(TDarkMatterCSF) :: this
        Type(TCubicSpline) :: Irreg
        real(dl) :: a, grho_phia2, rho1

        if (a< this%a_arr(2)) then
            rho1 = this%rho_phi(1)

            grho_phia2 = kappa * rho1 * a**4._dl/c**2

        else

            if (a< this%a_wkb) then 

                call Irreg%Init(this%a_arr, this%rho_phi)
                grho_phia2 = kappa * Irreg%Value(a) * a**4._dl/c**2
                
            else
                grho_phia2 = 3._dl *(1000._dl)**2._dl * 100._dl**2 * this%omphih2*a/c**2

            endif
        
        endif
        
        
    
        
    end function grho_phia2


    subroutine State_equation(this, a, w_a, cs2_a )
        class(TDarkMatterCSF) :: this
        Type(TCubicSpline) :: Irreg1
        Type(TCubicSpline) :: Irreg2
        real(dl) :: a, w_a, cs2_a, w0 


        if (a< this%a_arr(2)) then
            w_a = this%w_phi(1)

            cs2_a = w_a

        else

            if (a< this%a_wkb) then 

                call Irreg1%Init(this%a_arr, this%w_phi)
                call Irreg2%Init(this%a_arr, this%cs2)

                w_a = Irreg1%Value(a)
                cs2_a = Irreg2%Value(a)
                
            else
                w_a = 0._dl
                cs2_a = 0._dl
            endif
        
        endif
    
        
    end subroutine State_equation

    
    
end module DarkMatterInterface