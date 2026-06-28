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
        real(dl) :: Num_density = 0.0_dl       ! 共动粒子数密度 （正粒子数密度 - 反粒子数密度， 单位 cm^(-3)
        real(dl) :: omphih2 = 0.0_dl           ! 标量场的质量密度分数
        real(dl) :: dotR_ini = 0.0_dl          ! CSF 的径向速度
        real(dl) :: a_ini = 10._dl**(-5._dl)   ! 开始演化 CSF 的初始尺度因子
        real(dl) :: mtoH_osc = 10._dl          ! 开始振荡时的 m/H 
        real(dl) :: mtoH_wkb = 1000._dl        ! 切换WKB近似时的 m/H
        integer  :: N_ex                       ! 精确解数组的大小
        integer  :: N_tot                      ! 数组大小
        
        !!! 计算
        real(dl) :: mnorm                      ! 无量纲质量 mnorm = ma*eV/(c*hbar)*Mpc
        real(dl) :: Nnorm                      ! 无量纲共动粒子数密度
        real(dl) :: dotR_norm                  ! 无量纲径向速度
        real(dl) :: a_osc = 1.0_dl             ! 开始振荡的尺度因子
        real(dl) :: a_wkb = 1.0_dl             ! 转换到 WKB 近似的尺度因子

        !!!!!! 精确部分和 WKB 近似部分
        real(dl), allocatable :: a_arr(:)       ! 完整的宇宙尺度因子演化
        real(dl), allocatable :: H_arr(:)       ! 完整的宇宙学 Hubble 参数演化
        real(dl), allocatable :: grho_phi(:)    ! 无量纲 CSF 能量密度
        real(dl), allocatable :: dotgrho_phi(:) ! 无量纲 CSF 能量密度时间导数
        real(dl), allocatable :: gpres_phi(:)   ! 无量纲 CSF 压强
        real(dl), allocatable :: dotgpres_phi(:)! 无量纲 CSF 压强时间导数
        real(dl), allocatable :: grhon_phi(:)   ! 无量纲 CSF 粒子数密度
        real(dl), allocatable :: w_phi(:)       ! CSF 状态方程
        real(dl), allocatable :: cs2(:)         ! CSF 绝热声速
        real(dl), allocatable :: phi1(:)
        real(dl), allocatable :: phi2(:)
        real(dl), allocatable :: dotphi1(:)
        real(dl), allocatable :: dotphi2(:)

    contains
    procedure, nopass :: PythonClass => TDarkMatterCSF_PythonClass
    procedure, nopass :: SelfPointer => TDarkMatterCSF_SelfPointer
    procedure :: grho_phia2
    procedure :: dotgrho
    procedure :: gpres_phia2
    procedure :: dotgpres
    procedure :: get_uphi
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
        !!!! 能量密度插值 grhoa2 = 8 pi G *rho * a^4
        !!!! 早期宇宙采用 w=1和 w=-1 的标量场外延方法
        class(TDarkMatterCSF) :: this
        Type(TCubicSpline) :: Irreg
        real(dl) :: a, grho_phia2

        if (a< this%a_arr(1)) then
            if (this%w_phi(1)==1._dl) then 
                grho_phia2 = this%grho_phi(1)*this%a_arr(1)**2._dl*(a/this%a_arr(1))**(-2._dl)  ! w = 1 的相
            else
                grho_phia2 = this%grho_phi(1)*this%a_arr(1)**2._dl*(a/this%a_arr(1))**(4._dl)   ! w = -1 的相
            endif

        else
            if (a<this%a_wkb) then
                call Irreg%Init(this%a_arr(1:this%N_ex), this%grho_phi(1:this%N_ex))
                grho_phia2 = Irreg%Value(a)*a**2._dl
            else
                grho_phia2 = 3._dl *(1000._dl)**2._dl * 100._dl**2 * this%omphih2*a/c**2
            endif
        endif
        
        
    
        
    end function grho_phia2


    function dotgrho(this, a)
        !!! 插值能量密度时间导数 dotgrho = 8 pi G *rho' * a^2
        class(TDarkMatterCSF) :: this
        Type(TCubicSpline) :: Irreg
        real(dl) :: a, dotgrho

        if (a< this%a_arr(1)) then
            if (this%w_phi(1)==1._dl) then 
                dotgrho = this%dotgrho_phi(1)*(a/this%a_arr(1))**(-5._dl)
            else
                dotgrho = 0._dl
            endif
        else
            call Irreg%Init(this%a_arr(1:this%N_tot), this%dotgrho_phi(1:this%N_tot))
            dotgrho = Irreg%Value(a)

        endif

        
        if ( a<this%a_osc .and. abs(dotgrho)<1e-50_dl) then
            dotgrho = 0._dl
        endif
        
    
        
    end function dotgrho


    function gpres_phia2(this, a)
        class(TDarkMatterCSF) :: this
        Type(TCubicSpline) :: Irreg
        real(dl) :: a, gpres_phia2

        if (a< this%a_arr(1)) then
            if (this%w_phi(1)==1._dl) then 
                gpres_phia2 = this%gpres_phi(1)*this%a_arr(1)**2._dl*(a/this%a_arr(1))**(-2._dl)
            else
                gpres_phia2 = this%gpres_phi(1)*this%a_arr(1)**2._dl*(a/this%a_arr(1))**(4._dl)
            endif

        else
            if (a<this%a_wkb) then
                call Irreg%Init(this%a_arr(1:this%N_ex), this%gpres_phi(1:this%N_ex))
                gpres_phia2 = Irreg%Value(a)*a**2._dl
            else
                gpres_phia2 = 0.0_dl
            endif
        endif
        
    
        
    end function gpres_phia2


    function dotgpres(this, a)
        class(TDarkMatterCSF) :: this
        Type(TCubicSpline) :: Irreg
        real(dl) :: a, dotgpres

        if (a< this%a_arr(1)) then
            if (this%w_phi(1)==1._dl) then 
                dotgpres = this%dotgpres_phi(1)*(a/this%a_arr(1))**(-5._dl)
            else
                dotgpres = this%dotgpres_phi(1)*(a/this%a_arr(1))**(-5._dl)
            endif

        else
            call Irreg%Init(this%a_arr(1:this%N_tot), this%dotgpres_phi(1:this%N_tot))
            dotgpres = Irreg%Value(a)

        endif

        if ( a<this%a_osc .and. abs(dotgpres)<1e-60_dl) then
            dotgpres = 0._dl
        endif
    
        
    end function dotgpres


    subroutine get_uphi(this,  a, uphi)
        !!!! 需要修改，应该严格解析求解
        class(TDarkMatterCSF) :: this
        Type(TCubicSpline) :: Irreg1
        Type(TCubicSpline) :: Irreg2
        Type(TCubicSpline) :: Irreg3
        Type(TCubicSpline) :: Irreg4
        real(dl) :: a, uphi(4)

        if (a< this%a_arr(1)) then
            if (this%w_phi(1)==1._dl) then 
                uphi(1) = this%phi1(1)
                uphi(2) = this%phi2(1)
                uphi(3) = this%dotphi1(1)*(a/this%a_arr(1))**(0.45_dl)
                uphi(4) = this%dotphi2(1)*(a/this%a_arr(1))**(0.45_dl)
            else
                uphi(1) = this%phi1(1)
                uphi(2) = this%phi2(1)
                uphi(3) = this%dotphi1(1)*(a/this%a_arr(1))**(0.45_dl)
                uphi(4) = this%dotphi2(1)*(a/this%a_arr(1))**(-2._dl)
            endif
        else
            call Irreg1%Init(this%a_arr,this%phi1)
            call Irreg2%Init(this%a_arr,this%phi2)
            call Irreg3%Init(this%a_arr,this%dotphi1)
            call Irreg4%Init(this%a_arr,this%dotphi2)
            uphi(1) = Irreg1%Value(a)
            uphi(2) = Irreg2%Value(a)
            uphi(3) = Irreg3%Value(a)
            uphi(4) = Irreg4%Value(a)

        endif

    
        
    end subroutine get_uphi




    subroutine State_equation(this, a, w_a, cs2 )
        class(TDarkMatterCSF) :: this
        Type(TCubicSpline) :: Irreg1
        Type(TCubicSpline) :: Irreg2
        real(dl) :: a, w_a, cs2, w0 



        if (a< this%a_arr(1)) then
            w_a = this%w_phi(1)
        else

            if (a< this%a_wkb) then 

                call Irreg1%Init(this%a_arr(1:this%N_ex), this%w_phi(1:this%N_ex))
                w_a = Irreg1%Value(a)
            else
                w_a = 0._dl
            endif
        endif

        if (a< this%a_arr(1)) then
            cs2 = this%cs2(1)
        else

            if (a< this%a_wkb) then 

                call Irreg2%Init(this%a_arr(1:this%N_ex), this%cs2(1:this%N_ex))
                cs2 = Irreg2%Value(a)
            else
                cs2 = 0._dl
            endif
        endif
    
    
        
    end subroutine State_equation
    

    
    
end module DarkMatterInterface