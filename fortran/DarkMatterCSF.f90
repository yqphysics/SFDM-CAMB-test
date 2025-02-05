!!! Equations module allowing for the complex scalar field models
!!! Authors :: Qi Yang(杨祺)，  BoHua Li(李博华)


module DarkMatterCSF
    use MiscUtils
    use RangeUtils
    use constants
    use results
    use RKF45
    use MathUtils
    use DarkMatterInterface
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! 背景计算的类
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    private 
    type :: CSFBackground
        !!!!!!! a = a_wkb 时候的 CSF 的场量值
        real(dl) :: A_phi,  B_phi,  C_phi,  D_phi
        !!!!!!!! a = a_ini 的 CSF 的场量值
        real(dl) :: phi10, phi20, dotphi10, dotphi20
        real(dl) :: uphi0(4)        ! 约化场量初始条件

        !!!!!!! 精确计算得到的 CSF 场量
        real(dl), allocatable :: a_ex(:), H_ex(:)
        real(dl), allocatable :: phi1_ex(:), phi2_ex(:)
        real(dl), allocatable :: dotphi1_ex(:), dotphi2_ex(:)

        !!!!!! 数组
        integer :: N_ex, N_wkb, N_tot
        !!!!!! WKB 近似
        real(dl), allocatable :: a_wkb_arr(:)

    contains
    !!!!!!!!!!!!!!! CSF 性质的计算
    procedure :: get_rho
    procedure :: get_rho_arr
    procedure :: get_p
    procedure :: get_p_arr
    procedure :: get_w_phi
    procedure :: get_w_arr
    procedure :: get_cs2_arr
    procedure :: get_a_osc
    !!!!!! 背景演化所需的过程
    procedure :: Initial => CSFBackground_Initial
    !!!! KG 精确方程
    procedure :: Shooting => KGequations_Shooting
    procedure :: Connection
    procedure :: Calc_KGequations
    !!!! WKB 近似解
    procedure :: get_WKB_phi
    procedure :: get_WKB_phi_arr

    !!!! 计算 CSF 背景宇宙演化的状态方程
    procedure :: Calc_CSF_state
    procedure :: get_phi_arr

    end type CSFBackground

    public :: CSFBackground, CalCSFBackground

    
contains

    !!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!!!
    !!! CSF 性质计算
    !!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!!!


    subroutine get_rho(this, phiR, phiI, dotphiR, dotphiI, mu, rho_phi)
    !! 通过 phi 和 dotphi 计算能量密度 rho_phi
    
        class(CSFBackground) :: this
        real(dl), intent(in) :: phiR, phiI, dotphiR, dotphiI   ! phi 的实部，phi 的虚部，dotphi 的实部，dotphi 的虚部
        real(dl), intent(in) :: mu    ! 约化质量 mc^2/hbar
        real(dl), intent(out) :: rho_phi

        rho_phi = 0.5_dl*(dotphiR)**(2._dl) + 0.5_dl * mu**(2._dl)*(phiR)**(2._dl) + 0.5_dl*(dotphiI)**(2._dl) + 0.5_dl * mu**(2._dl)*(phiI)**(2._dl)
    
        
    end subroutine get_rho


    subroutine get_rho_arr(this, phiR_arr, phiI_arr, dotphiR_arr, dotphiI_arr, mu, rho_phi_arr)
    ! 计算 CSF 能量密度的数组 rho_phi_arr
        class(CSFBackground) :: this
        real(dl), intent(in) :: phiR_arr(:), phiI_arr(:), dotphiR_arr(:), dotphiI_arr(:)
        real(dl), intent(in) :: mu     ! 约化质量 mc^2/hbar
        real(dl) , intent(out) :: rho_phi_arr(:)
        integer :: N, ix

        !!!! 判断输入输入的标量场数组实部和虚部维度是否一致
        if (size(phiR_arr) /= size(phiI_arr)) then
            error stop '输入的标量场数组实部和虚部维度不一致！'
        endif

        if (size(dotphiR_arr) /= size(dotphiI_arr)) then
            error stop '输入的标量场的时间导数实部和虚部维度不一致！'
        endif

        !!! 计算 CSF 能量密度的数组
        N = size(phiR_arr)

        do ix = 1,N
            call this%get_rho(phiR_arr(ix), phiI_arr(ix), dotphiR_arr(ix), dotphiI_arr(ix),mu, rho_phi_arr(ix))
        enddo


    end subroutine get_rho_arr


    subroutine get_p(this, phiR, phiI, dotphiR, dotphiI, mu, p_phi)
        !! 通过 phi 和 dotphi 计算压强 p_phi
        
        class(CSFBackground) :: this
        real(dl), intent(in) :: phiR, phiI, dotphiR, dotphiI
        real(dl), intent(in) :: mu    ! 约化质量 mc^2/hbar
        real(dl), intent(out) :: p_phi

        p_phi = 0.5_dl*(dotphiR)**(2._dl) - 0.5_dl * mu**(2._dl)*(phiR)**(2._dl) + 0.5_dl*(dotphiI)**(2._dl) - 0.5_dl * mu**(2._dl)*(phiI)**(2._dl)
        
    end subroutine get_p

    subroutine get_p_arr(this, phiR_arr, phiI_arr, dotphiR_arr, dotphiI_arr, mu, p_phi_arr)
        !! 计算 CSF 压强数组 p_phi_arr
        class(CSFBackground) :: this
        real(dl), intent(in) :: phiR_arr(:), phiI_arr(:), dotphiR_arr(:), dotphiI_arr(:)
        real(dl), intent(in) :: mu    ! 约化质量 mc^2/hbar
        real(dl) , intent(out) :: p_phi_arr(:)
        integer :: N, ix

                !!!! 判断输入输入的标量场数组实部和虚部维度是否一致
        if (size(phiR_arr) /= size(phiI_arr)) then
            error stop '输入的标量场数组实部和虚部维度不一致！'
        endif

        if (size(dotphiR_arr) /= size(dotphiI_arr)) then
            error stop '输入的标量场的时间导数实部和虚部维度不一致！'
        endif

        N = size(phiR_arr)

        do ix = 1,N
            call this%get_p(phiR_arr(ix), phiI_arr(ix), dotphiR_arr(ix), dotphiI_arr(ix), mu, p_phi_arr(ix))
        enddo
        
    end subroutine get_p_arr


    subroutine get_w_phi(this,rho_phi,p_phi,w_phi)
        !! 计算状态方程 w_phi
        class(CSFBackground) :: this
        real(dl), intent(in) :: rho_phi , p_phi
        real(dl), intent(out) :: w_phi

        w_phi = p_phi/rho_phi

    end subroutine get_w_phi

    subroutine get_w_arr(this,rho_phi_arr,p_phi_arr,w_phi_arr)
        !! 计算 w_phi_arr
        class(CSFBackground) :: this
        real(dl), intent(in) :: rho_phi_arr(:), p_phi_arr(:)
        real(dl), intent(out) :: w_phi_arr(:)
        integer :: N, ix

        N = size(rho_phi_arr)

        do ix = 1,N
            call this%get_w_phi(rho_phi_arr(ix),p_phi_arr(ix),w_phi_arr(ix))
        enddo
    
        
    end subroutine get_w_arr


    subroutine get_cs2_arr(this,rho_phi_arr,p_phi_arr,cs2_arr)
        !! 通过插值计算绝热声速 cs^2 = dp/drho
        class(CSFBackground) :: this
        real(dl), intent(in) :: rho_phi_arr(:), p_phi_arr(:)
        real(dl), intent(out) :: cs2_arr(:)
        integer :: n, ix
        real(dl), allocatable, dimension(:) :: g, drho_arr, dp_arr

        n = size(rho_phi_arr)

        allocate(g(n))
        allocate(drho_arr(n))
        allocate(dp_arr(n))

        call splini(g,n)

        call splder(p_phi_arr,dp_arr,n, g)    ! 计算压强的差分的导数 dp_arr
        call splder(rho_phi_arr,drho_arr,n, g)   ! 计算能量密度的导数 drho_arr


        do ix = 1,n
            cs2_arr(ix) = dp_arr(ix)/drho_arr(ix)
        enddo
        
    end subroutine get_cs2_arr


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! 用 LambdaCDM 宇宙模型估算 a_osc 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function LambdaCDM_a_osc(CData, loga) result(fval)
        !!! 用于估算振荡的尺度因子的 LambdaCDM 宇宙模型， a 为尺度因子
        implicit none
        class(CAMBdata) :: CData
        real(dl) :: loga
        real(dl) :: factor = 1.0_dl/3.0_dl
        real(dl) :: a
        real(dl) :: H, mu
        real(dl) :: fval    ! 输出函数值

        !!!!!!!! 读取参数
        mu = CData%CP%DarkMatter%m_phi*eV/hbar

        ! write(*,*) "mu = ",mu

        
        !!!!!! 计算 Hubble 率
        a = 10.0_dl**loga
        H = WKB_confHubble(CData,a)/a

        fval =  H/mu - factor
    
        
    end function LambdaCDM_a_osc

    subroutine get_a_osc(this,CData)
        !!! 用 LambdaCDM 宇宙模型估算 a_osc 
        use MathUtils
        implicit none
        class(CSFBackground) :: this
        class(CAMBdata) :: CData
        real(dl) :: a_osc, a_wkb  ! 开始振荡的尺度因子，输出值
        !!!!!! 辅助变量
        real(dl) :: log_a1, log_a2, fa1, fa2, tol , x0, f0, n   ! 求解区间 (log_a1,log_a2) 及对应函数值 (fa1,fa2)
        integer :: iflag



        !!!! 计算初始区间和对应函数值
        log_a1 = -7.0_dl
        log_a2 = 0.0_dl
        fa1 = LambdaCDM_a_osc(CData, log_a1)
        fa2 = LambdaCDM_a_osc(CData, log_a2)
        tol = 1.0e-9_dl

        ! write(*,*) "fa1 = ", fa1
        ! write(*,*) "fa2 = ", fa2



        if (fa1<0.and.fa2<0) then
            a_osc = 10**(-7.0_dl)
            a_wkb = 100.0_dl*a_osc
        endif

        if (fa1>0.and.fa2>0) then
            a_osc = 0.5
            a_wkb = 2*a_osc
        endif

        !!! 用 Brent 方法找零点

        call brentq(CData,LambdaCDM_a_osc, log_a1, log_a2 , tol , x0 , f0 , iflag , fa1 , fa2 )
        a_osc = 10._dl**(x0)

        a_wkb = 10.0_dl*a_osc


        do while (a_wkb>1._dl)
            a_wkb = a_wkb/2.0_dl
        end do

        CData%CP%DarkMatter%a_osc = a_osc
        CData%CP%DarkMatter%a_wkb = a_wkb

    end subroutine get_a_osc



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! CSFBackground 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine CSFBackground_Initial(this,CData)
        !!!!! 检查 CSF背景条件是否符合要求，并返回 CSF 初始条件
        class(CSFBackground) :: this
        class(CAMBdata) :: CData
        real(dl) :: mu, Q, omphih2, dotR0, R0 , theta0, dottheta0 ,rho_phi0, a_ini             ! CSF 参数
        real(dl) :: phi10 , phi20, dotphi10, dotphi20                ! CSF 场变量
        real(dl) :: Qmax, kconst_phi, L    ! 辅助变量


        !!!!!!!! 读取参数
        mu = CData%CP%DarkMatter%m_phi*ev/hbar                    ! 约化质量
        Q = CData%CP%DarkMatter%Num_density                       ! 守恒荷，共动粒子数密度 Q, Q = Num_density, 单位为 m^-3
        omphih2 = CData%CP%DarkMatter%omphih2
        dotR0 = CData%CP%DarkMatter%dotR
        a_ini = CData%CP%DarkMatter%a_ini

        !!!!!!!! 计算 CSF 能量密度
        kconst_phi = (3.0_dl*c**2.0_dl)/(8.0_dl*const_pi*G)*(10.0_dl**(4.0_dl))*(1000/Mpc)**2.0_dl
        rho_phi0 = kconst_phi*omphih2                                           ! CSF 能量密度

        !!!!!!!! 判断初始条件
        Qmax = rho_phi0/(hbar*mu)                                                          ! Q的绝对值的最大值
 
        ! write(*,*) "Qmax =" , Qmax

        if (abs(Q) >= Qmax) then
            error stop "输入的 CSF 粒子数密度( Num_density )过大，使得计算的能量密度过大，不符合当前宇宙学限制！"
        endif



        !! 设置 theta0
        theta0 = 0.0_dl

        ! !!! 储存初始条件
        CData%CP%DarkMatter%theta0 = theta0
        CData%CP%DarkMatter%dotR0 = dotR0

    end subroutine CSFBackground_Initial

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! KG 方程精确解
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! KG 方程组和弗里德曼方程
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine CSFconfHubble(a, uphi, confHubble, CData)
        !!!! 精确解时期，用于计算共形 Hubble
        use DarkEnergyInterface
        implicit none
        real(dl), intent(in) :: a                        !!! 尺度因子 a
        real(dl), intent(in) :: uphi(4)                  !!! 约化 CSF 场量 uphi = sqrt(4*pi*G/c^2) [phi1, dphi1dtau, ph2 , dphi2dtau ] 
        real(dl), intent(out) :: confHubble
        class(CAMBdata) :: CData
        !!!! 辅助变量
        real(dl) :: detada
        real(dl) :: mu                                   !!! 约化质量， 单位 s^-1
        real(dl) :: grhoa2                               !!! 8*pi*G*rho*a**4.
        real(dl) :: grhov_t, grho_phi

        !!!!!! 读取数据
        mu = CData%CP%DarkMatter%m_phi*eV/hbar


        !!!!!! 计算 grho_phi = 8*pi*G/c^2*a^2*rho_phi
        grho_phi = a**(2._dl)*mu**(2._dl)*uphi(1)**(2._dl) + uphi(2)**(2._dl) + a**(2._dl)*mu**(2._dl)*uphi(3)**(2._dl) + uphi(4)**(2._dl)     

        !!! 获取暗能量
        call CData%CP%DarkEnergy%BackgroundDensityAndPressure(CData%grhov, a, grhov_t)

        !!!  计算 8*pi*G*rho*a**4.
        grhoa2 = CData%grho_no_de(a) +  grhov_t * a**(2._dl)
        grhoa2 = grhoa2 *(c/Mpc)**(2._dl) + grho_phi*a**(2._dl)

        if (grhoa2 <= 0) then
            call GlobalError('Universe stops expanding before today (recollapse not supported)', error_unsupported_params)
            detada = 0
        else
            detada = sqrt(3._dl / grhoa2)
        end if

        confHubble = 1._dl/(a*detada)
    end subroutine CSFconfHubble


    subroutine KGequations( neqn, a, uphi, duphida, CData)
        !!!!! CSF 的严格 Klein-Gordon 方程
        integer, intent(in) :: neqn                 ! 方程个数
        real(dl), intent(in) :: a                   ! 尺度因子的对数 a
        real(dl), intent(in) :: uphi(4)             ! 约化 CSF 场量 uphi = sqrt(4*pi*G/c^2) [phi1, dphi1dtau, ph2 , dphi2dtau ] 
        real(dl), intent(out) :: duphida(4)         ! 约化 CSF 场量 uphi 对尺度因子 a 的一阶导数 duphida
        class(CAMBdata) :: CData                      !
        !!!!! 辅助变量
        real(dl) :: confHubble                      ! 共形 Hubble 率
        real(dl) :: mu                              ! 约化质量

        !!!!!! 读取数据
        mu = CData%CP%DarkMatter%m_phi*eV/hbar

        !!!!! 计算精确的共形 Hubble

        call CSFconfHubble(a, uphi, confHubble, CData)


        !!!!! KG 方程组

        duphida(1) = uphi(2)/(a*confHubble)
        duphida(2) = - 2* uphi(2)/a - a*mu**2*uphi(1)/confHubble
        duphida(3) = uphi(4)/(a*confHubble)
        duphida(4) = - 2* uphi(4)/a - a*mu**2*uphi(3)/confHubble
        
    end subroutine KGequations


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! 打靶法寻找初始条件
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function Shooting_objective_function(CData, x) result(fval)
        !!! 打靶法的目标函数，输入初始条件 R0 和 dotR0，得到 a = a_wkb 时刻的计算目标函数 fval
        !!! 一般目标函数设置为 fval = (rho_wkb-rho_wkb_obj)/rho_wkb_obj ，rho_wkb 为 phi 在 a = a_wkb 时的能量密度和压强
        type(CSFBackground) :: CSFBack
        class(CAMBdata) :: CData                                                 ! 用于计算 CSF 的输入参数
        real(dl), intent(in) :: x                                                  ! 用于打靶的自变量
        real(dl) :: R0, dotR0, theta0, dottheta0
        real(dl) :: mu, Q, omphih2                                                 ! 约化质量， 共动粒子数密度
        real(dl) :: a_ini                                                          ! 开始积分的尺度因子 aini
        real(dl) :: a_wkb                                                          ! 转换到 WKB 近似的尺度因子 a_wkb
        integer, parameter :: neqn = 4                                             ! 微分方程个数
        integer, parameter :: n_step = 100                                         ! 区间数
        integer, parameter :: nlength = 101                                        ! 计算得到的点数
        real(dl) :: uphi0(4)                                                       ! 初始条件
        real(dl) :: phi10, phi20, dotphi10, dotphi20                               ! 初始条件
        type(optionD) :: option                                                    ! 计算精度
        !!!!!!!!!! 计算结果
        real(dl) :: a_phi(nlength)                                                 ! 计算得到的自变量 a_phi
        real(dl) :: uphi(neqn, nlength)                                            ! 计算得到的约化场量
        real(dl) :: phi1wkb, phi2wkb, dotphi1wkb, dotphi2wkb                       ! 衔接条件
        !!!!!!!!! 目标函数
        real(dl) :: fval, rho_wkb, rho_wkb_obj                                     ! 目标函数值及其相关参数
        !!!!!!!!! 辅助变量                                         

        !!!! 读取数据
        a_ini = CData%CP%DarkMatter%a_ini
        a_wkb = CData%CP%DarkMatter%a_wkb
        dotR0 =  CData%CP%DarkMatter%dotR0
        theta0 = CData%CP%DarkMatter%theta0
        mu = CData%CP%DarkMatter%m_phi*eV/hbar
        Q = CData%CP%DarkMatter%Num_density
        omphih2 = CData%CP%DarkMatter%omphih2


        !!! 设置初始条件
        R0 = 10._dl**(x)*(hbar*Q)/(2.0_dl*mu*a_ini**3.0_dl)
        dottheta0 = -(hbar*Q)/(2.0_dl*R0**2.0_dl*a_ini**3._dl)
        !!!!
        phi10 = sqrt(2._dl)*R0
        phi20 = 0._dl
        dotphi10 = sqrt(2._dl)*dotR0
        dotphi20 = sqrt(2._dl)*R0*dottheta0

        uphi0(1) = sqrt(4._dl*const_pi*G/c**(2._dl))*phi10
        uphi0(2) = sqrt(4._dl*const_pi*G/c**(2._dl))*dotphi10*a_ini
        uphi0(3) = sqrt(4._dl*const_pi*G/c**(2._dl))*phi20
        uphi0(4) = sqrt(4._dl*const_pi*G/c**(2._dl))*dotphi20*a_ini

        !!!!!!!!!!!! 设置微分方程精度
        option%abserr = 1.0e-6_dl
        option%flag = 1
        option%relerr = 1.0e-6_dl

        !!!!!!!!!!!! 从 a_ini 反演到 a_wkb 
        call ode45(KGequations, neqn , a_ini, a_wkb, n_step, uphi0, a_phi, uphi, option, CData)


        !!!!!!!!! 计算能量密度
        phi1wkb = sqrt(c**(2._dl)/(4._dl*const_pi*G))*uphi(1,nlength)
        dotphi1wkb = sqrt(c**(2._dl)/(4._dl*const_pi*G))*uphi(2,nlength)/a_wkb
        phi2wkb = sqrt(c**(2._dl)/(4._dl*const_pi*G))*uphi(3,nlength)
        dotphi2wkb = sqrt(c**(2._dl)/(4._dl*const_pi*G))*uphi(4,nlength)/a_wkb

        call CSFBack%get_rho(phi1wkb, phi2wkb, dotphi1wkb, dotphi2wkb, mu, rho_wkb)


        !!!!!!!! 计算目标能量密度和压强
        rho_wkb_obj = (3._dl*c**(2._dl))/(8._dl*const_pi*G)*omphih2*(a_wkb**(-3._dl))*(10.0_dl**4.0_dl)*1000._dl**(2.0_dl)/Mpc**(2._dl)

        ! write(*,*) 'rho_wkb_obj =', rho_wkb_obj
        ! write(*,*) 'rho_wkb =', rho_wkb
        ! fval = ((rho_wkb - rho_wkb_obj)/rho_wkb_obj)**(10._dl)
        fval = log10(rho_wkb/rho_wkb_obj)

    end function Shooting_objective_function


    subroutine KGequations_Shooting(this, CData)
        !!!! 打靶法
        use MathUtils
        implicit none
        class(CSFBackground) :: this
        class(CAMBdata) :: CData                                              ! 用于计算 CSF 的输入参数
        real(dl) :: mu, Q, a_ini
        real(dl) :: R0, dottheta0
        real(dl) :: ax,bx,fax,fbx,tol, fc, fd, xc
        real(dl) :: xzero, fzero
        integer :: iflag
        !!!! 读取数据
        a_ini = CData%CP%DarkMatter%a_ini
        mu = CData%CP%DarkMatter%m_phi*eV/hbar
        Q = CData%CP%DarkMatter%Num_density

        !!!! 确定区间范围
        ax = -10._dl
        bx = 10._dl
        fax = Shooting_objective_function(CData, ax)
        fbx = Shooting_objective_function(CData, bx)
        xc = (ax+bx)/2._dl
        fc = Shooting_objective_function(CData, xc)
        fd = (fbx - fax)/(bx - ax)

        if (fax*fbx>0) then
            if (fc<0) then
                ax = xc
                fax = Shooting_objective_function(CData, ax)
                do while (fax*fbx>0)
                    bx = bx + 5._dl
                    fbx = Shooting_objective_function(CData, bx)
                enddo 
            endif
    
            if (fc>0) then
                if (fd>0) then
                    bx = xc
                    fbx = Shooting_objective_function(CData, bx)
                    do while (fax*fbx>0)
                        ax = ax -5._dl
                        fax = Shooting_objective_function(CData, ax)
                    enddo
                else
                    ax = xc
                    fax = Shooting_objective_function(CData, ax)
                    do while (fax*fbx>0)
                        bx = bx +5._dl
                        fbx = Shooting_objective_function(CData, bx)
                    enddo
                endif
    
            endif    
        endif

        ! write(*,*) "fax =", fax 
        ! write(*,*) "fbx =", fbx

        tol = 1.0e-20_dl

        !!!! 二分法和打靶法
        call brentq(CData, Shooting_objective_function, ax, bx, tol, xzero, fzero, iflag)

        ! write(*,*) "xzero =", xzero
        ! write(*,*) "fzero =", fzero

        !!! 设置初始条件
        R0 = 10**(xzero)*(hbar*Q)/(2.0_dl*mu*a_ini**3.0_dl)
        dottheta0 = -(hbar*Q)/(2.0_dl*R0**2.0_dl*a_ini**3._dl)

        CData%CP%DarkMatter%R0 = R0
        CData%CP%DarkMatter%dottheta0 = dottheta0

    end subroutine KGequations_Shooting


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! 计算初始条件
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Connection(this,CData)
        !!!! 计算初始场量和WKB近似时期的场量
        class(CSFBackground) :: this
        class(CAMBdata) :: CData                                                 ! 用于计算 CSF 的输入参数
        real(dl) :: R0, dotR0, theta0, dottheta0
        real(dl) :: a_ini
        real(dl) :: phi10, phi20, dotphi10, dotphi20
        real(dl) :: uphi0(4)

        !!!! 读取数据
        a_ini = CData%CP%DarkMatter%a_ini
        R0 = CData%CP%DarkMatter%R0
        dotR0 = CData%CP%DarkMatter%dotR0
        theta0 = CData%CP%DarkMatter%theta0
        dottheta0 = CData%CP%DarkMatter%dottheta0

        !!!! 计算初始条件
        phi10 = sqrt(2._dl)*R0
        phi20 = 0._dl
        dotphi10 = sqrt(2._dl)*dotR0
        dotphi20 = sqrt(2._dl)*R0*dottheta0
        

        this%phi10 = phi10
        this%phi20 = phi20
        this%dotphi10 =  dotphi10
        this%dotphi20 = dotphi20


        uphi0(1) = sqrt(4._dl*const_pi*G/c**(2._dl))*phi10
        uphi0(2) = sqrt(4._dl*const_pi*G/c**(2._dl))*dotphi10*a_ini
        uphi0(3) = sqrt(4._dl*const_pi*G/c**(2._dl))*phi20
        uphi0(4) = sqrt(4._dl*const_pi*G/c**(2._dl))*dotphi20*a_ini

        this%uphi0 = uphi0
        

        
    end subroutine Connection


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! 精确求解 KG 方程
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine Calc_KGequations(this, CData)
        !!!!!!! 求解 KG 方程组
        class(CSFBackground) :: this
        class(CAMBdata) :: CData                                               ! 用于计算 CSF 的输入参数
        !!!!!!!!!! 微分方程的所需参数
        real(dl) :: a_ini                                                         ! 开始积分的尺度因子 aini
        real(dl) :: a_wkb                                                         ! 转换到 WKB 近似的尺度因子 a_wkb
        integer, parameter :: neqn = 4                                            ! 微分方程个数
        integer, parameter :: n_step = 4001                                        ! 区间数
        integer, parameter :: nlength = 4002                                       ! 计算得到的点数
        integer, parameter :: n_arr = 4000                                         ! 区间真实的点数
        real(dl) :: uphi0(4)                                                      ! 初始条件
        type(optionD) :: option                                                   ! 计算精度
        !!!!!!!!!! 计算结果
        real(dl) :: a_arr(nlength)                                                ! 计算得到的自变量 a_arr
        real(dl) :: uphi(neqn, nlength)                                           ! 计算得到的约化场量
        !!!!!!!!!! 输出结果
        real(dl) :: H_arr(n_arr)                                                  ! H_arr
        !!!!!!!!!! 辅助变量
        integer :: ix,jx
        real(dl) :: uphi_ix(neqn)
        real(dl) :: phi1wkb, phi2wkb, dotphi1wkb, dotphi2wkb
        real(dl) :: A_phi, B_phi, C_phi, D_phi
        real(dl) :: mu

        !!!!!!!!!!!! 初始条件设置

        a_ini = CData%CP%DarkMatter%a_ini
        a_wkb = CData%CP%DarkMatter%a_wkb
        mu = CData%CP%DarkMatter%m_phi*eV/hbar

        uphi0 = this%uphi0

        !!!!!!!!!!!! 设置微分方程精度
        option%abserr = 1.0e-6_dl
        option%flag = 1
        option%relerr = 1.0e-6_dl

        !!!!!!!!!!!! RK算法积分
        call ode45(KGequations, neqn , a_ini, a_wkb, n_step, uphi0, a_arr, uphi, option, CData)

        !!!!!!!!!!! 计算 CSF 场量
        allocate(this%phi1_ex(n_arr))
        allocate(this%phi2_ex(n_arr))
        allocate(this%dotphi1_ex(n_arr))
        allocate(this%dotphi2_ex(n_arr))
        allocate(this%a_ex(n_arr))
        allocate(this%H_ex(n_arr))


        do ix = 2,n_step
            jx = ix -1
            this%a_ex(jx) = a_arr(ix)
            !!!! 恢复 CSF 场量
            uphi_ix(1) = uphi(1, ix)
            uphi_ix(2) = uphi(2, ix)
            uphi_ix(3) = uphi(3, ix)
            uphi_ix(4) = uphi(4, ix)
            this%phi1_ex(jx) = sqrt(c**(2._dl)/(4._dl*const_pi*G))*uphi_ix(1)
            this%dotphi1_ex(jx) = sqrt(c**(2._dl)/(4._dl*const_pi*G))*uphi_ix(2)/a_arr(ix)
            this%phi2_ex(jx) = sqrt(c**(2._dl)/(4._dl*const_pi*G))*uphi_ix(3)
            this%dotphi2_ex(jx) = sqrt(c**(2._dl)/(4._dl*const_pi*G))*uphi_ix(4)/a_arr(ix)
            !!!! 计算 Hubble 率
            call CSFconfHubble(a_arr(ix), uphi_ix, H_arr(jx), CData)
            H_arr(jx) = H_arr(jx)/a_arr(ix)         !!! 从共形 Hubble 转换成 Hubble
            this%H_ex(jx) = H_arr(jx)*Mpc/1000._dl
        enddo

        phi1wkb = this%phi1_ex(n_arr)
        dotphi1wkb = this%dotphi1_ex(n_arr)
        phi2wkb = this%phi2_ex(n_arr)
        dotphi2wkb = this%dotphi2_ex(n_arr)

        ! write(*,*) "phi1wkb = ", phi1wkb
        ! write(*,*) "phi2wkb = ", phi2wkb
        ! write(*,*) "dotphi1wkb = ", dotphi1wkb
        ! write(*,*) "dotphi2wkb = ", dotphi2wkb

        !!!! 计算 WKB 近似系数
        this%A_phi = phi1wkb*a_wkb**(3.0_dl/2._dl)
        this%B_phi = dotphi1wkb*a_wkb**(3.0_dl/2._dl)/mu
        this%C_phi = phi2wkb*a_wkb**(3.0_dl/2._dl)
        this%D_phi = dotphi2wkb*a_wkb**(3.0_dl/2._dl)/mu

    end subroutine Calc_KGequations


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! WKB 近似解
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    function WKB_confHubble(CData,a) result(confHubble)
        !!!! 当 a >> a_osc, phi 的解用 WKB 近似解代替，Hubble 参数也满足 WKB 近似形似
        use DarkEnergyInterface
        implicit none
        class(CAMBdata) :: CData
        real(dl) :: a
        real(dl) :: detada, grhoa2, grhov_t, grho_phia2, omphih2
        real(dl) :: confHubble

        !!!! 计算 WKB 近似下的 grho_phia2 = 8*pi*G*rho_phi*a^4/c^2.
        omphih2 = CData%CP%DarkMatter%omphih2
        grho_phia2 = 3._dl*((1000.0_dl/Mpc)**2.0_dl)*(10.0_dl**4.0_dl)/(c**(2._dl))*omphih2* a  
        

        call CData%CP%DarkEnergy%BackgroundDensityAndPressure(CData%grhov, a, grhov_t)

        !  8*pi*G*rho*a**4.
        grhoa2 = CData%grho_no_de(a) +  grhov_t * a**(2._dl)
        grhoa2 = grhoa2 *(c/Mpc)**(2._dl) +  grho_phia2

        if (grhoa2 <= 0) then
            call GlobalError('Universe stops expanding before today (recollapse not supported)', error_unsupported_params)
            detada = 0
        else
            detada = sqrt(3._dl / grhoa2)
        end if

        confHubble = 1._dl/(a*detada)
        
    end function WKB_confHubble

    function Integrate_WKB_fun(CData,a) result(fval)
        !!!! 用于 WKB 相位因子积分的辅助函数
        class(CAMBdata) :: CData
        real(dl) :: a, confHubble, fval

        confHubble = WKB_confHubble(CData,a)

        fval = 1.0_dl/confHubble
    end function Integrate_WKB_fun



    subroutine get_WKB_phi(this, CData, a, phi1, phi2, dotphi1, dotphi2, A_phi, B_phi , C_phi, D_phi)
        !!!!!! 用于计算 WKB 近似下的 CSF 场量的值
        class(CSFBackground) :: this
        class(CAMBdata) :: CData
        real(dl), intent(in) :: a
        real(dl) :: mu, a_wkb, a0
        real(dl) :: A_phi, B_phi , C_phi, D_phi, S_phi_a, SinS, CosS, tol
        real(dl), intent(out) :: phi1, phi2, dotphi1, dotphi2

        !!!!! 读取数据
        mu = CData%CP%DarkMatter%m_phi*eV/hbar
        a_wkb = CData%CP%DarkMatter%a_wkb
        a0 = 1.0_dl


        if (a>a0) then
            error stop "输入的 a 应该小于 1!"
        endif

        if (a<a_wkb) then
            error stop "输入的 a < a_wkb, 这不符合 WKB 近似的条件！"
        endif
        !!!!! 利用 Romberg 算法计算相位因子 S
        tol = 1e-6_dl
        S_phi_a = mu*Integrate_Romberg(CData, Integrate_WKB_fun, a_wkb , a, tol)

        SinS = sin(S_phi_a)
        CosS = cos(S_phi_a)

        !!!!! 计算 CSF 的结果

        phi1 = a**(-3.0_dl/2.0_dl) *(A_phi* CosS + B_phi* SinS)
        dotphi1 = mu*a**(-3.0_dl/2.0_dl) * ( -A_phi*SinS + B_phi* CosS)
        phi2 = a**(-3.0_dl/2.0_dl) *(C_phi* CosS + D_phi* SinS)
        dotphi2 = mu*a**(-3.0_dl/2.0_dl) * ( -C_phi*SinS + D_phi* CosS)
        
    end subroutine get_WKB_phi


    subroutine get_WKB_phi_arr(this, CData, a_arr, phi1_arr, phi2_arr, dotphi1_arr, dotphi2_arr)
        !!!!!! 用于计算 WKB 近似下的 CSF 场量的数组
        class(CSFBackground) :: this
        class(CAMBdata) :: CData
        real(dl) :: a_arr(:), phi1_arr(:), phi2_arr(:), dotphi1_arr(:), dotphi2_arr(:)
        integer :: N, ix
        real(dl) :: A_phi, B_phi , C_phi, D_phi

        N = size(a_arr)
        A_phi = this%A_phi
        B_phi = this%B_phi
        C_phi = this%C_phi
        D_phi = this%D_phi

        do ix = 1,N
            call this%get_WKB_phi( CData, a_arr(ix), phi1_arr(ix), phi2_arr(ix), dotphi1_arr(ix), dotphi2_arr(ix), A_phi, B_phi , C_phi, D_phi)
        enddo



    end subroutine get_WKB_phi_arr



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! 计算 CSF 的能量密度，压强，状态方程，绝热声速，以及完整宇宙学背景演化
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Calc_CSF_state(this, CData, N_ex, N_wkb, N_tot)
        class(CSFBackground) :: this
        class(CAMBdata) :: CData
        integer :: N_ex, N_wkb, N_tot
        real(dl) :: mu, omphih2, h
        real(dl) :: rho_ex(N_ex), p_ex(N_ex), w_ex(N_ex), cs2_ex(N_ex)
        real(dl) :: phi1_ex(N_ex), phi2_ex(N_ex), dotphi1_ex(N_ex), dotphi2_ex(N_ex)
        integer :: ix

        !!!!!! 读取数据
        mu = CData%CP%DarkMatter%m_phi*eV/hbar
        omphih2 = CData%CP%DarkMatter%omphih2
        h = CData%CP%H0/100.0_dl

        !!! 读取精确结果
        phi1_ex = this%phi1_ex
        phi2_ex = this%phi2_ex
        dotphi1_ex = this%dotphi1_ex
        dotphi2_ex = this%dotphi2_ex



        !!!! 计算精确解的状态方程
        !!!!!!! 计算能量密度
        call this%get_rho_arr(phi1_ex, phi2_ex, dotphi1_ex, dotphi2_ex, mu, rho_ex)
        !!!!!!! 计算压强
        call this%get_p_arr(phi1_ex, phi2_ex, dotphi1_ex, dotphi2_ex, mu, p_ex)
        !!!!!!! 计算状态方程
        call this%get_w_arr( rho_ex, p_ex, w_ex)
        !!!!!!! 计算绝热声速
        call this%get_cs2_arr(rho_ex, p_ex, cs2_ex)





        !!!! 
        allocate(CData%CP%DarkMatter%a_arr(N_tot))
        allocate(CData%CP%DarkMatter%H_arr(N_tot))
        allocate(CData%CP%DarkMatter%rho_phi(N_tot))
        allocate(CData%CP%DarkMatter%p_phi(N_tot))
        allocate(CData%CP%DarkMatter%w_phi(N_tot))
        allocate(CData%CP%DarkMatter%cs2(N_tot))


        do ix = 1,N_tot
            if(ix<=N_ex) then
                CData%CP%DarkMatter%a_arr(ix) = this%a_ex(ix)
                CData%CP%DarkMatter%H_arr(ix) = this%H_ex(ix)
                CData%CP%DarkMatter%rho_phi(ix) = rho_ex(ix)
                CData%CP%DarkMatter%p_phi(ix) = p_ex(ix)
                CData%CP%DarkMatter%w_phi(ix) = w_ex(ix)
                CData%CP%DarkMatter%cs2(ix) = cs2_ex(ix)
            else
                CData%CP%DarkMatter%a_arr(ix) = this%a_wkb_arr(ix-N_ex)
                CData%CP%DarkMatter%H_arr(ix) = WKB_confHubble(CData,this%a_wkb_arr(ix-N_ex))/this%a_wkb_arr(ix-N_ex)*Mpc/1000.0_dl
                CData%CP%DarkMatter%rho_phi(ix) =(3*c**2)/(8*const_pi*G)* omphih2 * (this%a_wkb_arr(ix-N_ex))**(-3.0_dl)*(10.0_dl**4.0_dl)*1000._dl**(2.0_dl)/Mpc**2
                CData%CP%DarkMatter%p_phi(ix) = 0.0_dl
                CData%CP%DarkMatter%w_phi(ix) = 0.0_dl
                CData%CP%DarkMatter%cs2(ix) = 0.0_dl
            endif
        enddo

        
    
        
    end subroutine Calc_CSF_state





    subroutine get_phi_arr(this, CData, N_ex, N_wkb, N_tot)
        class(CSFBackground) :: this
        class(CAMBdata) :: CData
        integer :: N_ex, N_wkb, N_tot
        integer :: ix
        real(dl) :: a_ex(N_ex), a_wkb(N_wkb)
        real(dl) :: a_arr(N_tot)
        real(dl) :: phi1_ex(N_ex), phi2_ex(N_ex), dotphi1_ex(N_ex), dotphi2_ex(N_ex)
        real(dl) :: phi1_wkb(N_wkb), phi2_wkb(N_wkb), dotphi1_wkb(N_wkb), dotphi2_wkb(N_wkb)

        !!! 读取精确结果
        a_ex = this%a_ex
        phi1_ex = this%phi1_ex
        phi2_ex = this%phi2_ex
        dotphi1_ex = this%dotphi1_ex
        dotphi2_ex = this%dotphi2_ex

        !!!! 读取 WKB 近似范围
        a_wkb = this%a_wkb_arr
        call this%get_WKB_phi_arr(CData, a_wkb, phi1_wkb, phi2_wkb, dotphi1_wkb, dotphi2_wkb)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11111
        ! write(*,*) "phi1_wkb = ", phi1_wkb(1)
        ! write(*,*) "phi2_wkb = ", phi2_wkb(1)
        ! write(*,*) "dotphi1_wkb = ", dotphi1_wkb(1)
        ! write(*,*) "dotphi2_wkb = ", dotphi2_wkb(1)




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        !!! 
        allocate(CData%CP%DarkMatter%phi1_arr(N_tot))
        allocate(CData%CP%DarkMatter%phi2_arr(N_tot))
        allocate(CData%CP%DarkMatter%dotphi1_arr(N_tot))
        allocate(CData%CP%DarkMatter%dotphi2_arr(N_tot))

        do ix = 1,N_tot
            if (ix<=N_ex) then
                CData%CP%DarkMatter%phi1_arr(ix) = phi1_ex(ix)
                CData%CP%DarkMatter%phi2_arr(ix) = phi2_ex(ix)
                CData%CP%DarkMatter%dotphi1_arr(ix) = dotphi1_ex(ix)
                CData%CP%DarkMatter%dotphi2_arr(ix) = dotphi2_ex(ix)
            else
                CData%CP%DarkMatter%phi1_arr(ix) = phi1_wkb(ix-N_ex)
                CData%CP%DarkMatter%phi2_arr(ix) = phi2_wkb(ix-N_ex)
                CData%CP%DarkMatter%dotphi1_arr(ix) = dotphi1_wkb(ix-N_ex)
                CData%CP%DarkMatter%dotphi2_arr(ix) = dotphi2_wkb(ix-N_ex)
            endif
        enddo
        
    end subroutine get_phi_arr


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!  计算 CSF 的程序

    subroutine CalCSFBackground(CP,get_phi)
        !!!!!!! 计算 CSF 背景
        !!!!! 计算所需的类
        logical :: get_phi
        type(CSFBackground) :: CSFBack
        type (CAMBparams) :: CP
        type (CAMBdata) :: CData
        type(TRanges) :: Ra_wkb
        !!!! 计算的辅助参数
        real(dl) :: a0, a_wkb

        call CData%SetParams(CP)
        
        

        !!!!!!!!!!! 计算 a_osc

        !!! 获取 a_osc，并根据 a_osc 计算
        call CSFBack%get_a_osc(CData)    !!!! 获取 a_osc

        


        !!! 初始化
        call CSFBack%Initial(CData)
        !!! 打靶
        call CSFBack%Shooting(CData)
        call CSFBack%Connection(CData)
        call CSFBack%Calc_KGequations(CData)

        !!! 计算状态方程和场量
        !!!!! 初始化数组
        CSFBack%N_ex = size(CSFBack%a_ex)
        CSFBack%N_wkb = 200
        CSFBack%N_tot = CSFBack%N_ex + CSFBack%N_wkb
        allocate(CSFBack%a_wkb_arr(CSFBack%N_wkb))
        a_wkb = CData%CP%DarkMatter%a_wkb
        a0 = 1.0_dl
        call Ra_wkb%Init()
        call Ra_wkb%Add(a_wkb, a0, CSFBack%N_wkb-1,.true.)
        call Ra_wkb%GetArray()
        CSFBack%a_wkb_arr = Ra_wkb%points


        !!!! 计算状态方程
        call CSFBack%Calc_CSF_state(CData, CSFBack%N_ex, CSFBack%N_wkb, CSFBack%N_tot)


        !!!! 是否计算 phi

        if (get_phi) then
            call CSFBack%get_phi_arr(CData, CSFBack%N_ex, CSFBack%N_wkb, CSFBack%N_tot)
        endif

        !!!! 输出计算结果
        CP = CData%CP
        CP%DarkMatter%calculated_CSF = .True.
    end subroutine CalCSFBackground














    !!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc






end module DarkMatterCSF
