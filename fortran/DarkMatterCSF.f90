!!! Equations module allowing for the complex scalar field models
!!! Authors :: Qi Yang(杨祺)
!!!!!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!! 无量纲质量 mnorm = ma*eV/(c*hbar)*Mpc
!!! 无量纲场量 uphi = ( sqrt{4*pi*G/c^4}*phi1 ，sqrt{4*pi*G/c^4}*phi1', sqrt{4*pi*G/c^4}*phi2 ，sqrt{4*pi*G/c^4}*phi2' )
!!! tau 无量纲共形时间


module DarkMatterCSF
    use results
    use DarkMatterInterface
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! 背景计算的类
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    private 
    type :: CSFBackground
        !!!!!!! a = a_wkb 时候的 complex scalar field 的场量值
        real(dl) :: A_wkb,  B_wkb,  C_wkb,  D_wkb
        !!!!!!!! a = a_ini 的 complex scalar field  的场量值
        real(dl) :: uphi0(4)        ! 约化场量初始条件

        !!!!!!! 精确计算得到的 complex scalar field  场量
        real(dl), allocatable :: a_ex(:)
        real(dl), allocatable :: confH_ex(:)
        real(dl), allocatable :: uphi_ex(:,:)
        real(dl), allocatable :: grho_ex(:)
        real(dl), allocatable :: dotgrho_ex(:)
        real(dl), allocatable :: gpres_ex(:)
        real(dl), allocatable :: dotgpres_ex(:)
        real(dl), allocatable :: w_ex(:)
        real(dl), allocatable :: cs2_ex(:)

        !!!!!! WKB 近似得到的 complex scalar field 
        real(dl), allocatable :: a_wkb_arr(:)
        real(dl), allocatable :: confH_wkb(:)
        real(dl), allocatable :: uphi_wkb(:,:)
        real(dl), allocatable :: grho_wkb(:)
        real(dl), allocatable :: dotgrho_wkb(:)
        real(dl), allocatable :: gpres_wkb(:)
        real(dl), allocatable :: dotgpres_wkb(:)
        real(dl), allocatable :: w_wkb(:)
        real(dl), allocatable :: cs2_wkb(:)


    contains
    !!!!!!!!!!!!!!! complex scalar field 演化
    procedure :: Evolution
    !!!! 计算 complex scalar field 背景宇宙演化的状态方程
    procedure :: Full_fluid
    procedure :: Full_CSF
    !!!!!!!!!!!!!!! CSF 性质的计算
    procedure :: calc_dimensionless
    !!!!!!!!!! 能量密度
    procedure :: get_grho_phi
    procedure :: get_grho_arr
    procedure :: get_gdotrho_phi
    procedure :: get_gdotrho_arr
    !!!!!!!!!! 压强
    procedure :: get_gpres_phi
    procedure :: get_gpres_arr
    procedure :: get_gdotgpres_phi
    procedure :: get_gdotgpres_arr
    !!!!!!!!!! 粒子数密度
    procedure :: get_grhon_phi
    procedure :: get_grhon_arr
    !!!!!!!!!! 状态方程
    procedure :: get_w_phi
    procedure :: get_w_arr
    procedure :: get_cs2_arr
    !!!!!! 背景演化所需的过程
    !!!! KG 精确方程
    procedure :: calc_KG
    procedure :: Shooting
    procedure :: ExSolution
    procedure :: Exfluid
    !!!! WKB 近似解
    procedure :: get_WKB_phi
    procedure :: get_WKB_phi_arr
    procedure :: WKBfluid

    end type CSFBackground

    type :: Get_OSC
        !!!!！ 估算标量场的振荡特征
        type(CAMBdata) :: CData
        real(dl) :: mtoH = 1.0_dl
    contains
    procedure :: get_a_osc
    procedure :: mtoHfun
    procedure :: get_mtoH_a

    end type Get_OSC




    public :: CSFBackground, CalCSFBackground

    
contains
    
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!  计算 complex scalar field 的主程序
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    

    subroutine CalCSFBackground(CP,get_phi)
        !!!!!!! 计算 complex scalar field 背景
        implicit none
        type(Get_OSC) :: OSC
        type(CSFBackground) :: CSFBack
        type (CAMBparams) :: CP
        type (CAMBdata) :: CData
        !!!! 计算的辅助参数
        logical, intent(in) :: get_phi
        integer :: neqn,  n_step, N_ex, N_tot, N_wkb
        real(dl) :: a_wkb
        real(dl) :: ain, lna, dlna
        integer :: N, ix
        real(dl) :: uphi(4), rho

        call CData%SetParams(CP)


        !!!!!!!!!!! 计算无量纲量
        call CSFBack%calc_dimensionless(CData)
        
        !!!!!!!!!!! 计算 a_osc
        call OSC%get_a_osc(CData)
        !!!! 计算数组长度
        a_wkb = CData%CP%DarkMatter%a_wkb
        N_tot = CData%CP%DarkMatter%N_tot
        if (abs(a_wkb - 1._dl)<0.01_dl) then
            ! 当 a_wkb 接近 1 的时候，完全用精确解，不用 wkb 近似 
            CData%CP%DarkMatter%a_wkb = 1.0_dl
            N_ex = N_tot
            CData%CP%DarkMatter%N_ex = N_ex
        else
            call NumberSet( CData%CP%DarkMatter, N_tot, N_ex, N_wkb, neqn,  n_step)
            CData%CP%DarkMatter%N_ex = N_ex
        endif
        
        !!!!!!! 计算背景演化
        call CSFBack%Evolution(CData, neqn,  n_step, N_ex, N_wkb, N_tot)


        !!!! 是否计算 phi
        if (get_phi) then
            call CSFBack%Full_CSF(CData%CP%DarkMatter, CData, neqn, N_ex, N_wkb, N_tot)
        endif

        ! N = 10000
        ! lna = - 30._dl
        ! dlna = 30._dl/(N-1)
        ! open(unit=10, file='data.dat', status='replace', action='write')
        ! do ix = 1,N
        !     lna = - 30._dl + real(ix-1)*dlna
        !     ain = 10._dl**lna
        !     rho =  CData%CP%DarkMatter%gpres_phia2(ain)
        !     write(10, *) ain, rho
        ! enddo
        ! close(10)

        !!!! 输出计算结果
        CP = CData%CP
        CP%DarkMatter%calculated_CSF = .True.
    end subroutine CalCSFBackground


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !!!!!
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!!!
    !!! complex scalar field 背景演化
    !!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!!!

    subroutine NumberSet( CSF, N_tot, N_ex, N_wkb, neqn,  n_step)
        !!! 设置方程的参数
        class(TDarkMatterCSF) :: CSF
        integer , intent(in) :: N_tot
        real(dl) :: a_ini, a0, a_wkb
        real(dl) :: L_wkb, L_ex
        integer, intent(out) :: neqn,  n_step, N_wkb, N_ex

        a0 = 1._dl
        a_wkb = CSF%a_wkb
        a_ini = CSF%a_ini

        L_wkb = (log10(a0)-log10(a_wkb))/(log10(a0)-log10(a_ini))
        L_ex =  1- L_wkb


        neqn = 4                   ! KG 方程个数
        N_wkb = int(N_tot*L_wkb)   ! WKB 的数组个数
        N_ex = N_tot - N_wkb       ! 精确解数组数
        n_step = N_ex -1           ! 求解方程的中间步骤    
        
    end subroutine NumberSet


    subroutine Evolution(this, CData, neqn,  n_step, N_ex, N_wkb,  N_tot)
        !!! 背景演化
        implicit none
        class(CSFBackground) :: this
        class(CAMBdata) :: CData
        integer :: neqn, n_step, N_ex, N_wkb, N_tot
        real(dl) :: a_wkb 
        
        a_wkb = CData%CP%DarkMatter%a_wkb

        !!! 打靶
        call this%Shooting(CData)

        !!! 精确解
        call this%ExSolution(CData, neqn,  n_step, N_ex )
        call this%Exfluid(CData, neqn, N_ex)

        !!! 计算 WKB 近似的流体变量
        
        call this%WKBfluid(CData , N_wkb)

        !!! 获得全流体
        call this%Full_fluid(CData%CP%DarkMatter, neqn,  n_step, N_ex, N_wkb, N_tot)



    end subroutine Evolution


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! 计算 complex scalar field 的能量密度，压强，状态方程，绝热声速，以及完整宇宙学背景演化
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Full_fluid(this, CSF, neqn,  n_step, N_ex, N_wkb, N_tot)
        implicit none
        class(CSFBackground) :: this
        class(TDarkMatterCSF) :: CSF
        integer :: neqn,  n_step, N_ex, N_wkb, N_tot
        integer :: ix
        real(dl) :: a_wkb 

        a_wkb = CSF%a_wkb

        !!! 分配内存
        allocate(CSF%a_arr(N_tot))
        allocate(CSF%grho_phi(N_tot))
        allocate(CSF%dotgrho_phi(N_tot))
        allocate(CSF%gpres_phi(N_tot))
        allocate(CSF%dotgpres_phi(N_tot))
        allocate(CSF%w_phi(N_tot))
        allocate(CSF%cs2(N_tot))

        !!!! 将计算结果存入
        
        if (abs(a_wkb-1._dl)<1e-5_dl) then
            !!! 当 a_wkb 接近 1 的时候，完全用精确解代替
            CSF%a_arr = this%a_ex
            CSF%grho_phi = this%grho_ex
            CSF%dotgrho_phi = this%dotgrho_ex
            CSF%gpres_phi = this%gpres_ex
            CSF%dotgpres_phi = this%dotgpres_ex
            CSF%w_phi =  this%w_ex
            CSF%cs2 = this%cs2_ex
            CSF%cs2(1) = this%cs2_ex(2)
        else
            !!! 精确解
            CSF%a_arr(1:N_ex) = this%a_ex
            CSF%grho_phi(1:N_ex) = this%grho_ex
            CSF%dotgrho_phi(1:N_ex) = this%dotgrho_ex
            CSF%gpres_phi(1:N_ex) = this%gpres_ex
            CSF%dotgpres_phi(1:N_ex) = this%dotgpres_ex
            CSF%w_phi(1:N_ex) =  this%w_ex
            CSF%cs2(1:N_ex) = this%cs2_ex
            CSF%cs2(1) = this%cs2_ex(2)
            !!! WKB 近似
            CSF%a_arr(N_ex+1:N_tot) = this%a_wkb_arr
            CSF%grho_phi(N_ex+1:N_tot)= this%grho_wkb
            CSF%dotgrho_phi(N_ex+1:N_tot)= this%dotgrho_wkb
            CSF%gpres_phi(N_ex+1:N_tot) = this%gpres_wkb
            CSF%dotgpres_phi(N_ex+1:N_tot) = this%dotgpres_wkb
            CSF%w_phi(N_ex+1:N_tot) =  this%w_wkb
            CSF%cs2(N_ex+1:N_tot) = this%cs2_wkb
        endif


    end subroutine Full_fluid


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! 计算 complex scalar field 宇宙的完整演化
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine Full_CSF(this, CSF, CData, neqn, N_ex, N_wkb, N_tot)
        !!!! 计算 complex scalar field 宇宙的背景 Hubble 和 complex scalar field 场量演化
        class(CSFBackground) :: this
        class(TDarkMatterCSF) :: CSF
        class(CAMBdata) :: CData
        integer :: neqn, N_ex, N_wkb, N_tot
        integer :: ix
        real(dl) :: confH_wkb
        real(dl) :: uphi(neqn,N_wkb)
        real(dl) :: uphiT(neqn,N_tot)

        !!! 分配内存
        allocate(CSF%H_arr(N_tot))
        allocate(CSF%phi1(N_tot))
        allocate(CSF%dotphi1(N_tot))
        allocate(CSF%phi2(N_tot))
        allocate(CSF%dotphi2(N_tot))
        allocate(CSF%grhon_phi(N_tot))

        !!! 计算 WKB 近似的场
        call this%get_WKB_phi_arr( CData, neqn, N_wkb, this%a_wkb_arr, uphi)


        !!! 精确解
        CSF%H_arr(1:N_ex) = this%confH_ex/this%a_ex
        CSF%phi1(1:N_ex) = this%uphi_ex(1,:)
        CSF%dotphi1(1:N_ex) = this%uphi_ex(2,:)
        CSF%dotphi1(1) = this%uphi_ex(2,2)

        CSF%phi2(1:N_ex) = this%uphi_ex(3,:)
        CSF%dotphi2(1:N_ex) = this%uphi_ex(4,:)
        CSF%dotphi2(1) = this%uphi_ex(4,2)

        !!! WKB 近似解
        CSF%H_arr(N_ex+1:N_tot) = this%confH_wkb/this%a_wkb_arr
        CSF%phi1(N_ex+1:N_tot) = uphi(1,:)
        CSF%dotphi1(N_ex+1:N_tot) = uphi(2,:)
        CSF%phi2(N_ex+1:N_tot) = uphi(3,:)
        CSF%dotphi2(N_ex+1:N_tot) = uphi(4,:)

        !!! 计算粒子数密度
        uphiT(1,:) = CSF%phi1
        uphiT(2,:) = CSF%dotphi1
        uphiT(3,:) = CSF%phi2
        uphiT(4,:) = CSF%dotphi2
        call this%get_grhon_arr(neqn, N_tot, CSF%a_arr, uphiT, CSF%mnorm, CSF%Nnorm, CSF%grhon_phi)

        
        
    end subroutine Full_CSF




    !!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!!!
    !!! complex scalar field 性质计算
    !!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!!!


    subroutine calc_dimensionless(this, CData)
        !!! 计算无量纲量
        use constants, only : c, eV, Mpc, hbar, kappa
        implicit none
        class(CSFBackground) :: this
        class(CAMBdata) :: CData
        real(dl) :: a_ini, mnorm, Nnorm, dotR_ini, dotR_norm 
        real(dl), parameter :: metertoeV = 5067730.716156395  ! meter 转eV^-1
        real(dl), parameter :: sectoeV = eV/hbar             ! sec 转 eV^-1
        real(dl), parameter :: kgtoeV = c**2/eV              ! kg转eV
                  

        !!! 计算无量纲质量
        mnorm = CData%CP%DarkMatter%m_phi*eV/(c*hbar)*Mpc  ! 无量纲质量 mnorm = m_phi*eV/(c*hbar)*Mpc
        CData%CP%DarkMatter%mnorm = mnorm
        ! write(*,*) 'mnorm =', mnorm

        !!! 计算无量纲共动粒子数密度
        ! 无量纲共动粒子数密度 Nnorm = N*10^6*4*pi*G*hbar/c^3*Mpc
        Nnorm = CData%CP%DarkMatter%Num_density*10._dl**(6._dl)    ! 转换单位，cm^-3 转换为 m^-3
        Nnorm = 0.5_dl*Nnorm*kappa*hbar/c**3*Mpc
        CData%CP%DarkMatter%Nnorm = Nnorm
        ! write(*,*) 'Nnorm = ', Nnorm

        !!! 计算无量纲径向速度 dRdtau
        a_ini = CData%CP%DarkMatter%a_ini
        dotR_ini = CData%CP%DarkMatter%dotR_ini              ! 读取数据, eV^2
        dotR_norm  = sqrt(kappa/c**4._dl)*sqrt(metertoeV *sectoeV**2/kgtoeV)*a_ini*dotR_ini*Mpc   ! 无量纲化
        CData%CP%DarkMatter%dotR_norm = dotR_norm

        ! write(*,*) 'dotR_norm  =', dotR_norm 
    
        
    end subroutine calc_dimensionless


    !!!!!!!!!!!!!!!!!!!
    !!! 能量密度
    !!!!!!!!!!!!!!!!!!!
    subroutine get_grho_phi(this, a, uphi, mnorm, grho_phi)
    !! 计算无量纲物质密度 grho_phi
        implicit none
        class(CSFBackground) :: this
        real(dl), intent(in) :: a             ! 尺度因子
        real(dl), intent(in) :: uphi(4)       ! 无量纲场量
        real(dl), intent(in) :: mnorm         ! 无量纲质量 mnorm = m_phi*eV/(c*hbar)*Mpc
        real(dl), intent(out) :: grho_phi      ! grho_phi = kappa*a^2*rho_phi/c^2 * Mpc^2/c^2

        ! grho_phi = 8*pi*G/c^4*(|phi'|^2 + a^2*mnorm^2*|phi|^2)
        grho_phi = uphi(2)**2 + a**2*mnorm**2*uphi(1)**2 + uphi(4)**2 + a**2*mnorm**2*uphi(3)**2   
        
    end subroutine get_grho_phi


    subroutine get_grho_arr(this, neqn, N, a_arr, uphi, mnorm, grho_arr)
        !! 计算 complex scalar field 无量纲物质密度的数组 grho_arr
        implicit none
        class(CSFBackground) :: this
        integer,  intent(in) :: neqn, N                     ! 自变量个数, 数组大小
        real(dl), intent(in) :: a_arr(N)                    ! 尺度因子数组
        real(dl), intent(in) :: uphi(neqn, N)               ! 无量纲场量
        real(dl), intent(in) :: mnorm                       !无量纲质量 mnorm = m_phi*eV/(c*hbar)*Mpc
        real(dl) , intent(out) :: grho_arr(N)
        integer :: ix
        real(dl) :: uphia(neqn)
        
        
        !!! 计算 complex scalar field 能量密度的数组
        do ix = 1, N
            uphia(1) = uphi(1, ix)
            uphia(2) = uphi(2, ix)
            uphia(3) = uphi(3, ix)
            uphia(4) = uphi(4, ix)
            call this%get_grho_phi(a_arr(ix), uphia, mnorm, grho_arr(ix))
        enddo
        
        
    end subroutine get_grho_arr


    subroutine get_gdotrho_phi(this, a, uphi, confH, gdotrho_phi)
    !! 计算无量纲物质密度的时间导数 gdotrho_phi
        implicit none
        class(CSFBackground) :: this
        real(dl), intent(in) :: a             ! 尺度因子
        real(dl), intent(in) :: confH         ! 无量纲共形哈勃
        real(dl), intent(in) :: uphi(4)       ! 无量纲场量
        real(dl), intent(out) :: gdotrho_phi  ! gdotrho_phi = kappa*a^2*rho_phi'/c^2 * Mpc^2/c^2

        gdotrho_phi = -6._dl*confH*(uphi(2)**2 + uphi(4)**2)
        
    end subroutine get_gdotrho_phi


    subroutine get_gdotrho_arr(this, neqn, N, a_arr, uphi, confH_arr, gdotrho_arr)
        !! 计算 complex scalar field 无量纲物质密度时间导数的数组 gdotrho_arr
        implicit none
        class(CSFBackground) :: this
        integer,  intent(in) :: neqn, N                     ! 自变量个数, 数组大小
        real(dl), intent(in) :: a_arr(N)                    ! 尺度因子数组
        real(dl), intent(in) :: uphi(neqn, N)               ! 无量纲场量
        real(dl), intent(in) :: confH_arr(N)                     ! 无量纲共形 Hubble
        real(dl) , intent(out) :: gdotrho_arr(N)
        integer :: ix
        real(dl) :: uphia(neqn)
        
        
        !!! 计算 complex scalar field 能量密度的数组
        do ix = 1, N
            uphia(1) = uphi(1, ix)
            uphia(2) = uphi(2, ix)
            uphia(3) = uphi(3, ix)
            uphia(4) = uphi(4, ix)
            call this%get_gdotrho_phi(a_arr(ix), uphia, confH_arr(ix), gdotrho_arr(ix))
        enddo
        
        
    end subroutine get_gdotrho_arr

    !!!!!!!!!!!!!!!!!!!
    !!! 压强
    !!!!!!!!!!!!!!!!!!!
    subroutine get_gpres_phi(this, a, uphi, mnorm, gpres_phi)
        !! 计算无量纲压强 gpres_phi
        implicit none
        class(CSFBackground) :: this
        real(dl), intent(in) :: a             ! 尺度因子
        real(dl), intent(in) :: uphi(4)       ! 无量纲场量
        real(dl), intent(in) :: mnorm         ! 无量纲质量 mnorm = m_phi*eV/(c*hbar)*Mpc
        real(dl), intent(out) :: gpres_phi    ! grho_phi = kappa*a^2*rho_phi/c^2 * Mpc^2/c^2
        
        ! gpres_phi = 8*pi*G/c^4*(|phi'|^2 - a^2*mnorm^2*|phi|^2)
        gpres_phi = uphi(2)**2 - a**2*mnorm**2*uphi(1)**2 + uphi(4)**2 - a**2*mnorm**2*uphi(3)**2   
                
    end subroutine get_gpres_phi

    subroutine get_gpres_arr(this, neqn, N, a_arr, uphi, mnorm, gpres_arr)
        !! 计算 complex scalar field 无量纲压强的数组 gpres_arr
        implicit none
        class(CSFBackground) :: this
        integer,  intent(in) :: neqn, N                     ! 自变量个数, 数组大小
        real(dl), intent(in) :: a_arr(N)                    ! 尺度因子数组
        real(dl), intent(in) :: uphi(neqn, N)               ! 无量纲场量
        real(dl), intent(in) :: mnorm                       ! 无量纲质量 mnorm = m_phi*eV/(c*hbar)*Mpc
        real(dl) , intent(out) :: gpres_arr(N)
        integer :: ix
        real(dl) :: uphia(neqn)
        
        
        !!! 计算 complex scalar field 能量密度的数组
        do ix = 1, N
            uphia(1) = uphi(1, ix)
            uphia(2) = uphi(2, ix)
            uphia(3) = uphi(3, ix)
            uphia(4) = uphi(4, ix)
            call this%get_gpres_phi(a_arr(ix), uphia, mnorm, gpres_arr(ix))
        enddo
        
        
    end subroutine get_gpres_arr


    subroutine get_gdotgpres_phi(this, a, uphi, mnorm, confH, gdotgpres_phi)
        !! 计算无量纲压强的时间导数 gdotgpres_phi
            implicit none
            class(CSFBackground) :: this
            real(dl), intent(in) :: a                ! 尺度因子
            real(dl), intent(in) :: confH            ! 无量纲共形哈勃
            real(dl), intent(in) :: uphi(4)          ! 无量纲场量
            real(dl), intent(in) :: mnorm            ! 无量纲质量 mnorm = m_phi*eV/(c*hbar)*Mpc
            real(dl), intent(out) :: gdotgpres_phi   ! gdotgpres_phi  = kappa*a^2*p_phi'/c^2 * Mpc^2/c^2
    
            gdotgpres_phi  = -6._dl*confH*(uphi(2)**2 + uphi(4)**2) - 4._dl*mnorm**2*a**2*(uphi(1)*uphi(2)+uphi(3)*uphi(4))
            
    end subroutine get_gdotgpres_phi


    subroutine get_gdotgpres_arr(this, neqn, N, a_arr, uphi, mnorm, confH_arr, gdotgpres_arr)
        !! 计算 complex scalar field 无量纲压强时间导数的数组 gdotgpres_arr
        implicit none
        class(CSFBackground) :: this
        integer,  intent(in) :: neqn, N                     ! 自变量个数, 数组大小
        real(dl), intent(in) :: a_arr(N)                    ! 尺度因子数组
        real(dl), intent(in) :: uphi(neqn, N)               ! 无量纲场量
        real(dl), intent(in) :: mnorm                       ! 无量纲质量 mnorm = m_phi*eV/(c*hbar)*Mpc
        real(dl), intent(in) :: confH_arr(N)                ! 无量纲共形 Hubble
        real(dl) , intent(out) :: gdotgpres_arr(N)
        integer :: ix
        real(dl) :: uphia(neqn)
        
        
        do ix = 1, N
            uphia(1) = uphi(1, ix)
            uphia(2) = uphi(2, ix)
            uphia(3) = uphi(3, ix)
            uphia(4) = uphi(4, ix)
            call this%get_gdotgpres_phi(a_arr(ix), uphia, mnorm, confH_arr(ix), gdotgpres_arr(ix))
        enddo
        
        
    end subroutine get_gdotgpres_arr

    !!!!!!!!!!!!!!!!!!!
    !!! 粒子数密度
    !!!!!!!!!!!!!!!!!!!
    subroutine get_grhon_phi(this, a, uphi, mnorm, Nnorm, grhon_phi)
    !! 计算无量纲粒子数密度 grhon_phi
        implicit none
        class(CSFBackground) :: this
        real(dl), intent(in) :: a             ! 尺度因子
        real(dl), intent(in) :: uphi(4)       ! 无量纲场量
        real(dl), intent(in) :: mnorm         ! 无量纲质量 mnorm = m_phi*eV/(c*hbar)*Mpc
        real(dl), intent(in) :: Nnorm         !
        real(dl), intent(out) :: grhon_phi      ! grhon_phi = kappa*a^2*rho_phi/c^2 * Mpc^2/c^2

        ! grhon_phi = 8*pi*G*a^2*m_phi*rho_n*Mpc^2/c^2
        ! grhon_phi = -2._dl*mnorm*a*(uphi(1)*uphi(4)-uphi(3)*uphi(2))
        grhon_phi = 2._dl*mnorm*Nnorm/a
        
    end subroutine get_grhon_phi


    subroutine get_grhon_arr(this, neqn, N, a_arr, uphi, mnorm, Nnorm, grhon_arr)
        !! 计算 complex scalar field 无量纲粒子数密度的数组 grhon_arr
        implicit none
        class(CSFBackground) :: this
        integer,  intent(in) :: neqn, N                     ! 自变量个数, 数组大小
        real(dl), intent(in) :: a_arr(N)                    ! 尺度因子数组
        real(dl), intent(in) :: uphi(neqn, N)               ! 无量纲场量
        real(dl), intent(in) :: mnorm                       !无量纲质量 mnorm = m_phi*eV/(c*hbar)*Mpc
        real(dl), intent(in) :: Nnorm  
        real(dl) , intent(out) :: grhon_arr(N)
        integer :: ix
        real(dl) :: uphia(neqn)
        
        !!! 计算 complex scalar field 能量密度的数组
        do ix = 1, N
            uphia(1) = uphi(1, ix)
            uphia(2) = uphi(2, ix)
            uphia(3) = uphi(3, ix)
            uphia(4) = uphi(4, ix)
            call this%get_grhon_phi(a_arr(ix), uphia, mnorm, Nnorm, grhon_arr(ix))
        enddo
        
        
    end subroutine get_grhon_arr


    !!!!!!!!!!!!!!!!!!!
    !!! 状态方程
    !!!!!!!!!!!!!!!!!!!
    subroutine get_w_phi(this, grho_phi, gpres_phi, w_phi)
        !! 计算状态方程 w_phi
        class(CSFBackground) :: this
        real(dl), intent(in) :: grho_phi, gpres_phi
        real(dl), intent(out) :: w_phi

        w_phi = gpres_phi/grho_phi

    end subroutine get_w_phi

    subroutine get_w_arr(this, grho_arr, gpres_arr, w_phi_arr)
        !! 计算 w_phi_arr
        class(CSFBackground) :: this
        real(dl), intent(in) :: grho_arr(:), gpres_arr(:)
        real(dl), intent(out) :: w_phi_arr(:)
        integer :: N, ix

        N = size(grho_arr)

        do ix = 1,N
            call this%get_w_phi(grho_arr(ix), gpres_arr(ix), w_phi_arr(ix))
        enddo
    
    end subroutine get_w_arr


    subroutine get_cs2_arr(this, a_arr, neqn, N, uphi, mnorm, confH_arr, cs2_arr)
        !!! 计算绝热声速 cs2
        class(CSFBackground) :: this
        real(dl), intent(in) :: a_arr(N)
        integer,  intent(in) :: neqn, N                          ! 自变量个数
        real(dl), intent(in) :: uphi(neqn, N)                    ! 无量纲场量
        real(dl), intent(in) :: mnorm                            ! 无量纲质量 mnorm = ma*eV/(c*hbar)*Mpc
        real(dl), intent(in) :: confH_arr(N)                     ! 无量纲共形 Hubble
        integer :: ix
        real(dl) :: uphia(neqn)
        real(dl), intent(out) :: cs2_arr(N)
        real(dl) :: cs2max = 10

        do ix = 1, N
            uphia(1) = uphi(1,ix)
            uphia(2) = uphi(2,ix)
            uphia(3) = uphi(3,ix)
            uphia(4) = uphi(4,ix)

            cs2_arr(ix) = 1._dl + &
            (2._dl*mnorm**2._dl*a_arr(ix)**2._dl*( uphia(1)*uphia(2) + uphia(4)*uphia(3) ))/( 3._dl*confH_arr(ix)*( uphia(2)**2 + uphia(4)**2 ) )

            
        enddo
    end subroutine get_cs2_arr


    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!! CSFBackground 
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!! KG 方程精确解
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!! KG 方程组和弗里德曼方程
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine CSFconfHubble(a, uphi, confHubble, CData)
        !!!! 精确解时期，用于计算共形 Hubble
        use DarkEnergyInterface
        implicit none
        type(CSFBackground) :: CSFBack
        real(dl), intent(in) :: a                        !!! 尺度因子 a
        real(dl), intent(in) :: uphi(4)                  !!! complex scalar field 无量纲场量
        real(dl), intent(out) :: confHubble              !!! 无量纲共形 Hubble
        class(CAMBdata) :: CData
        !!!! 辅助变量
        real(dl) :: dtauda
        real(dl) :: mnorm                                !!! 无量纲质量 mnorm = ma*eV/(c*hbar)*Mpc
        real(dl) :: grhoa2                               !!! 8*pi*G*rho*a**4.
        real(dl) :: grhov_t, grho_phi

        !!!!!! 读取质量
        mnorm = CData%CP%DarkMatter%mnorm
        !!!!!! 计算 grho_phi = 8*pi*G/c^2*a^2*rho_phi
        call CSFBack%get_grho_phi(a, uphi, mnorm, grho_phi)   

        !!! 获取暗能量
        call CData%CP%DarkEnergy%BackgroundDensityAndPressure(CData%grhov, a, grhov_t)

        !!!  计算 8*pi*G*rho*a**4.
        grhoa2 = CData%grho_no_de(a) +  grhov_t * a**(2._dl)
        grhoa2 = grhoa2  + grho_phi*a**(2._dl)

        if (grhoa2 <= 0) then
            call GlobalError('Universe stops expanding before today (recollapse not supported)', error_unsupported_params)
            dtauda = 0
        else
            dtauda = sqrt(3._dl / grhoa2)
        end if

        confHubble = 1._dl/(a*dtauda)
    end subroutine CSFconfHubble


    subroutine KGequations( neqn, lna, uphi, duphidx, CData)
        !!!!! complex scalar field 的严格 Klein-Gordon 方程
        integer, intent(in) :: neqn                 ! 方程个数
        real(dl), intent(in) :: lna                 ! 演化自变量 lna = ln(a)
        real(dl) :: a                               ! 尺度因子 a
        real(dl) :: a2                              ! 尺度因子平方 a2 =a^2
        real(dl), intent(in) :: uphi(4)             ! complex scalar field  无量纲场量 uphi
        real(dl), intent(out) :: duphidx(4)         ! complex scalar field  无量纲场量 uphi 对 x 的一阶导数 duphidx
        class(CAMBdata) :: CData                    !
        !!!!! 辅助变量
        real(dl) :: confHubble                      ! 共形 Hubble 率
        real(dl) :: mnorm                           ! 无量纲质量 mnorm = ma*eV/(c*hbar)*Mpc

        a = exp(lna)
        a2 = a*a

        !!!!!! 读取数据
        mnorm = CData%CP%DarkMatter%mnorm

        !!!!! 计算精确的共形 Hubble
        call CSFconfHubble(a, uphi, confHubble, CData)


        !!!!! KG 方程组

        duphidx(1) = uphi(2)/confHubble
        duphidx(2) = - 2* uphi(2) - a2*mnorm**2*uphi(1)/confHubble
        duphidx(3) = uphi(4)/confHubble
        duphidx(4) = - 2* uphi(4) - a2*mnorm**2*uphi(3)/confHubble
        
    end subroutine KGequations

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! 精确求解 KG 方程
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine calc_KG(this, CData, lna_ini, lna_end, neqn, uphi0, n_step, nlength, a_arr, uphi)
        !!!!!!! 求解 KG 方程组
        use RKF45
        implicit none
        class(CSFBackground) :: this
        class(CAMBdata) :: CData                                               
        !!!!!!!!!! 微分方程的所需参数
        real(dl), intent(in) :: lna_ini               ! 开始积分的尺度因子 aini
        real(dl), intent(in) :: lna_end               ! 结束积分的尺度因子 a_end
        integer , intent(in) :: neqn                  ! 微分方程个数
        real(dl), intent(in) :: uphi0(neqn)           ! 初始条件
        integer , intent(in) :: n_step                ! 区间数
        integer , intent(in) :: nlength               ! 计算得到的点数
        
        type(optionD) :: option                       ! 计算精度
        !!!!!!!!!! 计算结果
        real(dl) :: lna_arr(nlength)                  ! 计算得到的自变量 lna_arr
        real(dl), intent(out) :: a_arr(nlength)       ! a_arr = exp(lna_arr)
        real(dl), intent(out) :: uphi(neqn, nlength)  ! 计算得到的约化场量
        !!!!!!!!! 辅助变量

        !!!!!!!!!!!! 设置微分方程精度
        option%abserr = 1.0e-10_dl
        option%flag = 1
        option%relerr = 1.0e-15_dl

        !!!!!!!!!!!! RK算法积分
        call oderk45(CData, KGequations, neqn , lna_ini, lna_end, n_step, uphi0, lna_arr, uphi, option )
        a_arr = exp(lna_arr)
    end subroutine calc_KG

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!!! 打靶法寻找初始条件
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function Shooting_objective_function(CData, x) result(fval)
        !!! 打靶法的目标函数，输入初始条件 R0 和 dotR0，得到 a = a_wkb 时刻的计算目标函数 fval
        implicit none
        type(CSFBackground) :: CSFBack
        class(CAMBdata) :: CData                                                   ! 用于计算 CSF 的输入参数
        real(dl), intent(in) :: x                                                  ! 用于打靶的自变量，R_ini = 10^x
        real(dl) :: mnorm, omphih2, Q,  dotR_norm                                  ! 无量纲质量，能量密度比例 ，共动粒子数密度
        real(dl) :: a_ini, lna_ini                                                 ! 开始积分的尺度因子 aini
        real(dl) :: a_wkb, lna_wkb                                                 ! 转换到 WKB 近似的尺度因子 a_wkb
        integer, parameter :: neqn = 4                                             ! 微分方程个数
        integer, parameter :: n_step = 1000                                        ! 区间数
        integer, parameter :: nlength = 1001                                       ! 计算得到的点数
        real(dl) :: uphi0(neqn)                                                    ! 初始条件
        !!!!!!!!!! 计算结果
        real(dl) :: a_arr(nlength)                                                 ! 计算得到的自变量 a_phi
        real(dl) :: uphi(neqn, nlength)                                            ! 计算得到的约化场量
        real(dl) :: uphi_wkb(neqn)                                                 ! 计算得到的最后精确值
        !!!!!!!!! 目标函数
        real(dl) :: fval, grhophi, grhophi_obj                                     ! 目标函数值及其相关参数
        !!!!!!!!! 辅助变量
        real(dl) :: R_ini, dotR_ini, theta_ini, dottheta_ini
        real(dl) :: h                                                    

        !!!! 读取数据
        h = CData%CP%H0/100._dl
        a_ini = CData%CP%DarkMatter%a_ini
        lna_ini = log(a_ini)
        a_wkb = CData%CP%DarkMatter%a_wkb
        lna_wkb = log(a_wkb)

        mnorm = CData%CP%DarkMatter%mnorm
        omphih2 = CData%CP%DarkMatter%omphih2
        Q = CData%CP%DarkMatter%Nnorm
        dotR_norm = CData%CP%DarkMatter%dotR_norm
        
        


        !!! 设置初始条件
        R_ini = 10**x
        dotR_ini = dotR_norm
        theta_ini = 0._dl
        dottheta_ini = -Q/(R_ini**2*a_ini**2)


        uphi0(1) = R_ini
        uphi0(2) = dotR_ini
        uphi0(3) = 0._dl
        uphi0(4) = R_ini*dottheta_ini 
        ! write(*,*) uphi0

        !!!!!!!! 解 KG 方程
        call CSFBack%calc_KG( CData, lna_ini, lna_wkb, neqn, uphi0, n_step, nlength, a_arr, uphi)

        !!!!!! 计算 grho_phi = 8*pi*G/c^2*a^2*rho_phi
        uphi_wkb(1) = uphi(1,nlength )
        uphi_wkb(2) = uphi(2,nlength )
        uphi_wkb(3) = uphi(3,nlength )
        uphi_wkb(4) = uphi(4,nlength )
        call CSFBack%get_grho_phi(a_arr(nlength), uphi_wkb, mnorm, grhophi)   

        ! write(*,*) 'grhophi=', grhophi

        !!!!!!!! 计算目标能量密度
        !!!! 计算 WKB 近似下的 grhophia2 = 8*pi*G*rho_phi*a^4/c^2.
        grhophi_obj = CData%grhocrit*omphih2/h**2/ a_wkb

        ! write(*,*) 'grhophi_obj =', grhophi_obj

        !!!!!! 计算打靶的目标函数

        fval = log10(grhophi/grhophi_obj)
        ! write(*,*) 'fval = ', fval

    end function Shooting_objective_function

    subroutine Shooting(this, CData)
        !!!! 打靶法
        use MathUtils
        implicit none
        class(CSFBackground) :: this
        class(CAMBdata) :: CData
        real(dl) :: ax, bx, fax, fbx, tol
        real(dl) :: xzero, fzero
        integer :: iflag
        !!! 辅助变量
        real(dl) :: a_ini, Q,  dotR_norm
        real(dl) :: R_ini, dotR_ini, theta_ini, dottheta_ini

        !!! 读数据
        a_ini = CData%CP%DarkMatter%a_ini
        Q = CData%CP%DarkMatter%Nnorm
        dotR_norm = CData%CP%DarkMatter%dotR_norm
        
        ! write(*,*) 'a_ini', a_ini


        !!!! 确定区间范围
        ax = 10
        bx = -10
        fax = Shooting_objective_function(CData, ax)
        fbx = Shooting_objective_function(CData, bx)


        ! write(*,*) "fax =", fax 
        ! write(*,*) "fbx =", fbx

        if (fax*fbx>0) then
            if(abs(fax)<abs(fbx)) then
                do while (fax*fbx>0)
                    ax = ax - 5._dl
                    fax = Shooting_objective_function(CData, ax)
                enddo 
            else
                do while (fax*fbx>0)
                    bx = bx + 5._dl
                    fbx = Shooting_objective_function(CData, bx)
                enddo 
            endif
        endif

        ! write(*,*) "fax =", fax 
        ! write(*,*) "fbx =", fbx

        !! 打靶精度
        tol = 1.0e-15_dl

        !!!! 二分法和打靶法
        call brentq(CData, Shooting_objective_function, ax, bx, tol, xzero, fzero, iflag)

        ! write(*,*) "xzero =", xzero
        ! write(*,*) "fzero =", fzero


        !!! 设置初始条件
        
        R_ini = 10**xzero
        dotR_ini = dotR_norm
        theta_ini = 0._dl
        dottheta_ini = -Q/(R_ini**2*a_ini**2)
        ! write(*,*) 'dottheta_ini=', dottheta_ini


        this%uphi0(1) = R_ini
        this%uphi0(2) = dotR_ini
        this%uphi0(3) = 0._dl
        this%uphi0(4) = R_ini*dottheta_ini 
        ! write(*,*) this%uphi0


    end subroutine Shooting

    subroutine ExSolution(this,  CData, neqn,  n_step, nlength )
        implicit none
        class(CSFBackground) :: this
        class(CAMBdata) :: CData
        integer :: neqn,  n_step, nlength
        real(dl) :: a_ini, lna_ini            ! 开始积分的尺度因子 aini
        real(dl) :: a_wkb, lna_wkb            ! 转换到 WKB 近似的尺度因子 a_wkb
        real(dl) :: mnorm                     ! 无量纲质量
        real(dl) :: uphi0(neqn)               ! 初始条件
        !!! 辅助变量
        integer :: ix

        !!! 读取数据
        mnorm = CData%CP%DarkMatter%mnorm
        a_ini = CData%CP%DarkMatter%a_ini
        lna_ini = log(a_ini)
        a_wkb = CData%CP%DarkMatter%a_wkb
        lna_wkb = log(a_wkb)

        !!! 设置初始条件
        uphi0 = this%uphi0
        ! write(*,*) this%uphi0

        !!! 分配内存
        allocate(this%a_ex(nlength))
        allocate(this%confH_ex(nlength))
        allocate(this%uphi_ex(neqn, nlength))
        call this%calc_KG(CData, lna_ini, lna_wkb, neqn, uphi0, n_step, nlength, this%a_ex, this%uphi_ex)


        do ix = 1, nlength
            call CSFconfHubble(this%a_ex(ix), this%uphi_ex(:,ix), this%confH_ex(ix), CData)
            ! write(*,*) this%confH_ex(ix)/this%a_ex(ix)
        enddo

        this%A_wkb = a_wkb**(1.5_dl)*this%uphi_ex(1,nlength)
        this%B_wkb = a_wkb**(0.5_dl)*this%uphi_ex(2,nlength)/mnorm
        this%C_wkb = a_wkb**(1.5_dl)*this%uphi_ex(3,nlength)
        this%D_wkb = a_wkb**(0.5_dl)*this%uphi_ex(4,nlength)/mnorm

    end subroutine ExSolution


    subroutine Exfluid(this,  CData, neqn, N_ex)
        !!! 计算精确流体
        implicit none
        class(CSFBackground) :: this
        class(CAMBdata) :: CData
        integer, intent(in) :: N_ex, neqn
        !!! 辅助变量
        integer :: ix
        real(dl):: mnorm

        !!! 读取数据
        mnorm = CData%CP%DarkMatter%mnorm
        !!! 分配内存
        allocate(this%grho_ex( N_ex))
        allocate(this%dotgrho_ex( N_ex))
        allocate(this%gpres_ex( N_ex))
        allocate(this%dotgpres_ex( N_ex))
        allocate(this%w_ex( N_ex))
        allocate(this%cs2_ex( N_ex))


        call this%get_grho_arr(neqn, N_ex, this%a_ex, this%uphi_ex, mnorm, this%grho_ex)
        call this%get_gdotrho_arr(neqn, N_ex, this%a_ex, this%uphi_ex, this%confH_ex, this%dotgrho_ex)
        call this%get_gpres_arr(neqn, N_ex, this%a_ex, this%uphi_ex, mnorm, this%gpres_ex)
        call this%get_gdotgpres_arr(neqn, N_ex, this%a_ex, this%uphi_ex, mnorm, this%confH_ex, this%dotgpres_ex)
        call this%get_w_arr(this%grho_ex, this%gpres_ex , this%w_ex)
        call this%get_cs2_arr(this%a_ex, neqn, N_ex, this%uphi_ex, mnorm, this%confH_ex, this%cs2_ex)
        
    end subroutine Exfluid


    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!! WKB 近似解
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function WKB_confHubble(CData,a) result(confHubble)
        !!!! 当 a >> a_osc, complex scalar field 的解用 WKB 近似解代替，Hubble 参数也满足 WKB 近似形似
        use DarkEnergyInterface
        implicit none
        class(CAMBdata) :: CData
        real(dl) :: a, h
        real(dl) :: dtauda, grhoa2, grhov_t, grhophia2, omphih2
        real(dl) :: confHubble                                     ! 无量纲共形 Hubble

        !!!! 计算 WKB 近似下的 grhophia2 = 8*pi*G*rho_phi*a^4/c^2.
        h = CData%CP%H0/100._dl
        omphih2 = CData%CP%DarkMatter%omphih2
        grhophia2 = CData%grhocrit*omphih2/h**2 * a 
        
        call CData%CP%DarkEnergy%BackgroundDensityAndPressure(CData%grhov, a, grhov_t)

        !  8*pi*G*rho*a**4.
        grhoa2 = CData%grho_no_de(a) +  grhov_t * a**(2._dl) + grhophia2

        if (grhoa2 <= 0) then
            call GlobalError('Universe stops expanding before today (recollapse not supported)', error_unsupported_params)
            dtauda = 0
        else
            dtauda = sqrt(3._dl / grhoa2)
        end if

        confHubble = 1._dl/(a*dtauda)
        
    end function WKB_confHubble


    function Integrate_WKB_fun(CData,a) result(fval)
        !!!! 用于 WKB 相位因子积分的辅助函数
        class(CAMBdata) :: CData
        real(dl) :: a, confHubble, fval

        confHubble = WKB_confHubble(CData,a)
        fval = CData%CP%DarkMatter%mnorm/confHubble
    end function Integrate_WKB_fun



    subroutine get_WKB_phi(this, CData, a, neqn, uphi)
        !!!!!! 用于计算 WKB 近似下的 complex scalar field 场量的值
        use MathUtils, only : Integrate_Romberg
        implicit none
        class(CSFBackground) :: this
        class(CAMBdata) :: CData
        real(dl), intent(in) :: a
        integer , intent(in) :: neqn
        real(dl) :: S_a, SinS, CosS, tol
        real(dl), intent(out) :: uphi(neqn)  
        real(dl) :: A_wkb, B_wkb, C_wkb, D_wkb

        !!! 判断是否能使用 WKB 近似
        if (a>1) then
            error stop "输入的 a 应该小于 1!"
        endif
        if (a<CData%CP%DarkMatter%a_wkb) then
            error stop "输入的 a < a_wkb, 这不符合 WKB 近似的条件！"
        endif
        !!!!! 利用 Romberg 算法计算相位因子 S_a
        tol = 1e-6_dl
        S_a = Integrate_Romberg(CData, Integrate_WKB_fun, a_wkb , a, tol)

        SinS = sin(S_a)
        CosS = cos(S_a)

        !!!!! 计算 axion 的结果
        A_wkb = this%A_wkb
        B_wkb = this%B_wkb
        C_wkb = this%C_wkb
        D_wkb = this%D_wkb

        uphi(1) = a**(-1.5_dl)*(A_wkb* CosS + B_wkb* SinS)
        uphi(2) = CData%CP%DarkMatter%mnorm*a**(-0.5_dl)*( -A_wkb*SinS + B_wkb* CosS)
        uphi(3) = a**(-1.5_dl)*(C_wkb* CosS + D_wkb* SinS)
        uphi(4) = CData%CP%DarkMatter%mnorm*a**(-0.5_dl)*( -C_wkb*SinS + D_wkb* CosS)
        
    end subroutine get_WKB_phi


    subroutine get_WKB_phi_arr(this, CData, neqn, N_wkb, a_arr, uphi)
        !!!!!! 用于计算 WKB 近似下的 complex scalar field 场量的数组
        implicit none
        class(CSFBackground) :: this
        class(CAMBdata) :: CData
        integer , intent(in) :: neqn, N_wkb
        real(dl), intent(in) :: a_arr(N_wkb)
        real(dl), intent(out) :: uphi(neqn, N_wkb)
        integer :: ix

        do ix = 1, N_wkb
            call this%get_WKB_phi(CData, a_arr(ix), neqn, uphi(:,ix))
        enddo
    end subroutine get_WKB_phi_arr


    subroutine WKBfluid(this, CData , N_wkb)
        !!! 计算 WKB 近似的流体
        use RangeUtils
        implicit none
        class(CSFBackground) :: this
        class(CAMBdata) :: CData
        type(TRanges) :: Ra_wkb
        integer , intent(in) :: N_wkb
        !!! 辅助变量
        integer :: ix
        real(dl) :: a_wkb, a0, omphih2, h

        !!! 设置数据
        h = CData%CP%H0/100._dl
        omphih2 = CData%CP%DarkMatter%omphih2
        a_wkb = CData%CP%DarkMatter%a_wkb
        a0 = 1._dl
        
        !!! 分配内存
        allocate(this%a_wkb_arr(N_wkb))
        allocate(this%confH_wkb(N_wkb))
        allocate(this%grho_wkb(N_wkb))
        allocate(this%dotgrho_wkb(N_wkb))
        allocate(this%gpres_wkb(N_wkb))
        allocate(this%dotgpres_wkb(N_wkb))
        allocate(this%w_wkb(N_wkb))
        allocate(this%cs2_wkb(N_wkb))

        !!!! 得到 a_wkb 数组
        call Ra_wkb%Init()
        call Ra_wkb%Add(a_wkb, a0, N_wkb,.true.)
        call Ra_wkb%GetArray()
        this%a_wkb_arr = Ra_wkb%points(2:N_wkb+1)

        !!! 计算 WKB 近似的 comfHubble
        do ix = 1, N_wkb
            this%confH_wkb(ix) = WKB_confHubble(CData ,this%a_wkb_arr(ix))
        enddo
        !!! 无压流体近似
        this%gpres_wkb = 0._dl
        this%w_wkb = 0._dl
        this%cs2_wkb = 0._dl
        this%grho_wkb = CData%grhocrit*omphih2/h**2/this%a_wkb_arr
        this%dotgrho_wkb = -3._dl*this%grho_wkb*this%confH_wkb
        this%dotgpres_wkb = 0._dl
        
    end subroutine WKBfluid


    !!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!!!
    !!! 估算振荡的性质 Get_OSC
    !!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! 用 LambdaCDM 宇宙模型估算 a_osc 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine get_a_osc(this,  CData)
        implicit none
        class(Get_OSC) :: this
        class(CAMBdata) :: CData
        real(dl) :: a_osc, a_wkb
        real(dl) :: mtoH


        mtoH = CData%CP%DarkMatter%mtoH_osc
        call this%get_mtoH_a( CData, mtoH, a_osc) 
        mtoH = CData%CP%DarkMatter%mtoH_wkb
        call this%get_mtoH_a( CData, mtoH, a_wkb) 
        CData%CP%DarkMatter%a_osc = a_osc
        CData%CP%DarkMatter%a_wkb = a_wkb
        ! write(*,*) a_osc, a_wkb
    
        
    end subroutine get_a_osc

    
    function mtoHfun(this, loga) result(fval)
        !!! 用于估算振荡的尺度因子的 LambdaCDM 宇宙模型， a 为尺度因子
        implicit none
        class(Get_OSC) :: this
        type(CAMBdata) :: CData
        real(dl) :: loga
        real(dl) :: mtoH
        real(dl) :: a, H
        real(dl) :: mnorm   ! 无量纲质量
        real(dl) :: fval    ! 输出函数值



        !!!!!!!! 读取参数
        mtoH = this%mtoH
        CData = this%CData
        mnorm = CData%CP%DarkMatter%mnorm

        ! write(*,*) "mu = ",mu

        
        !!!!!! 计算 Hubble 率
        a = 10.0_dl**loga
        H = WKB_confHubble(CData,a)/a

        ! write(*,*) 'H0 = ', WKB_confHubble(CData,1._dl)*Mpc/1000

        fval =    log(mnorm/H/(mtoH)) 

        ! write(*,*) fval
    
        
    end function mtoHfun



    subroutine get_mtoH_a(this, CData, mtoH, a_obj)
        !!! 从 mtoH 估算对应的 a 的值，采用 Lambda CDM 宇宙模型
        use MathUtils
        implicit none
        class(Get_OSC) :: this
        type(CAMBdata) :: CData
        real(dl), intent(in) :: mtoH          ! 输入的 m 和 H 的比
        real(dl), intent(out) :: a_obj        ! 输出对应 mtoH 的估算的 a
        !!!!!! 辅助变量
        real(dl) :: log_a1, log_a2, fa1, fa2, tol , x0, f0, n   ! 求解区间 (log_a1,log_a2) 及对应函数值 (fa1,fa2)
        integer :: iflag

        this%mtoH = mtoH
        this%CData = CData

        !!!! 计算初始区间和对应函数值
        log_a1 = -10.0_dl
        log_a2 = 0.0_dl
        fa1 = this%mtoHfun( log_a1)
        fa2 = this%mtoHfun( log_a2)
        tol = 1.0e-9_dl

        ! write(*,*) "fa1 = ", fa1
        ! write(*,*) "fa2 = ", fa2

        if (fa1>0.and.fa2>0) then
            do while (fa1>0)
                log_a1 = log_a1 - 5._dl
                fa1 = this%mtoHfun( log_a1)
            end do
        endif

        if (fa1<0.and.fa2<0) then
            a_obj = 1._dl
            return 
        endif

        ! write(*,*) "fa1 = ", fa1
        ! write(*,*) "fa2 = ", fa2

        !!! 用 Brent 方法找零点
        if (fa1*fa2<0._dl) then
            call brentq(this, mtoHfun, log_a1, log_a2 , tol , x0 , f0 , iflag , fa1 , fa2 )
            a_obj = 10._dl**(x0)
        endif
    
        
    end subroutine get_mtoH_a











    !!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc






end module DarkMatterCSF
