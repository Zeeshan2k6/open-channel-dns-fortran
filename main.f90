module params
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  
  ! Global Parameters (Now Variables)
  integer :: nx, ny, nz
  integer :: max_steps  ! Previously nsteps
  integer :: stats_freq, checkpoint_freq
  real(rk) :: dt, re_tau, nu, dpdx
  
  ! Domain Size (Now Variables)
  real(rk) :: lx, ly, lz
  
  ! Restart & Config Flags
  logical :: restart_enabled
  logical :: restart_reset_stats
  logical :: write_geometry_flag
  character(len=64) :: restart_filename
  
  ! Constants
  integer, parameter :: poisson_iters = 80
  
  ! Roughness Parameters (Now Variables)
  real(rk) :: rough_amp
  
  real(rk), parameter :: rough_wx = 2.0_rk
  real(rk), parameter :: rough_wz = 2.0_rk
  
  ! Bead Parameters (Now Variables)
  real(rk) :: bead_radius
  real(rk) :: bead_layer_height
  integer :: bead_nx
  integer :: bead_nz
  
contains
  subroutine read_input()
    integer :: i, ios
    character(len=128) :: line, key, val_str
    
    ! Default Values (Based on D=1 dimensionless scaling)
    ! Assuming H/D ~ 10 for open channel
    nx = 128; ny = 129; nz = 128
    lx = 40.0_rk * acos(-1.0_rk)     ! ~125.66 (40*pi)
    ly = 10.0_rk                     ! H = 10*D
    lz = (40.0_rk / 3.0_rk) * acos(-1.0_rk) ! ~41.88
    
    dt = 0.005_rk                    ! Increased dt for larger domain
    max_steps = 1000
    re_tau = 180.0_rk
    dpdx = -1.0_rk                   ! Adjust driving force if needed for new scale
    stats_freq = 50
    checkpoint_freq = 1000
    restart_enabled = .false.
    restart_reset_stats = .false.
    write_geometry_flag = .true.
    restart_filename = ""
    
    ! Default Bead Geometry (Packed Bed, D=1.0)
    bead_radius = 0.5_rk             ! D = 1.0
    bead_layer_height = 1.2_rk       ! Center of first layer approx D/2 + roughness
    bead_nx = 126
    bead_nz = 42
    
    ! Roughness (Default 0 for smooth beads)
    rough_amp = 0.0_rk
    
    open(unit=10, file="simulation.in", status="old", iostat=ios)
    if (ios /= 0) then
      print *, "WARNING: simulation.in not found! Using defaults."
      nu = 1.0_rk / re_tau
      return
    end if
    
    do
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit
      
      ! Strip comments (everything after '!')
      i = index(line, '!')
      if (i > 0) line = line(1:i-1)
      
      if (len_trim(line) == 0) cycle
      
      i = index(line, '=')
      if (i > 0) then
        key = trim(adjustl(line(1:i-1)))
        val_str = trim(adjustl(line(i+1:)))
        
        ! Convert key to UPPERCASE for robustness
        call to_upper(key)
        
        select case (key)
        ! Grid Points
        case ('NX')
          read(val_str, *) nx
        case ('NY')
          read(val_str, *) ny
        case ('NZ')
          read(val_str, *) nz
        
        ! Domain Size
        case ('LX')
          read(val_str, *) lx
        case ('LY')
          read(val_str, *) ly
        case ('LZ')
          read(val_str, *) lz
          
        ! Time & Physics
        case ('DT')
          read(val_str, *) dt
        case ('MAX_STEPS')
          read(val_str, *) max_steps
        case ('RE_TAU')
          read(val_str, *) re_tau
        case ('DPDX')
          read(val_str, *) dpdx
          
        ! Output
        case ('STATS_INTERVAL')
          read(val_str, *) stats_freq
        case ('CHECKPOINT_INTERVAL')
          read(val_str, *) checkpoint_freq
        case ('WRITE_GEOMETRY')
          read(val_str, *) write_geometry_flag

        ! Restart
        case ('ENABLE_RESTART')
          read(val_str, *) restart_enabled
        case ('RESTART_FILE')
          read(val_str, *) restart_filename  
          if (restart_filename(1:1) == '"') restart_filename = restart_filename(2:len_trim(restart_filename)-1)
        case ('RESET_STATS')
          read(val_str, *) restart_reset_stats
        
        ! Bead Geometry Input
        case ('BEAD_NX')
          read(val_str, *) bead_nx
        case ('BEAD_NZ')
          read(val_str, *) bead_nz
        case ('BEAD_RADIUS')
          read(val_str, *) bead_radius
        case ('BEAD_HEIGHT')
          read(val_str, *) bead_layer_height
        case ('ROUGH_AMP')
          read(val_str, *) rough_amp
        end select
      end if
    end do
    close(10)
    
    nu = 1.0_rk / re_tau
    print *, "Configuration Loaded:"
    print *, "  Grid:", nx, "x", ny, "x", nz
    print *, "  Size:", lx, "x", ly, "x", lz
    print *, "  Beads:", bead_nx, "x", bead_nz, " Radius:", bead_radius
    print *, "  Roughness:", rough_amp
    print *, "  dt:", dt, " Max Steps:", max_steps
  end subroutine read_input

  subroutine to_upper(str)
    character(*), intent(inout) :: str
    integer :: i, ic
    do i = 1, len_trim(str)
      ic = ichar(str(i:i))
      if (ic >= ichar('a') .and. ic <= ichar('z')) then
        str(i:i) = char(ic - 32)
      end if
    end do
  end subroutine to_upper
end module params

module grid
  use params, only: rk, nx, ny, nz, lx, ly, lz
  implicit none
  real(rk), allocatable, dimension(:) :: x
  real(rk), allocatable, dimension(:) :: y
  real(rk), allocatable, dimension(:) :: z
  real(rk) :: dx, dy, dz
contains
  subroutine init_grid()
    integer :: i
    
    if (nx > 0) allocate(x(nx))
    if (ny > 0) allocate(y(ny))
    if (nz > 0) allocate(z(nz))
    
    dx = lx / real(nx, rk)
    dy = ly / real(ny - 1, rk)
    dz = lz / real(nz, rk)
    
    do i = 1, nx
      x(i) = (i - 1) * dx
    end do
    do i = 1, ny
      y(i) = (i - 1) * dy
    end do
    do i = 1, nz
      z(i) = (i - 1) * dz
    end do
  end subroutine init_grid
end module grid

module field
  use params, only: rk, nx, ny, nz
  implicit none
  ! Current Fields
  real(rk), allocatable, dimension(:,:,:) :: u, v, w, p
  real(rk), allocatable, dimension(:,:,:) :: p_new ! For Jacobi Solver
  
  ! Statistics
  real(rk), allocatable, dimension(:,:,:) :: sum_u, sum_v, sum_w, sum_p
  real(rk), allocatable, dimension(:,:,:) :: sum_uu, sum_vv, sum_ww, sum_pp
  real(rk), allocatable, dimension(:,:,:) :: sum_uv, sum_uw, sum_vw
  integer :: stats_count
  
  ! Geometry
  logical, allocatable, dimension(:,:,:) :: is_solid
  
  ! Intermediate Fields (for RHS)
  real(rk), allocatable, dimension(:,:,:) :: rhs_u, rhs_v, rhs_w, div_star
  
contains
  subroutine allocate_fields()
    allocate(u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), p(nx,ny,nz))
    allocate(p_new(nx,ny,nz))
    
    allocate(sum_u(nx,ny,nz), sum_v(nx,ny,nz), sum_w(nx,ny,nz), sum_p(nx,ny,nz))
    allocate(sum_uu(nx,ny,nz), sum_vv(nx,ny,nz), sum_ww(nx,ny,nz), sum_pp(nx,ny,nz))
    allocate(sum_uv(nx,ny,nz), sum_uw(nx,ny,nz), sum_vw(nx,ny,nz))
    
    allocate(is_solid(nx,ny,nz))
    allocate(rhs_u(nx,ny,nz), rhs_v(nx,ny,nz), rhs_w(nx,ny,nz), div_star(nx,ny,nz))
  end subroutine allocate_fields

  subroutine initialize_fields()
    u = 0.0_rk; v = 0.0_rk; w = 0.0_rk; p = 0.0_rk; p_new = 0.0_rk
    
    sum_u = 0.0_rk; sum_v = 0.0_rk; sum_w = 0.0_rk; sum_p = 0.0_rk
    sum_uu = 0.0_rk; sum_vv = 0.0_rk; sum_ww = 0.0_rk; sum_pp = 0.0_rk
    sum_uv = 0.0_rk; sum_uw = 0.0_rk; sum_vw = 0.0_rk
    
    stats_count = 0
    is_solid = .false.
  end subroutine initialize_fields
  
  subroutine accumulate_statistics()
    stats_count = stats_count + 1
    !$OMP PARALLEL DO COLLAPSE(3)
    do k=1,nz
      do j=1,ny
        do i=1,nx
             sum_u(i,j,k) = sum_u(i,j,k) + u(i,j,k)
             sum_v(i,j,k) = sum_v(i,j,k) + v(i,j,k)
             sum_w(i,j,k) = sum_w(i,j,k) + w(i,j,k)
             sum_p(i,j,k) = sum_p(i,j,k) + p(i,j,k)
             
             sum_uu(i,j,k) = sum_uu(i,j,k) + u(i,j,k)**2
             sum_vv(i,j,k) = sum_vv(i,j,k) + v(i,j,k)**2
             sum_ww(i,j,k) = sum_ww(i,j,k) + w(i,j,k)**2
             sum_pp(i,j,k) = sum_pp(i,j,k) + p(i,j,k)**2
             
             sum_uv(i,j,k) = sum_uv(i,j,k) + u(i,j,k)*v(i,j,k)
             sum_uw(i,j,k) = sum_uw(i,j,k) + u(i,j,k)*w(i,j,k)
             sum_vw(i,j,k) = sum_vw(i,j,k) + v(i,j,k)*w(i,j,k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine accumulate_statistics
end module field

! --- [Skipping modules bc, numerics, scalar as they are mostly unchanged for stats logic] ---
! But wait, 'write_time_averaged_stats' is in 'io_utils' which is further down.
! I must target 'module field' first, then 'io_utils' separately because they are far apart.
! This replacement block handles 'module field'.


module bc
  use params, only: rk, nx, ny, nz, rough_amp, rough_wx, rough_wz, bead_radius, bead_layer_height, bead_nx, bead_nz
  use grid, only: x, y, z, dx, dz
  use field, only: u, v, w, p, is_solid, p_new ! Needs p_new for BCs
  implicit none
contains
  real(rk) function roughness_height(ix, iz)
    integer, intent(in) :: ix, iz
    roughness_height = rough_amp * (sin(rough_wx * x(ix)) * sin(rough_wz * z(iz)))
  end function roughness_height

  logical function in_bead(ix, iy, iz)
    integer, intent(in) :: ix, iy, iz
    integer :: ibx, ibz
    real(rk) :: xc, zc, r2
    in_bead = .false.
    if (y(iy) > bead_layer_height) return
    do ibx = 1, bead_nx
      xc = (real(ibx, rk) - 0.5_rk) * (x(nx) + dx) / real(bead_nx, rk)
      do ibz = 1, bead_nz
        zc = (real(ibz, rk) - 0.5_rk) * (z(nz) + dz) / real(bead_nz, rk)
        r2 = (x(ix) - xc)**2 + (z(iz) - zc)**2 + (y(iy) - bead_radius)**2
        if (r2 <= bead_radius**2) then
          in_bead = .true.
          return
        end if
      end do
    end do
  end function in_bead

  subroutine build_immersed_boundary()
    integer :: i, j, k
    real(rk) :: h
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          h = roughness_height(i, k)
          is_solid(i, j, k) = (y(j) <= h) .or. in_bead(i, j, k)
        end do
      end do
    end do
  end subroutine build_immersed_boundary

  subroutine apply_immersed_boundary()
    integer :: i, j, k
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (is_solid(i, j, k)) then
            u(i, j, k) = 0.0_rk
            v(i, j, k) = 0.0_rk
            w(i, j, k) = 0.0_rk
          end if
        end do
      end do
    end do
  end subroutine apply_immersed_boundary

  subroutine apply_velocity_bc()
    integer :: i, k
    
    !$OMP PARALLEL DO PRIVATE(i,k)
    do k = 1, nz
      do i = 1, nx
        ! No-slip at Bottom
        u(i, 1, k) = 0.0_rk
        v(i, 1, k) = 0.0_rk
        w(i, 1, k) = 0.0_rk
        
        ! Free-slip at Top
        u(i, ny, k) = u(i, ny-1, k) ! dU/dy = 0
        v(i, ny, k) = 0.0_rk        ! V = 0
        w(i, ny, k) = w(i, ny-1, k) ! dW/dy = 0
      end do
    end do
    !$OMP END PARALLEL DO
    call apply_immersed_boundary()
  end subroutine apply_velocity_bc

  subroutine apply_pressure_bc()
    integer :: i, k
    do k = 1, nz
      do i = 1, nx
        p(i, 1, k) = p(i, 2, k)
        p(i, ny, k) = p(i, ny - 1, k)
      end do
    end do
  end subroutine apply_pressure_bc
end module bc

module numerics
  use params, only: rk, nx, ny, nz, dt, nu, dpdx, poisson_iters
  use grid, only: dx, dy, dz
  use field, only: u, v, w, p, p_new, rhs_u, rhs_v, rhs_w, div_star, is_solid
  use bc, only: apply_velocity_bc
  implicit none
contains
  integer function ip(i, n)
    integer, intent(in) :: i, n
    ip = i + 1
    if (ip > n) ip = 1
  end function ip

  integer function im(i, n)
    integer, intent(in) :: i, n
    im = i - 1
    if (im < 1) im = n
  end function im

  subroutine compute_tentative()
    integer :: i, j, k
    integer :: ipx, imx, ipz, imz
    real(rk) :: dudx, dudy, dudz
    real(rk) :: dvdx, dvdy, dvdz
    real(rk) :: dwdx, dwdy, dwdz
    real(rk) :: lapu, lapv, lapw

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,ipx,imx,ipz,imz,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,lapu,lapv,lapw) COLLAPSE(2)
    do k = 1, nz
      do j = 2, ny - 1
        do i = 1, nx
          ipx = ip(i, nx)
          imx = im(i, nx)
          ipz = ip(k, nz)
          imz = im(k, nz)

          if (is_solid(i,j,k)) then
             rhs_u(i,j,k) = 0.0_rk
             rhs_v(i,j,k) = 0.0_rk
             rhs_w(i,j,k) = 0.0_rk
             cycle
          endif

          dudx = (u(ipx, j, k) - u(imx, j, k)) / (2.0_rk * dx)
          dudy = (u(i, j + 1, k) - u(i, j - 1, k)) / (2.0_rk * dy)
          dudz = (u(i, j, ipz) - u(i, j, imz)) / (2.0_rk * dz)

          dvdx = (v(ipx, j, k) - v(imx, j, k)) / (2.0_rk * dx)
          dvdy = (v(i, j + 1, k) - v(i, j - 1, k)) / (2.0_rk * dy)
          dvdz = (v(i, j, ipz) - v(i, j, imz)) / (2.0_rk * dz)

          dwdx = (w(ipx, j, k) - w(imx, j, k)) / (2.0_rk * dx)
          dwdy = (w(i, j + 1, k) - w(i, j - 1, k)) / (2.0_rk * dy)
          dwdz = (w(i, j, ipz) - w(i, j, imz)) / (2.0_rk * dz)

          lapu = (u(ipx, j, k) - 2.0_rk * u(i, j, k) + u(imx, j, k)) / (dx * dx)
          lapu = lapu + (u(i, j + 1, k) - 2.0_rk * u(i, j, k) + u(i, j - 1, k)) / (dy * dy)
          lapu = lapu + (u(i, j, ipz) - 2.0_rk * u(i, j, k) + u(i, j, imz)) / (dz * dz)

          lapv = (v(ipx, j, k) - 2.0_rk * v(i, j, k) + v(imx, j, k)) / (dx * dx)
          lapv = lapv + (v(i, j + 1, k) - 2.0_rk * v(i, j, k) + v(i, j - 1, k)) / (dy * dy)
          lapv = lapv + (v(i, j, ipz) - 2.0_rk * v(i, j, k) + v(i, j, imz)) / (dz * dz)

          lapw = (w(ipx, j, k) - 2.0_rk * w(i, j, k) + w(imx, j, k)) / (dx * dx)
          lapw = lapw + (w(i, j + 1, k) - 2.0_rk * w(i, j, k) + w(i, j - 1, k)) / (dy * dy)
          lapw = lapw + (w(i, j, ipz) - 2.0_rk * w(i, j, k) + w(i, j, imz)) / (dz * dz)

          rhs_u(i, j, k) = u(i, j, k) + dt * (-(u(i, j, k) * dudx + v(i, j, k) * dudy + w(i, j, k) * dudz) + nu * lapu - dpdx)
          rhs_v(i, j, k) = v(i, j, k) + dt * (-(u(i, j, k) * dvdx + v(i, j, k) * dvdy + w(i, j, k) * dvdz) + nu * lapv)
          rhs_w(i, j, k) = w(i, j, k) + dt * (-(u(i, j, k) * dwdx + v(i, j, k) * dwdy + w(i, j, k) * dwdz) + nu * lapw)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    
    ! Boundary Updates (Parallel safe if k-loop)
    !$OMP PARALLEL DO PRIVATE(i,k)
    do k=1,nz
       do i=1,nx
          rhs_u(i,1,k) = 0.0_rk; rhs_v(i,1,k)=0.0_rk; rhs_w(i,1,k)=0.0_rk
          rhs_u(i,ny,k) = rhs_u(i,ny-1,k); rhs_v(i,ny,k)=0.0_rk; rhs_w(i,ny,k)=rhs_w(i,ny-1,k)
       end do
    end do
    !$OMP END PARALLEL DO
    
    ! Update Tentative Velocity
    u = rhs_u; v = rhs_v; w = rhs_w
  end subroutine compute_tentative

  subroutine solve_pressure()
    integer :: iter, i, j, k
    integer :: ipx, imx, ipz, imz
    real(rk) :: coeff
    
    coeff = 0.5_rk / (1.0_rk/dx**2 + 1.0_rk/dy**2 + 1.0_rk/dz**2)
    
    ! Compute Divergence (RHS for Poisson)
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,ipx,imx,ipz,imz) COLLAPSE(2)
    do k=1,nz
       do j=2,ny-1
          ipz = ip(k, nz)
          imz = im(k, nz)
          do i=1,nx
             ipx = ip(i, nx)
             imx = im(i, nx)
             
             div_star(i,j,k) = ( (u(ipx,j,k)-u(imx,j,k))/(2*dx) + &
                                 (v(i,j+1,k)-v(i,j-1,k))/(2*dy) + &
                                 (w(i,j,ipz)-w(i,j,imz))/(2*dz) ) / dt
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    ! Solve Poisson Equation (Jacobi Iteration)
    do iter = 1, poisson_iters
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,ipx,imx,ipz,imz) COLLAPSE(2)
       do k=1,nz
          do j=2,ny-1
             ipz = ip(k, nz)
             imz = im(k, nz)
             do i=1,nx
                ipx = ip(i, nx)
                imx = im(i, nx)
                
                p_new(i,j,k) = coeff * ( &
                             (p(ipx,j,k)+p(imx,j,k))/dx**2 + &
                             (p(i,j+1,k)+p(i,j-1,k))/dy**2 + &
                             (p(i,j,ipz)+p(i,j,imz))/dz**2 - div_star(i,j,k) )
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       
       ! Apply Pressure BCs to p_new
       !$OMP PARALLEL DO PRIVATE(i,k)
       do k=1,nz
          do i=1,nx
             p_new(i,1,k) = p_new(i,2,k)
             p_new(i,ny,k) = p_new(i,ny-1,k)
          end do
       end do
       !$OMP END PARALLEL DO
       
       ! Update
       p = p_new
    end do
  end subroutine solve_pressure

  subroutine project_velocity()
    integer :: i, j, k
    integer :: ipx, imx, ipz, imz
    
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,ipx,imx,ipz,imz) COLLAPSE(2)
    do k=1,nz
       do j=2,ny-1
          ipz = ip(k, nz)
          imz = im(k, nz)
          do i=1,nx
             ipx = ip(i, nx)
             imx = im(i, nx)
             
             if (is_solid(i,j,k)) cycle
             
             u(i,j,k) = u(i,j,k) - dt * (p(ipx,j,k)-p(imx,j,k))/(2.0_rk*dx)
             v(i,j,k) = v(i,j,k) - dt * (p(i,j+1,k)-p(i,j-1,k))/(2.0_rk*dy)
             w(i,j,k) = w(i,j,k) - dt * (p(i,j,ipz)-p(i,j,imz))/(2.0_rk*dz)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    call apply_velocity_bc()
  end subroutine project_velocity
end module numerics

module scalar
  use params, only: rk, nx, ny, nz, dt
  use grid, only: dx, dy, dz
  use field, only: u, v, w, is_solid
  implicit none
  
  real(rk), allocatable :: dye(:,:,:), new_dye(:,:,:)
  
contains
  subroutine init_scalar()
    allocate(dye(nx, ny, nz), new_dye(nx, ny, nz))
    dye = 0.0_rk; new_dye = 0.0_rk
  end subroutine init_scalar
  
  subroutine inject_dye()
    integer :: k
    ! Continuous Injection Source: Line source at x=2.0 (approx index 10-20), near bed
    ! Adjust indices based on your LX. If Lx=24, dx~0.06. x=2.0 is i~32.
    integer :: i_src, j_src
    
    i_src = max(1, min(nx-2, int(2.0_rk / dx)))
    ! Inject in middle of fluid column (above beads)
    real(rk) :: y_src
    
    y_src = bead_layer_height + 0.5_rk * (ly - bead_layer_height)
    
    i_src = max(1, min(nx-2, int(2.0_rk / dx)))
    j_src = max(2, min(ny-2, int(y_src / dy)))
    
    ! Inject in a cross-stream line (z=2 to z=10)
    !$OMP PARALLEL DO PRIVATE(k)
    do k = 1, nz
       if (k > nz/4 .and. k < 3*nz/4) then
         dye(i_src:i_src+2, j_src:j_src+4, k) = 1.0_rk
       endif
    end do
    !$OMP END PARALLEL DO
  end subroutine inject_dye
  
  subroutine update_scalar()
    integer :: i, j, k
    real(rk) :: flux_x, flux_y, flux_z
    
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,flux_x,flux_y,flux_z) COLLAPSE(2)
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
           if (is_solid(i,j,k)) cycle
           
           ! Streamwise (u > 0 assumed)
           flux_x = u(i,j,k) * (dye(i,j,k) - dye(i-1,j,k)) / dx
           
           ! Vertical (simple central or upwind depending on sign)
           if (v(i,j,k) > 0.0_rk) then
             flux_y = v(i,j,k) * (dye(i,j,k) - dye(i,j-1,k)) / dy
           else
             flux_y = v(i,j,k) * (dye(i,j+1,k) - dye(i,j,k)) / dy
           endif
           
           ! Spanwise
           if (w(i,j,k) > 0.0_rk) then
             flux_z = w(i,j,k) * (dye(i,j,k) - dye(i,j,k-1)) / dz
           else
             flux_z = w(i,j,k) * (dye(i,j,k+1) - dye(i,j,k)) / dz
           endif
           
           new_dye(i,j,k) = dye(i,j,k) - dt * (flux_x + flux_y + flux_z)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    
    dye = new_dye
    
    ! --- Boundary Conditions ---
    ! Inlet (x=1): Clean water
    !$OMP PARALLEL DO PRIVATE(k)
    do k=1,nz
         dye(1,:,k) = 0.0_rk
    end do
    !$OMP END PARALLEL DO
    
    ! Outlet (x=nx): Zero Gradient (Outflow)
    !$OMP PARALLEL DO PRIVATE(k)
    do k=1,nz
         dye(nx,:,k) = dye(nx-1,:,k)
    end do
    !$OMP END PARALLEL DO
    
    ! Walls (y=1, y=ny): No flux (copy neighbor)
    !$OMP PARALLEL DO PRIVATE(i,k)
    do k=1,nz
       do i=1,nx
          dye(i,1,k) = dye(i,2,k)
          dye(i,ny,k) = dye(i,ny-1,k)
       end do
    end do
    !$OMP END PARALLEL DO
    
    ! Periodic Z
    dye(:, :, 1) = dye(:, :, nz-1)
    dye(:, :, nz) = dye(:, :, 2)
    
  end subroutine update_scalar
end module scalar

module io_utils
  use params, only: rk, nx, ny, nz, dt, bead_nx, bead_nz, lx, lz, bead_layer_height
  use grid, only: x, y, z
  use field, only: u, v, w, p, sum_u, sum_v, sum_w, sum_uu, sum_vv, sum_ww, &
                   sum_p, sum_pp, sum_uv, sum_uw, sum_vw, stats_count, is_solid
  use scalar, only: dye
  implicit none
contains
  subroutine write_flow_vtk(step)
    integer, intent(in) :: step
    integer :: i, j, k
    character(len=64) :: filename
    
    write(filename, '("dump/flow_", I0.6, ".vtk")') step
    open(unit=40, file=filename, status='replace')
    
    write(40, '(A)') '# vtk DataFile Version 3.0'
    write(40, '(A)') 'DNS Flow Field'
    write(40, '(A)') 'ASCII'
    write(40, '(A)') 'DATASET STRUCTURED_POINTS'
    write(40, '(A, 3(I5, 1x))') 'DIMENSIONS', nx, ny, nz
    write(40, '(A)') 'ORIGIN 0.0 0.0 0.0'
    ! Assuming uniform spacing for VTK simplicity, or use RECTILINEAR_GRID for stretched
    write(40, '(A, 3(E12.5, 1x))') 'SPACING', x(2)-x(1), y(2)-y(1), z(2)-z(1)
    
    write(40, '(A, I0)') 'POINT_DATA ', nx*ny*nz
    
    ! Write Velocity Vectors
    write(40, '(A)') 'VECTORS velocity float'
    do k=1,nz
      do j=1,ny
        do i=1,nx
          write(40, '(3(E12.5, 1x))') u(i,j,k), v(i,j,k), w(i,j,k)
        end do
      end do
    end do
    
    ! Write Dye Scalar
    write(40, '(A)') 'SCALARS dye float 1'
    write(40, '(A)') 'LOOKUP_TABLE default'
    do k=1,nz
      do j=1,ny
        do i=1,nx
          write(40, '(E12.5)') dye(i,j,k)
        end do
      end do
    end do
    
    close(40)
    print *, "Flow viz written: ", trim(filename)
  end subroutine write_flow_vtk

  subroutine write_time_averaged_stats()
    integer :: j
    real(rk) :: m_u, m_v, m_w, m_p
    real(rk) :: m_uu, m_vv, m_ww, m_pp
    real(rk) :: m_uv, m_uw, m_vw
    real(rk) :: u_rms, v_rms, w_rms, p_rms
    real(rk) :: uv, uw, vw ! Reynolds Stresses <u'v'> etc
    real(rk) :: factor
    character(len=200) :: header

    if (stats_count == 0) return

    factor = 1.0_rk / real(stats_count, rk)

    open(unit=10, file="result/time_averaged_stats.dat", status="replace")
    
    header = '      y           <U>          <V>          <W>          <P>        ' // &
             'u_rms        v_rms        w_rms        p_rms       ' // &
             '<u''v''>       <u''w''>       <v''w''>'
    write(10, '(A)') trim(header)
    
    do j = 1, ny
      m_u = sum(sum_u(:, j, :)) * factor / real(nx * nz, rk)
      m_v = sum(sum_v(:, j, :)) * factor / real(nx * nz, rk)
      m_w = sum(sum_w(:, j, :)) * factor / real(nx * nz, rk)
      m_p = sum(sum_p(:, j, :)) * factor / real(nx * nz, rk)

      m_uu = sum(sum_uu(:, j, :)) * factor / real(nx * nz, rk)
      m_vv = sum(sum_vv(:, j, :)) * factor / real(nx * nz, rk)
      m_ww = sum(sum_ww(:, j, :)) * factor / real(nx * nz, rk)
      m_pp = sum(sum_pp(:, j, :)) * factor / real(nx * nz, rk)
      
      m_uv = sum(sum_uv(:, j, :)) * factor / real(nx * nz, rk)
      m_uw = sum(sum_uw(:, j, :)) * factor / real(nx * nz, rk)
      m_vw = sum(sum_vw(:, j, :)) * factor / real(nx * nz, rk)

      u_rms = sqrt(max(0.0_rk, m_uu - m_u**2))
      v_rms = sqrt(max(0.0_rk, m_vv - m_v**2))
      w_rms = sqrt(max(0.0_rk, m_ww - m_w**2))
      p_rms = sqrt(max(0.0_rk, m_pp - m_p**2))
      
      ! Reynolds Stresses <u'v'> = <uv> - <u><v>
      ! Note: This formulation is correct for time-averaging
      uv = m_uv - m_u*m_v
      uw = m_uw - m_u*m_w
      vw = m_vw - m_v*m_w

      write(10, '(12(1x,E12.5))') y(j), m_u, m_v, m_w, m_p, &
                                  u_rms, v_rms, w_rms, p_rms, &
                                  uv, uw, vw
    end do
    close(10)
  end subroutine write_time_averaged_stats

  subroutine write_ideal_beads()
    integer :: i, k
    real(rk) :: xc, zc, space_x, space_z
    character(len=64) :: filename
    
    filename = "dump/beads_ideal.vtk"
    open(unit=21, file=filename, status="replace")
    
    write(21, '(A)') '# vtk DataFile Version 3.0'
    write(21, '(A)') 'Ideal Bead Centers'
    write(21, '(A)') 'ASCII'
    write(21, '(A)') 'DATASET POLYDATA'
    write(21, '(A, I0, A)') 'POINTS ', bead_nx * bead_nz, ' float'
    
    space_x = lx / real(bead_nx, rk)
    space_z = lz / real(bead_nz, rk)
    
    do k = 1, bead_nz
      do i = 1, bead_nx
        ! Calculate Center (Assuming simple grid packing for now)
        xc = (real(i, rk) - 0.5_rk) * space_x
        zc = (real(k, rk) - 0.5_rk) * space_z
        
        write(21, '(3(1x,E12.5))') xc, bead_layer_height, zc
      end do
    end do
    
    close(21)
    print *, "Ideal Bead Centers written to ", trim(filename)
    print *, "--> Open in ParaView and use 'Glyph' Filter (Sphere) to visualize."
  end subroutine write_ideal_beads

  subroutine write_geometry()
    integer :: i, j, k, n_solid
    character(len=64) :: filename
    
    filename = "dump/geometry.vtk"
    open(unit=20, file=filename, status="replace")
    
    ! 1. Count Solid Points
    n_solid = 0
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (is_solid(i, j, k)) n_solid = n_solid + 1
        end do
      end do
    end do
    
    ! 2. Write VTK Header
    write(20, '(A)') '# vtk DataFile Version 3.0'
    write(20, '(A)') 'Immersed Boundary Geometry'
    write(20, '(A)') 'ASCII'
    write(20, '(A)') 'DATASET POLYDATA'
    write(20, '(A, I0, A)') 'POINTS ', n_solid, ' float'
    
    ! 3. Write Points
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (is_solid(i, j, k)) then
            write(20, '(3(1x,E12.5))') x(i), y(j), z(k)
          end if
        end do
      end do
    end do
    
    close(20)
    print *, "Geometry point cloud written to ", trim(filename)
  end subroutine write_geometry

  subroutine write_checkpoint(step)
    integer, intent(in) :: step
    character(len=64) :: filename
    
    ! Ensure output directory exists (handled by Makefile usually)
    ! Filename format: restart/checkpoint_010000.bin
    write(filename, '("restart/checkpoint_", I0.6, ".bin")') step
    open(unit=30, file=filename, form='unformatted', status='replace')
    
    ! 1. Metadata
    write(30) step, dt
    
    ! 2. Instantaneous Fields
    write(30) u, v, w, p
    
    ! 3. Statistics Accumulators
    write(30) stats_count
    write(30) sum_u, sum_v, sum_w, sum_p
    write(30) sum_uu, sum_vv, sum_ww, sum_pp
    write(30) sum_uv, sum_uw, sum_vw
    
    close(30)
    print *, "Checkpoint written: ", trim(filename)
  end subroutine write_checkpoint

  subroutine read_checkpoint(filename, start_step)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: start_step
    real(rk) :: saved_dt
    
    open(unit=31, file=filename, form='unformatted', status='old')
    
    ! 1. Metadata
    read(31) start_step, saved_dt
    
    ! Check Dt consistency
    if (abs(saved_dt - dt) > 1.0e-10) then
       print *, "WARNING: Time step changed from checkpoint!"
    end if
    
    ! 2. Instantaneous Fields
    read(31) u, v, w, p
    
    ! 3. Statistics Accumulators
    read(31) stats_count
    read(31) sum_u, sum_v, sum_w
    read(31) sum_uu, sum_vv, sum_ww
    
    close(31)
    print *, "Checkpoint loaded: ", trim(filename), " starting at step ", start_step
  end subroutine read_checkpoint

  function itoa(i) result(str)
    integer, intent(in) :: i
    character(len=32) :: str
    write(str, '(I0)') i
  end function itoa
end module io_utils

program open_channel_dns
  use params, only: rk, dt, read_input, restart_enabled, restart_filename, &
                    stats_freq, checkpoint_freq, max_steps, restart_reset_stats, write_geometry_flag
  use grid, only: init_grid, dx, dy, dz
  use field, only: allocate_fields, initialize_fields, accumulate_statistics, u, &
                   sum_u, sum_v, sum_w, sum_uu, sum_vv, sum_ww, stats_count
  use bc, only: apply_velocity_bc, build_immersed_boundary
  use numerics, only: compute_tentative, solve_pressure, project_velocity
  use io_utils, only: write_time_averaged_stats, write_geometry, write_checkpoint, read_checkpoint, &
                      write_flow_vtk, write_ideal_beads
  use scalar, only: init_scalar, inject_dye, update_scalar
  use ieee_arithmetic
  implicit none
  integer :: step
  integer :: start_step
  
  real(rk) :: u_max, cfl

  call read_input()
  call init_grid()
  call allocate_fields()
  call init_scalar() ! Allocate Dye Field
  
  ! Check for Restart
  if (restart_enabled) then
    call read_checkpoint(trim(restart_filename), start_step)
    
    if (restart_reset_stats) then
       print *, "RESETTING STATISTICS as requested in input file."
       sum_u = 0.0_rk; sum_v = 0.0_rk; sum_w = 0.0_rk
       sum_uu = 0.0_rk; sum_vv = 0.0_rk; sum_ww = 0.0_rk
       stats_count = 0
    endif
    
    start_step = start_step + 1
  else
    print *, "Starting FRESH simulation."
    call initialize_fields()
    start_step = 1
  endif

  call build_immersed_boundary()
  
  if (write_geometry_flag) then
     call write_geometry()     ! Point Cloud (Rough/Voxel)
     call write_ideal_beads()  ! Smooth Spheres (Visualization)
  endif
  
  call apply_velocity_bc()

  print *, "Step   Max_U      CFL"
  flush(6)
  
  do step = start_step, max_steps
    
    ! --- Scalar Transport (Dye) ---
    call inject_dye()
    call update_scalar()
    
    ! --- Fluid Solver ---
    call compute_tentative()
    call solve_pressure()
    call project_velocity()

    u_max = maxval(abs(u))
    cfl = u_max * dt / min(dx, min(dy, dz))
    
    ! Console Output
    if (mod(step, 100) == 0) then
      print '(I6, 1x, E12.5, 1x, F8.4)', step, u_max, cfl
      flush(6)
    end if

    ! Stability Check
    if (cfl > 1.5 .or. ieee_is_nan(u_max)) then
      print *, "FATAL ERROR: Simulation unstable (CFL > 1.5 or NaN)"
      print *, "Step:", step, " Max U:", u_max, " CFL:", cfl
      stop
    end if

    ! Statistics
    if (restart_enabled .or. step > 100) then
      call accumulate_statistics()
    end if

    if (mod(step, stats_freq) == 0) then
      call write_time_averaged_stats()
    end if
    
    ! Checkpoint & 3D Visualization
    if (mod(step, checkpoint_freq) == 0) then
      call write_checkpoint(step)
      call write_flow_vtk(step)  ! Write full 3D flow viz
    end if
  end do
  
  ! Final checkpoint
  if (mod(max_steps, checkpoint_freq) /= 0) then
    call write_checkpoint(max_steps)
    call write_flow_vtk(max_steps)
  end if
  
end program open_channel_dns
