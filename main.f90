module params
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, parameter :: nx = 32
  integer, parameter :: ny = 33
  integer, parameter :: nz = 32
  real(rk), parameter :: lx = 4.0_rk * acos(-1.0_rk)
  real(rk), parameter :: ly = 1.0_rk
  real(rk), parameter :: lz = (4.0_rk / 3.0_rk) * acos(-1.0_rk)
  real(rk), parameter :: re_tau = 180.0_rk
  real(rk), parameter :: nu = 1.0_rk / re_tau
  real(rk), parameter :: dt = 5.0e-4_rk
  integer, parameter :: nsteps = 10000
  integer, parameter :: poisson_iters = 80
  real(rk), parameter :: dpdx = -1.0_rk
  real(rk), parameter :: rough_amp = 0.02_rk
  real(rk), parameter :: rough_wx = 2.0_rk
  real(rk), parameter :: rough_wz = 2.0_rk
end module params

module grid
  use params, only: rk, nx, ny, nz, lx, ly, lz
  implicit none
  real(rk), dimension(nx) :: x
  real(rk), dimension(ny) :: y
  real(rk), dimension(nz) :: z
  real(rk) :: dx, dy, dz
contains
  subroutine init_grid()
    integer :: i
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
  real(rk), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)
  real(rk), allocatable :: p(:,:,:)
  real(rk), allocatable :: ut(:,:,:), vt(:,:,:), wt(:,:,:)
  real(rk), allocatable :: rhs(:,:,:)
  
  ! Statistics arrays
  real(rk), allocatable :: sum_u(:,:,:), sum_v(:,:,:), sum_w(:,:,:)
  real(rk), allocatable :: sum_uu(:,:,:), sum_vv(:,:,:), sum_ww(:,:,:)
  integer :: stats_count

contains
  subroutine allocate_fields()
    allocate(u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz))
    allocate(p(nx, ny, nz))
    allocate(ut(nx, ny, nz), vt(nx, ny, nz), wt(nx, ny, nz))
    allocate(rhs(nx, ny, nz))
    
    ! Allocate statistics arrays
    allocate(sum_u(nx, ny, nz), sum_v(nx, ny, nz), sum_w(nx, ny, nz))
    allocate(sum_uu(nx, ny, nz), sum_vv(nx, ny, nz), sum_ww(nx, ny, nz))
  end subroutine allocate_fields

  subroutine initialize_fields()
    u = 0.0_rk
    v = 0.0_rk
    w = 0.0_rk
    p = 0.0_rk
    
    ! Initialize statistics
    sum_u = 0.0_rk; sum_v = 0.0_rk; sum_w = 0.0_rk
    sum_uu = 0.0_rk; sum_vv = 0.0_rk; sum_ww = 0.0_rk
    stats_count = 0
  end subroutine initialize_fields
  
  subroutine accumulate_statistics()
    stats_count = stats_count + 1
    sum_u = sum_u + u
    sum_v = sum_v + v
    sum_w = sum_w + w
    sum_uu = sum_uu + u**2
    sum_vv = sum_vv + v**2
    sum_ww = sum_ww + w**2
  end subroutine accumulate_statistics
end module field

module bc
  use params, only: rk, nx, ny, nz, rough_amp, rough_wx, rough_wz
  use grid, only: x, y, z
  use field, only: u, v, w, p
  implicit none
contains
  real(rk) function roughness_height(ix, iz)
    integer, intent(in) :: ix, iz
    roughness_height = rough_amp * (sin(rough_wx * x(ix)) * sin(rough_wz * z(iz)))
  end function roughness_height

  subroutine apply_velocity_bc()
    integer :: i, k
    real(rk) :: h
    do k = 1, nz
      do i = 1, nx
        h = roughness_height(i, k)
        if (y(1) <= h) then
          u(i, 1, k) = 0.0_rk
          v(i, 1, k) = 0.0_rk
          w(i, 1, k) = 0.0_rk
        else
          u(i, 1, k) = 0.0_rk
          v(i, 1, k) = 0.0_rk
          w(i, 1, k) = 0.0_rk
        end if
      end do
    end do

    do k = 1, nz
      do i = 1, nx
        v(i, ny, k) = 0.0_rk
        u(i, ny, k) = u(i, ny - 1, k)
        w(i, ny, k) = w(i, ny - 1, k)
      end do
    end do
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
  use field, only: u, v, w, p, ut, vt, wt, rhs
  use bc, only: apply_velocity_bc, apply_pressure_bc
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

    do k = 1, nz
      ipz = ip(k, nz)
      imz = im(k, nz)
      do j = 2, ny - 1
        do i = 1, nx
          ipx = ip(i, nx)
          imx = im(i, nx)

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

          ut(i, j, k) = u(i, j, k) + dt * (-(u(i, j, k) * dudx + v(i, j, k) * dudy + w(i, j, k) * dudz) + nu * lapu + dpdx)
          vt(i, j, k) = v(i, j, k) + dt * (-(u(i, j, k) * dvdx + v(i, j, k) * dvdy + w(i, j, k) * dvdz) + nu * lapv)
          wt(i, j, k) = w(i, j, k) + dt * (-(u(i, j, k) * dwdx + v(i, j, k) * dwdy + w(i, j, k) * dwdz) + nu * lapw)
        end do
      end do
    end do
  end subroutine compute_tentative

  subroutine solve_pressure()
    integer :: iter, i, j, k
    integer :: ipx, imx, ipz, imz
    real(rk) :: coef
    coef = 2.0_rk / (dx * dx) + 2.0_rk / (dy * dy) + 2.0_rk / (dz * dz)

    do k = 1, nz
      ipz = ip(k, nz)
      imz = im(k, nz)
      do j = 2, ny - 1
        do i = 1, nx
          ipx = ip(i, nx)
          imx = im(i, nx)
          rhs(i, j, k) = ((ut(ipx, j, k) - ut(imx, j, k)) / (2.0_rk * dx) + &
                          (vt(i, j + 1, k) - vt(i, j - 1, k)) / (2.0_rk * dy) + &
                          (wt(i, j, ipz) - wt(i, j, imz)) / (2.0_rk * dz)) / dt
        end do
      end do
    end do

    do iter = 1, poisson_iters
      do k = 1, nz
        ipz = ip(k, nz)
        imz = im(k, nz)
        do j = 2, ny - 1
          do i = 1, nx
            ipx = ip(i, nx)
            imx = im(i, nx)
            p(i, j, k) = ((p(ipx, j, k) + p(imx, j, k)) / (dx * dx) + &
                          (p(i, j + 1, k) + p(i, j - 1, k)) / (dy * dy) + &
                          (p(i, j, ipz) + p(i, j, imz)) / (dz * dz) - rhs(i, j, k)) / coef
          end do
        end do
      end do
      call apply_pressure_bc()
    end do
  end subroutine solve_pressure

  subroutine project_velocity()
    integer :: i, j, k
    integer :: ipx, imx, ipz, imz
    do k = 1, nz
      ipz = ip(k, nz)
      imz = im(k, nz)
      do j = 2, ny - 1
        do i = 1, nx
          ipx = ip(i, nx)
          imx = im(i, nx)
          u(i, j, k) = ut(i, j, k) - dt * (p(ipx, j, k) - p(imx, j, k)) / (2.0_rk * dx)
          v(i, j, k) = vt(i, j, k) - dt * (p(i, j + 1, k) - p(i, j - 1, k)) / (2.0_rk * dy)
          w(i, j, k) = wt(i, j, k) - dt * (p(i, j, ipz) - p(i, j, imz)) / (2.0_rk * dz)
        end do
      end do
    end do
    call apply_velocity_bc()
  end subroutine project_velocity
end module numerics

module io_utils
  use params, only: rk, nx, ny, nz
  use grid, only: y
  use field, only: u, v, w, sum_u, sum_v, sum_w, sum_uu, sum_vv, sum_ww, stats_count
  implicit none
contains
  subroutine write_time_averaged_stats()
    integer :: j
    real(rk) :: mean_u, mean_v, mean_w
    real(rk) :: mean_uu, mean_vv, mean_ww
    real(rk) :: u_rms, v_rms, w_rms
    real(rk) :: factor
    
    if (stats_count == 0) return
    
    factor = 1.0_rk / real(stats_count, rk)
    
    open(unit=10, file="output/time_averaged_stats.dat", status="replace")
    write(10, '(A)') '      y           <U>          <V>          <W>        u_rms        v_rms        w_rms'
    do j = 1, ny
      ! Calculate planar averages of the TIME-AVERAGED fields
      mean_u = sum(sum_u(:, j, :)) * factor / real(nx * nz, rk)
      mean_v = sum(sum_v(:, j, :)) * factor / real(nx * nz, rk)
      mean_w = sum(sum_w(:, j, :)) * factor / real(nx * nz, rk)
      
      mean_uu = sum(sum_uu(:, j, :)) * factor / real(nx * nz, rk)
      mean_vv = sum(sum_vv(:, j, :)) * factor / real(nx * nz, rk)
      mean_ww = sum(sum_ww(:, j, :)) * factor / real(nx * nz, rk)

      ! Calculate RMS: sqrt(<u^2> - <u>^2)
      u_rms = sqrt(max(0.0_rk, mean_uu - mean_u**2))
      v_rms = sqrt(max(0.0_rk, mean_vv - mean_v**2))
      w_rms = sqrt(max(0.0_rk, mean_ww - mean_w**2))

      write(10, '(7(1x,F12.6))') y(j), mean_u, mean_v, mean_w, u_rms, v_rms, w_rms
    end do
    close(10)
  end subroutine write_time_averaged_stats

  function itoa(i) result(str)
    integer, intent(in) :: i
    character(len=32) :: str
    write(str, '(I0)') i
  end function itoa
end module io_utils

program open_channel_dns
  use params, only: nsteps
  use grid, only: init_grid
  use field, only: allocate_fields, initialize_fields, accumulate_statistics
  use bc, only: apply_velocity_bc
  use numerics, only: compute_tentative, solve_pressure, project_velocity
  use io_utils, only: write_time_averaged_stats
  implicit none
  integer :: step

  call init_grid()
  call allocate_fields()
  call initialize_fields()
  call apply_velocity_bc()

  do step = 1, nsteps
    call compute_tentative()
    call solve_pressure()
    call project_velocity()
    
    ! Accumulate statistics starting after initial transient (e.g., step 100)
    if (step > 100) then
      call accumulate_statistics()
    end if
    
    if (mod(step, 100) == 0) then
      call write_time_averaged_stats()
    end if
  end do
end program open_channel_dns