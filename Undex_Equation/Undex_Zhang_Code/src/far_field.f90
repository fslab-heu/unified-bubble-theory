    module far_field
        use DG, only: DG_mesh
        implicit none
        type(DG_mesh) :: far_DG ! nvar = 2, 1 for p-p0, 2 for vel
        real :: k1, k2, k3,  pm
        ! eos state for water at p0 with the same entropy as the shock front
        real :: rho0s, c0s, ss
        

    contains

    subroutine ini_far_field()
        use model
        use near_field, only: near_DG
        use polynomial
        use DG
        implicit none
        integer :: i,j,k, m
        real :: rho, vel , p, xg, xs
        real :: basis(npol)
        real :: lambda, c, v
        real :: rg, rhs(npol)
        real :: ein

        
        ! find L
        do i=near_DG%ncell, 1, -1
            rho = near_DG%unk(1,1,i)
            vel = near_DG%unk(1,2,i)/rho
            ein = near_DG%unk(1,3,i) - 0.5*rho*vel**2
            p = water%press(rho, ein)
            c = water%sound_spd(rho, p)
            v = vel + c - dRs
            if(v<0)then
                exit
            endif
        end do
        lambda = 5*(1.0-near_DG%node(i))*L

        ! initialize far field mesh
        call far_DG%init_DG(2, 50)
        ! project mesh data from nearfield to farfield
        do i = 1, far_DG%ncell
            rhs = 0
            do j = 1, ngauss
                ! physical position
                rg = (0.5*far_DG%cell_length(i)*(gauss_quad%pos(j) + 1.0) + &
                    far_DG%node(i))*lambda + Rs - lambda
                ! position in near_DG
                xg = (rg - Rb)/L
                ! find the cell in near_DG
                do k = 1, near_DG%ncell
                    if((near_DG%node(k)<=xg).and.(near_DG%node(k+1)>=xg)) then
                        exit
                    end if
                end do
                if(k>near_DG%ncell) k = near_DG%ncell
                ! zetag in [-1,1]
                xs = (xg - near_DG%node(k))/near_DG%cell_length(k)
                do m = 1, npol
                    basis(m) = leg_ply(m-1)%evaluate(2.0*xs-1.0)
                end do
                rho = dot_product((near_DG%unk(:,1,k)), basis)
                vel = dot_product((near_DG%unk(:,2,k)), basis)/rho
                p = water%press(rho, dot_product((near_DG%unk(:,3,k)), basis)-0.5*rho*vel**2) - p0
                rhs = rhs + gauss_quad%weight(j)*far_DG%cell_length(i)/2.0*p*basis_gauss(:,j)
            end do
            far_DG%unk(:,1,i) = matmul(far_DG%imass(:,:,i), rhs)
        end do
        L = lambda
        dL = 0
        call update_far_field(far_DG)
        call press2vel()
        ! output far field files
        !!open(101, file='output/far_field_grid.dat', status='replace', action='write', &
        !!    form='formatted')
        !!do i=1, far_DG%ncell
        !!    write(101,'(2E20.10)', advance = 'no') far_DG%node(i),  far_DG%node(i+1)
        !!end do
        !!close(101)
        !!open(101, file='output/far_field.dat', status='replace', action='write', &
        !!    form='formatted')
        
        call output_far_field(far_DG)
        return
    end subroutine ini_far_field

    subroutine press2vel()
        use model, only: us, L, Rs
        use DG, only: gauss_quad, ibasis_gauss, npol, ngauss, basis_gauss, basis_face
        implicit none
        ! convert pressure to velocity based on ESA method
        ! local
        integer :: i, j
        real :: ip(ngauss), ujump, unk(npol)
        real :: rhs(npol), rg(ngauss)

        far_DG%unk(:,2,:) = 0.0
        far_DG%unk(:,2,:) = far_DG%unk(:,1,:) / (rho0s*c0s)
        far_DG%unk(1,2,:) = far_DG%unk(1,2,:) + us - pm/(rho0s*c0s)
        ujump = 0.0
        ! calculate the integration
        do i=far_DG%ncell, 1, -1
            rg = 0.5*far_DG%cell_length(i)*(gauss_quad%pos + 1.0) + &
                far_DG%node(i)
            rg = rg*L + Rs - L
            ! integration of p
            unk = far_DG%unk(:,1,i)
            unk(1) = unk(1) - pm
            ip = (matmul(transpose(ibasis_gauss), unk) + ujump)&
                /rg*far_DG%cell_length(i)/2.0*L
            rhs = 0
            do j=1, npol
                rhs(j) = dot_product(ip, basis_gauss(j,:)*gauss_quad%weight)*&
                            far_DG%cell_length(i)/2.0
            enddo
            unk = matmul(far_DG%imass(:,:,i), rhs)
            far_DG%unk(:,2,i) = far_DG%unk(:,2,i) - unk/rho0s/c0s
            ujump = ujump + (far_DG%unk(1,1,i)-pm)*&
                    far_DG%cell_center(i)/2.0*L
        enddo
        ! debug, output vel_esa
        !!open(105, file='output/vel_esa.dat', status='replace', action='write', &
        !!    form='formatted')
        !!do i=1, far_DG%ncell
        !!    write(105,'(2E20.10)', advance = 'no') &
        !!        dot_product(far_DG%unk(:,2,i), basis_face(:,1)), &
        !!        dot_product(far_DG%unk(:,2,i), basis_face(:,2))
        !!end do
        !!close(105)

    end subroutine
    


    subroutine update_far_field(mesh)
        use model
        use DG
        implicit none
        type(DG_mesh), intent(inout) :: mesh
        ! local
        pm = dot_product(mesh%unk(:,1,mesh%ncell), basis_face(:,2))! + p0
        ps = pm + p0
        call update_shock_front()

        ss = water%entropy(ps, rhos)
        rho0s = water%dens_from_entropy(ss, p0)
        c0s = water%sound_spd(rho0s, p0)

        k1 = us + c0 - dRs - &
            pm/(rho0s*c0s)*(1.0-(1+water%gamma)/4.0*pm/(rho0s*c0s**2))
        k2 = 1/(rho0s*c0s)*(water%gamma+1)/2.0
        k3 = -(1+water%gamma)/4.0/(rho0s**2*c0s**3)

        Rs = Rs + dRs*dt
        if(Rs-dt*dRs<dist.and.Rs>=dist)then
            t_arrive = t + (dist - Rs)/dRs
            ! call output_far_field(mesh, t_arrive)
        endif
        
    end subroutine update_far_field
    

    subroutine update_dt(mesh)
        use model
        use DG, only: npol
        implicit none
        type(DG_mesh), intent(in) :: mesh  
        ! local
        real dt0
        integer :: i
        real :: v
        v = us + cs - dRs
        if(v<0.0)then
            print*, 'Warning: negative wave speed in far field!'
            stop
        else
            dt = minval(mesh%cell_length)/v/(npol+1)*0.5
        end if
        if(Rs<dist.and.Rs+dt*dRs>=dist)then
            dt = min(dt, L/dRs*0.01)
        end if
        if(Rs>dist.and.Rs-L<dist)then
            dt = min(dt, 0.01*L/dRs)
        end if
        if(t-t_arrive+dt>tend)then
            dt = tend - t + t_arrive
        endif
        dt = 0.2*dt
        return
    end subroutine update_dt

    subroutine flux(U, F, zeta)
        use model
        implicit none
        real, intent(in):: U(:), zeta
        real, intent(out):: F(:)
        ! local
        real :: rho, p, vel, c, v
        F(1) = k1*U(1) + k2*U(1)**2/2.0 + k3*U(1)**3/3.0
        F(2) = 0.0
        return
    end subroutine flux

    subroutine source(U, S, r)
        use model
        implicit none
        real, intent(in):: U(:), r
        real, intent(out):: S(:)
        ! local
        real :: rho, p, c, vel, v
        p = U(1)
        vel = U(2)
        rho = water%dens_from_entropy(ss, p + p0)
        c = water%sound_spd(rho, p + p0)
        v = vel + c - dRs
        S(1) = -2*rho*c**2*vel/r*(1+ c/(v-2*c))
        S(2) = 0
        ! S=0
        return
    end subroutine source

    subroutine num_flux(UL, UR, F, zeta)
        use model, only: L, dL, dRb, water
        implicit none
        real, intent(in):: UL(:), UR(:), zeta
        real, intent(out):: F(:)
        ! local
        real FL(2), FR(2)
        real alpha
        real vl, vr
        call flux(UL, FL, zeta)
        call flux(UR, FR, zeta)
        ! calculate the maximum wave speed
        vl = k1 + k2*UL(1) + k3*UL(1)**2
        vr = k1 + k2*UR(1) + k3*UR(1)**2
        alpha = max(abs(vl), abs(vr))
        F = (FL + FR) / 2.0 + &
                     alpha * (UL - UR) / 2.0
        F(2) = 0.0
        return
    end subroutine num_flux

    subroutine bc_left(UR, F, zeta)
        use model, only: L, dL, dRb, pg
        implicit none
        real, intent(in):: UR(:), zeta
        real, intent(out):: F(:)
        ! local

        call flux(UR, F, zeta)

    end subroutine bc_left

    subroutine bc_right(UL, F, zeta)
        use model, only: rho0, p0, dRs, L, E0
        implicit none
        real, intent(in):: UL(:), zeta
        real, intent(out):: F(:)
        ! local
        call flux(UL, F, zeta)

    end subroutine bc_right


    subroutine output_far_field(mesh, tcurrent)
        use model
        use DG
        use polynomial
        implicit none
        type(DG_mesh), intent(in) :: mesh
        real, optional, intent(in) :: tcurrent
        ! local
        integer :: i
        real :: uedge(1, 2)
        real :: xs, basis(npol), usensor(3)
        real :: rho, vel, p
        integer :: cell_id
        real :: tprint
        
        if(present(tcurrent)) then
            tprint = tcurrent
        else
            tprint = t
        end if

        ! output far field pressure
        do i = 1, mesh%ncell
            uedge(:,1) = matmul(transpose(mesh%unk(:,:,i)), basis_face(:,1))
            uedge(:,2) = matmul(transpose(mesh%unk(:,:,i)), basis_face(:,2))
            write(101,'(2E20.10)', advance='no') uedge(1,1), uedge(1,2)
        end do
        write(101,*)

        ! output sensor data
        if(Rs>dist.and.Rs-L<dist)then ! sensor is in the fluid region
            do i=1, mesh%ncell
                if((mesh%node(i)*L+Rs-L)<=dist .and. &
                    (mesh%node(i+1)*L+Rs-L)>=dist) then
                    cell_id = i
                    exit
                end if
            end do
            if(cell_id>mesh%ncell.or. cell_id<1)then
                return
            endif
            xs = (dist - mesh%node(cell_id)*L - Rs + L)/mesh%cell_length(cell_id)/L
            do i=1, npol
                basis(i) = leg_ply(i-1)%evaluate(2.0*xs-1.0)
            end do
            p = dot_product((mesh%unk(:,1,cell_id)), basis)
            write(200, '(2E15.6)') tprint- t_arrive, p
            t_shock_end = tprint - t_arrive
            continue
        endif
        
        write(202, '(3E15.6)') tprint, Rs, pm
        return
    end subroutine output_far_field

    subroutine calculate_shock2bubble()
        use model, only: t, water, p0, L, Rs, c0, dRs &
            , shock2bubble, n_s2b
        use DG, only: ngauss, gauss_quad, basis_gauss
        ! use polynomial, only: 
        implicit none
        real :: Tk(ngauss), vel(ngauss), p(ngauss), r(ngauss)
        integer :: i, j
        real :: rho
        n_s2b = n_s2b + 1
        shock2bubble(1, n_s2b) = t + Rs/c0  ! offset time
        shock2bubble(2, n_s2b) = 0.0
        do i=1, far_DG%ncell
            r = 0.5*far_DG%cell_length(i)*(gauss_quad%pos + 1.0) + &
                far_DG%node(i)
            r = r*L + Rs - L
            p = matmul(transpose(basis_gauss), far_DG%unk(:,1,i))
            vel = matmul(transpose(basis_gauss), far_DG%unk(:,2,i))
            
            do j=1, ngauss
                rho = water%dens_from_entropy(ss, p(j) + p0)
                shock2bubble(2, n_s2b) = shock2bubble(2, n_s2b) + &
                    rho*vel(j)**2/r(j)*&
                    gauss_quad%weight(j)*far_DG%cell_length(i)/2.0*L
            end do
        enddo
        shock2bubble(2, n_s2b) = shock2bubble(2, n_s2b)*c0*c0/(c0+dRs)
        return
    end subroutine
    end module far_field