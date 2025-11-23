    module near_field
        use DG, only: DG_mesh
        implicit none
        type(DG_mesh) :: near_DG_1e, near_DG
    contains

    subroutine update_near_field(mesh)
        use model
        use DG
        implicit none
        type(DG_mesh), intent(inout) :: mesh
        ! local
        real ul(3), ur(3)
        real vel, rho
        ul = matmul(transpose(mesh%unk(:,:,1)), basis_face(:,1))
        ur = matmul(transpose(mesh%unk(:,:,mesh%ncell)), basis_face(:,2))
        dRb = ul(2)/ul(1)
        
        Pg = explosive%press_adb(&
            Pg0, explosive%dens0, explosive%dens0*Rb0**3/Rb**3)
        rho = ur(1)
        vel = ur(2)/rho
        Ps = water%press(rho, ur(3)-0.5*rho*vel**2)
        ! pm = ps-p0
        call update_shock_front()
        Rb = Rb + dRb*dt
        if(Rs<dist.and.Rs+dRs*dt>=dist)then
            print*, 'Sensor is reached by the shock front at time ', t + (dist - Rs)/dRs
            t_arrive = t + (dist - Rs)/dRs - t_offset
            write(200, '(2E15.6)') 0.0, 0.0
            write(200, '(2E15.6)') t_offset, ps - p0
        endif
        Rs = Rs + dRs*dt
        L = Rs - Rb
        dL = dRs - dRb

    end subroutine update_near_field
    

    subroutine update_dt(mesh)
        use model
        use DG, only: npol
        implicit none
        type(DG_mesh), intent(in) :: mesh  
        ! local
        real dt0
        integer :: i
        real :: c, vel, uhat, rho, E, chat1, chat2, umax, p
        dt0 = 1e10
        do i = 1, mesh%ncell
            rho = mesh%unk(1,1,i)
            E = mesh%unk(1,3,i)/rho
            vel = mesh%unk(1,2,i)/rho
            uhat = (vel - dL*mesh%cell_center(i) - dRb)/L
            p = water%press(rho, (E - 0.5*vel**2)*rho)
            c = water%sound_spd(rho, p)
            chat1 = (c - dL*mesh%cell_center(i) - dRb)/L
            chat2 = (c + dL*mesh%cell_center(i) + dRb)/L
            umax = max(abs(uhat+chat1), abs(uhat+chat2))
            dt0 = min(dt0, mesh%cell_length(i)/umax)
        end do
        dt = dt0*CFL/(npol+1)
        if(t+dt>tend) dt = tend - t+1e-8
        return
    end subroutine update_dt

    subroutine init_condi_near_field(ncell, nvar)
        use model, only: rho0, water, us, Ps, rhos, pi
        use DG, only: basis_face, basis_gauss, dbasis_gauss, ngauss, npol, gauss_quad
        ! use fsintegration
        implicit none
        integer, intent(in):: ncell, nvar
        integer :: i,j
        real :: rhs(npol), xg, ug

        call near_DG_1e%init_DG(nvar, 1)
        call near_DG%init_DG(nvar, ncell)

        ! set initial condition
        near_DG_1e%unk = 0.0
        near_DG_1e%unk(1,1,:) = rhos
        near_DG_1e%unk(1,2,:) = rhos*us
        near_DG_1e%unk(1,3,:) = water%internal_energy(Ps, rhos) + &
            0.5*rhos*us**2

        ! open a file for output
        open(101, file='output/near_field_grid.dat', status='replace', action='write', &
            form='formatted')
        do i=1, near_DG%ncell
            write(101,'(2E20.10)', advance = 'no') near_DG%node(i),  near_DG%node(i+1)
        end do
        close(101)
        open(101, file='output/near_field_rho.dat', status='replace', action='write', &
            form='formatted')
        open(102, file='output/near_field_v.dat', status='replace', action='write', &
            form='formatted')
        open(103, file='output/near_field_p.dat', status='replace', action='write', &
            form='formatted')
        return
    end subroutine init_condi_near_field

    subroutine flux(U, F, zeta)
        use model
        implicit none
        real, intent(in):: U(:), zeta
        real, intent(out):: F(:)
        ! local variables
        real :: rho, vel, p, uhat

        rho = U(1)
        vel = U(2)/rho
        p = water%press(rho, U(3)-0.5*rho*vel**2)
        uhat = (vel - zeta*dL - dRb)/L
        F(1) = rho * uhat
        F(2) = rho * vel*uhat + p/L
        F(3) = U(3)*uhat + vel*p/L
        return
    end subroutine flux

    subroutine source(U, S, r)
        use model
        implicit none
        real, intent(in):: U(:), r
        real, intent(out):: S(:)
        ! local
        real :: rho, vel, p
        rho = U(1)
        vel = U(2)/rho
        p = water%press(rho, U(3)-0.5*rho*vel**2)
        S(1) = -2.0*U(1)*vel/r
        S(2) = -2.0*U(2)*vel/r
        S(3) = -2.0*(U(3)+p)*vel/r
        S = S-U*dL/L
        return
    end subroutine source

    subroutine num_flux(UL, UR, F, zeta)
        use model, only: L, dL, dRb, water
        implicit none
        real, intent(in):: UL(:), UR(:), zeta
        real, intent(out):: F(:)
        ! local
        real FL(3), FR(3)
        real alpha, CL, CR
        real :: rhoL, velL, pL, uhatL
        real :: rhoR, velR, pR, uhatR
        ! left state

        rhoL = UL(1)
        velL = UL(2)/rhoL
        pL = water%press(rhoL, UL(3)-0.5*rhoL*velL**2)
        CL = water%sound_spd(rhoL, pL)
        uhatL = (velL - zeta*dL - dRb)/L
        FL(1) = rhoL * uhatL
        FL(2) = rhoL * velL*uhatL + pL/L
        FL(3) = UL(3)*uhatL + velL*pL/L
        ! right state
        rhoR = UR(1)
        velR = UR(2)/rhoR
        pR = water%press(rhoR, UR(3)-0.5*rhoR*velR**2)
        CR = water%sound_spd(rhoR, pR)
        uhatR = (velR - zeta*dL - dRb)/L
        FR(1) = rhoR * uhatR
        FR(2) = rhoR * velR*uhatR + pR/L
        FR(3) = UR(3)*uhatR + velR*pR/L
        ! LF flux
        alpha = max(abs(uhatL)+CL/L, abs(uhatR)+CR/L)
        F = 0.5 * (FL + FR) + 0.5 * alpha * (UL - UR)
        return
    end subroutine num_flux

    subroutine bc_left(UR, F, zeta)
        use model, only: L, dL, dRb, pg
        implicit none
        real, intent(in):: UR(:), zeta
        real, intent(out):: F(:)
        ! local
        real rhoR, velR, uhatR
        rhoR = UR(1)
        velR = UR(2)/rhoR
        uhatR = (velR - zeta*dL - dRb)/L
        F(1) = 0
        F(2) = pg/L
        F(3) = dRb*pg/L

    end subroutine bc_left

    subroutine bc_right(UL, F, zeta)
        use model, only: rho0, p0, dRs, L, E0
        implicit none
        real, intent(in):: UL(:), zeta
        real, intent(out):: F(:)
        ! local
        real rhoL, velL, pL
        rhoL = UL(1)
        velL = UL(2)/rhoL
        F(1) = -rho0*dRs/L
        F(2) = p0/L
        F(3) = -dRs*E0/L
    end subroutine bc_right


    subroutine output_near_field(mesh, tcurrent)
        use model
        use DG
        use polynomial
        implicit none
        type(DG_mesh), intent(in) :: mesh
        real, optional, intent(in) :: tcurrent
        ! local
        integer :: i
        real :: uedge(3, 2), pedge(2)
        real :: xs, basis(npol), usensor(3)
        real :: rho, vel, p
        integer :: cell_id
        real :: tprint

        if(present(tcurrent)) then
            tprint = tcurrent
        else
            tprint = t
        end if

        ! output field data
        do i = 1, mesh%ncell
            uedge(:,1) = matmul(transpose(mesh%unk(:,:,i)), basis_face(:,1))
            uedge(:,2) = matmul(transpose(mesh%unk(:,:,i)), basis_face(:,2))
            pedge(1) = water%press(uedge(1,1), uedge(3,1)-0.5*uedge(2,1)**2/uedge(1,1))
            pedge(2) = water%press(uedge(1,2), uedge(3,2)-0.5*uedge(2,2)**2/uedge(1,2))
            write(101,'(2E20.10)', advance='no') uedge(1,1), uedge(1,2)
            write(102,'(2E20.10)', advance='no') uedge(2,1)/uedge(1,1), uedge(2,2)/uedge(1,2)
            write(103,'(2E20.10)', advance='no') pedge(1)-p0, pedge(2)-p0
        end do
        write(101,*)
        write(102,*)
        write(103,*)

        ! output sensor data
        if(Rs>dist.and.Rb<dist)then ! sensor is in the fluid region
            do i=1, mesh%ncell
                if((mesh%node(i)*L+Rb)<=dist .and. &
                    (mesh%node(i+1)*L+Rb)>=dist) then
                    cell_id = i
                    exit
                end if
            end do
            xs = (dist - mesh%node(cell_id)*L - Rb)/mesh%cell_length(cell_id)/L
            do i=1, npol
                basis(i) = leg_ply(i-1)%evaluate(2.0*xs-1.0)
            end do
            usensor = matmul(transpose(mesh%unk(:,:,cell_id)), basis)
            rho = usensor(1)
            vel = usensor(2)/rho
            p = water%press(rho, usensor(3)-0.5*rho*vel**2)
            t_shock_end = tprint - t_arrive
            write(200, '(2E15.6)') t_shock_end, p - p0
            continue
        elseif(Rb>=dist) then   ! sensor is in the bubble region
            t_shock_end = tprint - t_arrive
            write(200, '(2E15.6)') t_shock_end, pg - p0
            continue
        endif
        
        ! bubble data
        write(201, '(3E15.6)') tprint, Rb, 0.0
        ! shock wave data
        write(202, '(3E15.6)') tprint, Rs, ps - p0
        return
    end subroutine output_near_field

    subroutine near_shock2bubble()
        use model, only: shock2bubble, n_s2b, Rs, Ps, t, cs, Rb, L
        use DG, only: ngauss, gauss_quad, basis_gauss
        implicit none
        ! local
        integer :: i
        real :: xg(ngauss),  ug(3, ngauss), rg(ngauss), rho(ngauss), vel(ngauss)
        ! code

        n_s2b = n_s2b + 1
        shock2bubble(1,n_s2b) = t + Rs/cs
        shock2bubble(2,n_s2b) = 0.0
        do i=1, near_DG%ncell / 2
            xg = 0.5*near_DG%cell_length(i)*(gauss_quad%pos + 1.0) + &
                near_DG%node(i)
            rg = xg*L+Rb
            ug = matmul(transpose(near_DG%unk(:,:,i)), basis_gauss)
            rho = ug(1,:)
            vel = ug(2,:)/rho
            shock2bubble(2,n_s2b) = shock2bubble(2,n_s2b) + &
                sum(rho*vel**2 * &
                gauss_quad%weight(1:ngauss)) * &
                near_DG%cell_length(i)/2.0*L
        enddo
    end subroutine near_shock2bubble

    end module near_field