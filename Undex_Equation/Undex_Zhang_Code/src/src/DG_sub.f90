    submodule (DG) DG_sub
        implicit none

    contains
    module subroutine init_DG(this, nvar_inp, ncell_inp)
        use polynomial
        use fsintegration
        use linear_equation
        class(DG_mesh), intent(inout) :: this
        integer, intent(in) :: nvar_inp, ncell_inp
        integer :: i,j, k, ierr
        real :: val(ngauss)
        real :: mass0(npol,npol), dx

        ! set parameters
        this%nvar = nvar_inp
        this%ncell = ncell_inp

        ! allocate arrays
        allocate(this%unk(npol,this%nvar,this%ncell))
        allocate(this%flux(this%nvar,this%ncell+1))
        allocate(this%unk_bak(npol,this%nvar,this%ncell))
        allocate(this%imass(npol,npol,this%ncell),&
            this%val_face(2,this%nvar,this%ncell+1),&
            this%node(this%ncell+1),&
            this%cell_length(this%ncell),&
            this%cell_center(this%ncell))
        

        ! Initialize mesh (uniform mesh in [0,1])
        do i = 1, this%ncell+1
            this%node(i) = real(i-1)/real(this%ncell)
        end do
        do i = 1, this%ncell
            this%cell_length(i) = this%node(i+1) - this%node(i)
            this%cell_center(i) = 0.5*(this%node(i+1)+this%node(i))
        end do
        ! Initialize mass matrix
        do i = 1, this%ncell
            dx = this%cell_length(i)
            mass0 = 0
            do j = 1, npol
                do k=1,npol
                    ! Initialize mass matrix
                    val(1:ngauss) = basis_gauss(j,1:ngauss) * &
                                    basis_gauss(k,1:ngauss)
                    mass0(j,k) = gauss_quad%gl_integrate(val(1:ngauss), dx)
                end do
            end do
            ! calculate invert
            call inv(this%imass(1:npol,1:npol,i), mass0(1:npol,1:npol), ierr)
        end do
    end subroutine init_DG



    module subroutine march_DG(this, dt, flux, num_flux, source, &
            bc_left, bc_right)
        use model, only: L, dL, Rb, Rs
        class(DG_mesh), intent(inout) :: this
        procedure(flux_func) :: flux
        procedure(num_flux_func) :: num_flux
        procedure(source_func) :: source
        procedure(bc_L_func) :: bc_left
        procedure(bc_R_func) :: bc_right
        real, intent(in) :: dt
        ! local
        integer :: i,j,k,m
        real :: rhs(npol), ug(this%nvar, ngauss), &
            sg(this%nvar, ngauss), flxg(this%nvar, ngauss)
        real :: rg(ngauss),zetag(ngauss), dunk(npol)
        ! calculate face values
        do i=1, this%ncell
            do m=1, this%nvar
                this%val_face(2,m,i) = 0.0
                this%val_face(1,m,i+1) = 0.0
                do j=1, npol
                    this%val_face(2,m,i) = this%val_face(2,m,i) + &
                        basis_face(j,1)*this%unk(j,m,i)
                    this%val_face(1,m,i+1) = this%val_face(1,m,i+1) + &
                        basis_face(j,2)*this%unk(j,m,i)
                end do
            end do
        end do
        ! apply boundary conditions
        call bc_left(this%val_face(2,:,1), this%flux(:,1), 0.0)
        call bc_right(this%val_face(1,:,this%ncell+1), this%flux(:,this%ncell+1), 1.0)

        ! calculate numerical flux at faces
        do i=2, this%ncell
            call num_flux(this%val_face(1,1:this%nvar,i), &
                          this%val_face(2,1:this%nvar,i), &
                          this%flux(1:this%nvar,i), this%node(i))
        end do

        ! calculate derivative and march
        do i=1, this%ncell
            ! internal flux
            ug = matmul(transpose(this%unk(:,:,i)), basis_gauss)
            zetag = 0.5*this%cell_length(i)*(gauss_quad%pos + 1.0) + &
                this%node(i)
            rg = zetag*L+Rs - L
            
            do j=1, ngauss
                ! flux term
                call flux(ug(:,j), flxg(:,j), zetag(j))
                ! source term
                call source(ug(:,j), sg(:,j), rg(j))
            enddo
            do j=1, this%nvar
                ! boundary flux
                rhs = this%flux(j,i)*basis_face(:,1) - &
                    this%flux(j,i+1)*basis_face(:,2)
                ! intelnal flux
                rhs(:) = rhs(:) + &
                    matmul(dbasis_gauss, flxg(j,:)*&
                                gauss_quad%weight(1:ngauss))*&
                                this%cell_length(i)/2.0/&
                                this%cell_length(i)*2  ! scale dbasis to cell length
                ! source term
                rhs = rhs + &
                    matmul(basis_gauss, sg(j,:)*&
                                gauss_quad%weight(1:ngauss))*this%cell_length(i)/2.0
                ! update unknowns
                dunk = matmul(this%imass(:,:,i), rhs)
                this%unk(:,j,i) = this%unk(:,j,i) + &
                    dt*dunk
            enddo
        enddo

    end subroutine march_DG

    module function integrate(this, core) result(res)
        use model, only: L, Rs
        implicit none
        class(DG_mesh), intent(in):: this
        procedure(core_func):: core
        ! output
        real res(this%nvar)
        ! local
        integer i
        real :: ug(this%nvar, ngauss), rg(ngauss),zetag(ngauss)
        real :: integrant(this%nvar,ngauss)
        ! code
        res = 0
        do i=1, this%ncell
            zetag = 0.5*this%cell_length(i)*(gauss_quad%pos + 1.0) + &
                this%node(i)
            rg = zetag*L+Rs - L
            ug = matmul(transpose(this%unk(:,:,i)),basis_gauss)
            ! call core_func(ug, rg, this%nvar, integrant)
        enddo
        res = res*L/2.0
        return
    end function

    end submodule DG_sub