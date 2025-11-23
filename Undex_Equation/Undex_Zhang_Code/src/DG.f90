    module DG
        use fsintegration
        implicit none
        ! solution variables
        type DG_mesh
        integer :: nvar, ncell, nface
        real, allocatable :: unk(:,:,:), &  ! Unknown variables, npol*nvar*ncell
                          flux(:,:), &    ! numerical Flux variables, nvar*nface
                          unk_bak(:,:,:),&  ! Backup of unknown variables, npol*nvar*ncell
                          imass(:,:,:),&   ! Mass matrix, npol*npol*ncell
                          val_face(:,:,:),&! face values, nvar*2*ncell+1
                          node(:),&  ! Node coordinates, ncell+1
                          cell_length(:),&  ! Cell lengths, ncell
                            cell_center(:)  ! Cell centers, ncell
            contains
            procedure :: init_DG
            procedure :: march_DG

        end type DG_mesh

        ! basis
        integer, parameter :: npol=3, ngauss=3
        real,ALLOCATABLE:: basis_gauss(:,:), &  ! npols*ngauss
                            basis_face(:,:),&     ! npols*2
                            dbasis_gauss(:,:),&   ! npols*ngauss
                            basis_center(:),&  ! npols
                            ibasis_gauss(:,:)  ! npols*ngauss
        type(gauss_legendre):: gauss_quad


        interface
            subroutine core_func(ug, rg, nvar, ans)
                implicit none
                integer, intent(in):: nvar
                real, intent(in):: ug(:, :), rg(:)
                real :: ans(nvar)
            end subroutine
            module function integrate(this, core) result(res)
                implicit none
                class(DG_mesh), intent(in):: this
                procedure(core_func):: core
                real res(this%nvar)
            end function
        end interface

        interface 
            subroutine flux_func(U, flux, r)
                real, intent(in):: U(:), r
                real, intent(out):: flux(:)
            end subroutine flux_func
            subroutine source_func(U,source, r)
                real, intent(in):: U(:), r
                real, intent(out):: source(:)
            end subroutine source_func
            subroutine num_flux_func(UL, UR,num_flux, r) 
                real, intent(in):: UL(:), UR(:), r
                real, intent(out):: num_flux(:)
            end subroutine num_flux_func
            subroutine bc_L_func(UR, F, r)
                real, intent(in):: UR(:), r
                real, intent(out):: F(:)
            end subroutine bc_L_func
            subroutine bc_R_func(UL, F, r)
                real, intent(in):: UL(:), r
                real, intent(out):: F(:)
            end subroutine bc_R_func
        end interface

        interface
        module subroutine init_DG(this, nvar_inp, ncell_inp)
            class(DG_mesh), intent(inout) :: this
            integer, intent(in) :: nvar_inp, ncell_inp
        end subroutine init_DG
        module subroutine march_DG(this, dt, flux, num_flux, source, &
            bc_left, bc_right)
        class(DG_mesh), intent(inout) :: this
        procedure(flux_func) :: flux
        procedure(num_flux_func) :: num_flux
        procedure(source_func) :: source
        procedure(bc_L_func) :: bc_left
        procedure(bc_R_func) :: bc_right
        real, intent(in) :: dt
        end subroutine march_DG
        end interface


        contains

        subroutine project_1e_to_multi(mesh_1e, mesh_multi)
        use model
        use polynomial
        implicit none
        type(DG_mesh), intent(in) :: mesh_1e
        type(DG_mesh), intent(inout) :: mesh_multi
        integer :: i, j, k
        real :: r, val(3, ngauss), dx, zetag(ngauss), basis(npol)
        real :: rhs(npol, 3)
        ! project 1e mesh to multi mesh
        do i = 1, mesh_multi%ncell
            zetag = 0.5*mesh_multi%cell_length(i)*(gauss_quad%pos + 1.0) + &
                mesh_multi%node(i)
            zetag = zetag * 2 - 1.0
            rhs = 0
            do k = 1, ngauss
                do j=1, npol
                    basis(j) = leg_ply(j-1)%evaluate(zetag(k))
                end do
                val(:,k) = matmul(transpose(mesh_1e%unk(:,:,1)),basis)
                rhs(:,1) = rhs(:,1) + val(1,k)*gauss_quad%weight(k) * basis_gauss(:,k)
                rhs(:,2) = rhs(:,2) + val(2,k)*gauss_quad%weight(k) * basis_gauss(:,k)
                rhs(:,3) = rhs(:,3) + val(3,k)*gauss_quad%weight(k) * basis_gauss(:,k)
            end do
            rhs = rhs * 0.5 * mesh_multi%cell_length(i)
            mesh_multi%unk(:,:,i) = matmul(mesh_multi%imass(:,:,i), rhs)
        enddo
        return
        end subroutine project_1e_to_multi


        subroutine init_basis()
            use polynomial
            use fsintegration
            implicit none
            integer :: i,j,k
            
            call initiate_legendre(4)

            ! transfer basis integration 
            ! from zeta to 1 instead of -1 to zeta
            do i=1, npol
                leg_ply(i-1)%ia(0) = leg_ply(i-1)%a(0)*2.0
                do j=1, leg_ply(i-1)%order
                    leg_ply(i-1)%ia(j) = -leg_ply(i-1)%ia(j)
                end do
            end do

            allocate(basis_gauss(npol,ngauss), &
                     basis_face(npol,2), &
                     dbasis_gauss(npol,ngauss),&
                     basis_center(npol),&
                     ibasis_gauss(npol,ngauss))

            ! Initialize basis functions
            call gauss_quad%gl_init_legendre(ngauss)
            do i = 1, npol
                do j = 1, ngauss
                    basis_gauss(i,j) = leg_ply(i-1)%evaluate(gauss_quad%pos(j))
                    dbasis_gauss(i,j) = leg_ply(i-1)%derivative(gauss_quad%pos(j))
                    ibasis_gauss(i,j) = leg_ply(i-1)%integrate(gauss_quad%pos(j))
                end do
                basis_center(i) = leg_ply(i-1)%evaluate(0.0)
            end do

            ! Initialize face basis functions
            do i = 1, npol
                basis_face(i,1) = leg_ply(i-1)%evaluate(-1.0)   ! left face
                basis_face(i,2) = leg_ply(i-1)%evaluate(1.0)    ! right face
            end do

        end subroutine init_basis


    end module