 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
 !!!!! CopyRight by Liu Yunlong   Email:yunlong_liu@hrbeu.edu.cn !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !

module linear_equation

   interface
      subroutine lie(A,X,B)
         integer i,j,k,N,m
         real A(:,:),X(:),B(:)
      end subroutine
   end interface
   interface gmres
      module procedure gmres_multi
      module procedure gmres_single
   end interface
contains
   subroutine gmres_multi(A, B, X, restart, tol, maxiter, ierr)
      ! GMRES solver for AX = B with multiple right-hand sides
      ! Inputs:
      !   A        : real matrix (n,n)
      !   B        : right-hand side matrix (n, nrhs)
      !   restart  : restart parameter (integer)
      !   tol      : tolerance for convergence (real)
      !   maxiter  : maximum number of iterations (integer)
      ! In/out:
      !   X        : initial guess (input), solution (output) (n, nrhs)
      ! Output:
      !   ierr     : 0 if converged, 1 otherwise (optional, size nrhs)
      implicit none
      real, intent(in) :: A(:,:), B(:,:)
      real, intent(inout) :: X(:,:)
      integer, intent(in) :: restart, maxiter
      real, intent(in) :: tol
      integer, optional, intent(out) :: ierr(:)

      integer :: n, nrhs, i
      integer, allocatable :: ierr_local(:)
      n = size(A,1)
      nrhs = size(B,2)
      if (present(ierr)) then
         if (size(ierr) /= nrhs) stop 'ierr size mismatch'
      else
         allocate(ierr_local(nrhs))
      end if

      do i = 1, nrhs
         call gmres_single(A, B(:,i), X(:,i), restart, tol, maxiter, &
            ierr_local(i))
      end do

      if (present(ierr)) ierr = ierr_local

      deallocate(ierr_local)
   end subroutine gmres_multi
   subroutine gmres_single(A, b, x, restart, tol, maxiter, ierr)
      ! GMRES solver for Ax = b
      ! Inputs:
      !   A        : real matrix (n,n)
      !   b        : right-hand side vector (n)
      !   restart  : restart parameter (integer)
      !   tol      : tolerance for convergence (real)
      !   maxiter  : maximum number of iterations (integer)
      ! In/out:
      !   x        : initial guess (input), solution (output)
      ! Output:
      !   ierr     : 0 if converged, 1 otherwise (optional)
      implicit none
      real, intent(in) :: A(:,:), b(:)
      real, intent(inout) :: x(:)
      integer, intent(in) :: restart, maxiter
      real, intent(in) :: tol
      integer, optional, intent(out) :: ierr

      integer :: n, i, j, k, m, iter, itot
      real :: beta, resid, temp
      real, allocatable :: V(:,:), H(:,:), cs(:), sn(:), e1(:), s(:), y(:)
      real, allocatable :: w(:), r(:)
      n = size(b)
      m = restart
      allocate(V(n,m+1), H(m+1,m), cs(m), sn(m), e1(m+1), s(m+1), y(m), w(n), r(n))

      x = 0.0
      r = b - matmul(A, x)
      beta = norm2(r)
      resid = beta
      if (present(ierr)) ierr = 1
      if (resid < tol) then
         if (present(ierr)) ierr = 0
         deallocate(V,H,cs,sn,e1,s,y,w,r)
         return
      end if

      itot = 0
      do while (itot < maxiter)
         V(:,1) = r / beta
         e1 = 0.0
         e1(1) = beta
         H = 0.0
         s = e1
         do k = 1, m
            w = matmul(A, V(:,k))
            do j = 1, k
               H(j,k) = dot_product(w, V(:,j))
               w = w - H(j,k) * V(:,j)
            end do
            H(k+1,k) = norm2(w)
            if (H(k+1,k) /= 0.0) then
               V(:,k+1) = w / H(k+1,k)
            else
               V(:,k+1) = 0.0
            end if

            ! Apply Givens rotations
            do i = 1, k-1
               temp = cs(i)*H(i,k) + sn(i)*H(i+1,k)
               H(i+1,k) = -sn(i)*H(i,k) + cs(i)*H(i+1,k)
               H(i,k) = temp
            end do

            ! Compute i-th rotation
            call givens(H(k,k), H(k+1,k), cs(k), sn(k))
            temp = cs(k)*s(k)
            s(k+1) = -sn(k)*s(k)
            s(k) = temp

            temp = cs(k)*H(k,k) + sn(k)*H(k+1,k)
            H(k+1,k) = 0.0
            H(k,k) = temp

            resid = abs(s(k+1))
            if (resid < tol) exit
         end do

         ! Solve upper triangular system
         y(1:k) = s(1:k)
         do i = k, 1, -1
            y(i) = y(i) / H(i,i)
            if (i > 1) y(1:i-1) = y(1:i-1) - H(1:i-1,i) * y(i)
         end do

         x = x + matmul(V(:,1:k), y(1:k))
         r = b - matmul(A, x)
         resid = norm2(r)
         itot = itot + k
         if (resid < tol) then
            if (present(ierr)) ierr = 0
            exit
         end if
      end do

      deallocate(V,H,cs,sn,e1,s,y,w,r)
   contains
      pure real function norm2(v)
         real, intent(in) :: v(:)
         norm2 = sqrt(sum(v*v))
      end function norm2

      subroutine givens(a, b, c, s)
         real, intent(in) :: a, b
         real, intent(out) :: c, s
         real :: r
         if (b == 0.0) then
            c = 1.0
            s = 0.0
         else if (abs(b) > abs(a)) then
            r = a / b
            s = 1.0 / sqrt(1.0 + r*r)
            c = s*r
         else
            r = b / a
            c = 1.0 / sqrt(1.0 + r*r)
            s = c*r
         end if
      end subroutine givens
   end subroutine gmres_single

   function inv33(A)
      real inv33(3,3)
      real A(3,3),det
      real a1122,a1123,a1221,a1223,a1322,a1321

      a1122 = a(1,1)*a(2,2)
      a1123 = a(1,1)*a(2,3)
      a1221 = a(1,2)*a(2,1)
      a1223 = a(1,2)*a(2,3)
      a1322 = a(1,3)*a(2,2)
      a1321 = a(1,3)*a(2,1)
      det = a1122*A(3,3)-a1123*A(3,2)-&
         a1221*A(3,3)+a1223*A(3,1)-&
         a1322*A(3,1)+a1321*A(3,2)

      inv33=0
      inv33(1,1)=A(2,2)*A(3,3)-A(2,3)*A(3,2)
      inv33(1,2)=A(1,3)*A(3,2)-A(1,2)*A(3,3)
      inv33(1,3)=a1223-a1322
      inv33(2,1)=A(2,3)*A(3,1)-A(2,1)*A(3,3)
      inv33(2,2)=A(1,1)*A(3,3)-A(1,3)*A(3,1)
      inv33(2,3)=A(1,3)*A(2,1)-a1123
      inv33(3,1)=A(2,1)*A(3,2)-A(2,2)*A(3,1)
      inv33(3,2)=A(1,2)*A(3,1)-A(1,1)*A(3,2)
      inv33(3,3)=a1122-a1221

      inv33=inv33/det
      return
   end function

   function inv22(A)
      real inv22(2,2)
      real A(2,2),det
      det = A(1,1)*A(2,2)-A(1,2)*A(2,1)
      inv22(1,:) = [A(2,2),-A(1,2)]
      inv22(2,:) = [-A(2,1),A(1,1)]
      inv22 = inv22/det
      return
   end function

   subroutine inv_lie(A0,X,tp)
      real,intent(in):: A0(:,:)
      real,intent(inout)::x(:,:)
      integer,optional:: tp
      integer i,j,k,N,m
      real det,max
      real, allocatable::mid(:),A(:,:),B(:,:)
      logical linv
      N=size(A0,1)
      allocate (A(N,N))
      allocate(mid(N))
      allocate(B(N,N))
      A = A0
      linv = .true.
      if(present(tp))then
         if(tp.ne.0)then
            linv = .false.
         endif
      endif
      det=1
      B = 0
      if (.not.linv)then
         B=X
      else
         forall(i=1:N)
            B(i,i)=1
         end forall
      endif
      do k=1,N-1
         max=0
         do i=k,N
            if(abs(A(i,k))>max)then
               m=i
               max=abs(A(i,k))
            endif
         enddo
         if(m.ne.k)then
            mid(k:N)=A(m,K:N);A(m,K:N)=A(k,K:N);A(k,K:N)=mid(K:N)
            mid=B(m,:);B(m,:)=B(k,:);B(k,:)=mid
            det=-det
         endif
         do i=k+1,N
            A(i,k)=A(i,k)/A(k,k)
            do j=k+1,N
               A(i,j)=A(i,j)-A(i,k)*A(k,j)
            enddo
            B(i,:)=B(i,:)-A(i,k)*B(k,:)
            det=det*A(k,k)
         enddo
      enddo
      B(N,:)=B(N,:)/A(N,N)
      do k=1,N-1
         i=N-k
         do j=1,N
            B(i,j)=(B(i,j)-sum(A(i,i+1:N)*B(i+1:N,j)))/A(i,i)
         enddo
      enddo
      det=det*A(N,N)
      X=B
      deallocate(A,mid,B)
   end subroutine inv_lie

   subroutine inv(mat_out,mat_in,sta)
      implicit none
      integer i,j,k,N
      integer,optional::sta
      real mat_in(:,:),mat_out(:,:)
      real ierror
      integer l_inv,info
      real t_inv,temp3,temp
      integer,allocatable:: pivot(:)
      real,allocatable:: work(:)
      N=size(mat_in(:,1))
      allocate(pivot(N),work(N))
      ierror = 0
      if(present(sta))sta=1
      mat_out= mat_in
      info = 0
      do k = 1, N-1
         l_inv = k
         do i = k+1, N  
            if ( abs ( mat_out(i,k) ) > abs ( mat_out(l_inv,k) ) ) then
               l_inv = i
            end if
         end do
         pivot(k) = l_inv
         if ( abs(mat_out(l_inv,k)) < 1e-10 ) then 
            print*, mat_out(l_inv,k)
            info = k
            write ( *, *) '  Zero pivot on step ', info
            if(present(sta))sta=-1
            return
         end if
         if ( l_inv /= k ) then   
            temp3=mat_out(k,k)
            mat_out(k,k)=mat_out(l_inv,k)
            mat_out(l_inv,k)=mat_out(k,k)
         end if
         mat_out(k+1:N,k) =-mat_out(k+1:N,k) / mat_out(k,k)
         do j = k+1, N
            if ( l_inv /= k ) then 
               temp3=mat_out(k,j)
               mat_out(k,j)=mat_out(l_inv,j)
               mat_out(l_inv,j)=mat_out(k,j)
            end if
            mat_out(k+1:N,j)=mat_out(k+1:N,j)+&
               mat_out(k+1:N,k)*mat_out(k,j)
         end do
      end do
      pivot(N) = N
      if ( abs(mat_out(N,N)) < 1e-10 ) then
         info = N
         write ( *, * ) '  Zero pivot on step ', info
         if(present(sta))sta=-1
         return
      end if
      if ( ierror /= 0 ) then
         return
      end if
      do k = 1, N
         mat_out(k,k) = 1.0E+00 / mat_out(k,k)
         mat_out(1:k-1,k) = -mat_out(1:k-1,k) * mat_out(k,k)
         do j = k + 1, N
            temp = mat_out(k,j)
            mat_out(k,j) = 0.0E+00
            mat_out(1:k,j) = mat_out(1:k,j) + temp * mat_out(1:k,k)
         end do
      end do
      do k = N - 1, 1, -1
         work(k+1:N) = mat_out(k+1:N,k)
         mat_out(k+1:N,k) = 0.0E+00
         do j = k + 1, N
            mat_out(1:N,k) = mat_out(1:N,k) + mat_out(1:N,j) * work(j)
         end do
         if ( pivot(k) /= k ) then
            do i = 1, N
               temp3=mat_out(i,pivot(k))
               mat_out(i,pivot(k))=mat_out(i,k)
               mat_out(i,k)=mat_out(i,pivot(k))
            end do
         end if
      end do
   end subroutine inv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine Gauss_jordan(A,ANS,S)
      implicit none
      real A(:,:),S(:),ANS(:)
      real,allocatable:: B(:,:)
      integer::i,N
      N=size(A,1)
      allocate(B(N,N))
      B=A
      ANS=S
      call Upper(B,ANS,N)
      call Lower(B,ANS,N)
      forall(i=1:N)
         ANS(i)=ANS(i)/B(i,i)
      end forall
      return
   end subroutine Gauss_jordan
!!!-------------------------------
   Subroutine Upper(M,S,N)
      implicit none
      integer::N
      real M(N,N),S(N)
      integer:: i,j
      real E
      do i=1,N-1
         do j=i+1,N
            E=M(J,I)/M(i,i)
            M(J,i:N)=M(J,i:N)-M(i,i:N)*E
            S(J)=S(J)-S(i)*E
         enddo
      enddo
      return
   end subroutine upper
!!----------------------------------
   Subroutine Lower(M,S,N)
      implicit none
      integer::N
      real M(N,N),S(N)
      integer:: i,j
      real E
      do i=N,2,-1
         do j=i-1,1,-1
            E=M(J,I)/M(i,i)
            M(J,i:N)=M(J,i:N)-M(i,i:N)*E
            S(J)=S(J)-S(i)*E
         enddo
      enddo
      return
   end subroutine Lower
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine GS_FAST(A,X,rhs)
      ! solve linear equations with Gauss-Seidel Method
      integer num,i,j,k,mid,m,jiasu,N 
      integer,parameter::max_inter=50,Nop=5
      real:: A(:,:),X(:),rhs(:),err
      real:: mid_b,max
      real,allocatable::RES(:,:),X_bak(:,:)
      real,allocatable:: S(:,:),F(:),cons(:)&
         ,X0(:),mid_A(:),invv(:,:)
      integer, allocatable:: sor(:)
      integer bubble_count,sta,II,JJ,Nacc
      N=size(A,1)
      allocate (RES(max_inter,N))
      allocate (X_bak(max_inter,N))
      allocate (S(max_inter+1,max_inter+1))
      allocate (invv(max_inter+1,max_inter+1))
      allocate (F(max_inter+1))
      allocate (cons(max_inter+1))
      allocate (X0(N))
      allocate(mid_A(N))
      allocate(sor(N))
      sta=1
      jiasu=1
      Nacc=0
      X0=X
      k=0
      err=10
      X_bak=0
      RES=0
      do i=1,N
         sor(i)=i
      enddo
      X0=0
      do while(k<max_inter.and.err>1e-13)
         k=k+1
         err=0
         do i=1,N
            X(i)=rhs(i)
            do j=1,i-1
               X(i)=X(i)-A(i,j)*X(j)
            enddo
            do j=i+1,N
               X(i)=X(i)-A(i,j)*X0(j)
            enddo
            X(i)=X(i)/A(i,i)
         enddo
         RES(k,:)=rhs-matmul(A,X)
         X_bak(k,:)=X

         if(jiasu==1.and.k>Nop)then
            S=0
            do i=1,Nop
               II=i-Nop+k
               do j=1,Nop
                  JJ=j-Nop+k
                  do m=1,N
                     S(i,j)=S(i,j)+RES(II,m)*RES(JJ,m)
                  enddo
               enddo
            enddo
            S(Nop+1,1:Nop)=1;S(1:Nop,Nop+1)=0.5;
            cons(1:Nop)=0;cons(Nop+1)=1;

!		  if (k>20)then  !!!
!		    call inv(invv(1:k+1,1:k+1),S(1:k+1,1:k+1),sta)
!		    if (sta==1)then
!		        F(1:k+1)=matmul(invv(1:k+1,1:k+1),cons(1:k+1))
!		    endif
!
!		  else
            call lie(S(1:Nop+1,1:Nop+1),F(1:Nop+1),cons(1:Nop+1))
!		  endif
            if (sum(F(1:Nop+1))<100000)then
               RES(k,:)=matmul(F(1:Nop),RES(k-Nop+1:k,:))
               X(:)=matmul(F(1:Nop),X_bak(k-Nop+1:k,:))
               Nacc=Nacc+1
            endif
            X_bak(k,:)=X
         endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!c
         err=maxval(abs(RES(k,:)))
!        if (sum(abs(X))==0) then
!             err=0
!        else
!             err=err/sum(abs(X))
!        endif
         X0=X
      enddo
      deallocate (RES,S,F,cons,X_bak)
      if (N>100.and.k>20)write(*,*) 'iteratal=',k,&
         'GS method Max Erro=',err
   end subroutine GS_FAST
   subroutine chase(A,ans,B)
      real A(:,:),ans(:),B(:)
      integer i,j,k,N
      real,allocatable::alpha(:),beta(:),y(:)
      N=size(A(:,1))
      allocate(beta(N))
      allocate(y(N))
      beta(1)=A(1,3)/A(1,2)
      do i=2,N-1
         beta(i)=A(i,3)/(A(i,2)-A(i,1)*beta(i-1))
      enddo
      y(1)=B(1)/A(1,2)
      do i=2,N
         y(i)=(B(i)-A(i,1)*y(i-1))/(A(i,2)-A(i,1)*beta(i-1))
      enddo
      ans(N)=y(N)
      do i=1,N-1
         j=N-i
         ans(j)=y(j)-beta(i)*ans(j+1)
      enddo
   end subroutine chase
end module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine lie(A,X,B)
!---------------------
   integer i,j,k,N,m
   real A(:,:),X(:),B(:),det,max
   real, allocatable::mid(:),A0(:,:),B0(:)
   N=size(A,1)
   allocate (A0(N,N))
   allocate(mid(N))
   allocate(B0(N))
   det=1
   A0=A
   B0=B
   do k=1,N-1
      max=0
      do i=k,N
         if(abs(A(i,k))>max)then
            m=i
            max=abs(A(i,k))
         endif
      enddo
      if(m.ne.k)then
         mid(k:N)=A(m,K:N);A(m,K:N)=A(k,K:N);A(k,K:N)=mid(K:N)
         mid(1)=B(m);B(m)=B(k);B(k)=mid(1)
         det=-det
      endif
      do i=k+1,N
         A(i,k)=A(i,k)/A(k,k)
         do j=k+1,N
            A(i,j)=A(i,j)-A(i,k)*A(k,j)
         enddo
         B(i)=B(i)-A(i,k)*B(k)
         det=det*A(k,k)

      enddo
   enddo
   B(N)=B(N)/A(N,N)
   do k=1,N-1
      i=N-k
      B(i)=(B(i)-sum(A(i,i+1:N)*B(i+1:N)))/A(i,i)
   enddo
   det=det*A(N,N)
   X=B
   A=A0
   B=B0
end subroutine lie
