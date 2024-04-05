    module linear_equation
    
    contains
    	subroutine inv_lie(A0,X,tp)
!---------------------列主元消去法  tp=1表示求解线性方程组，tp=0表示求逆
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
    
	subroutine lie(A,X,B)
!!DEC$ ATTRIBUTES DLLEXPORT :: lie
!---------------------列主元消去法
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
    end module