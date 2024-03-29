    subroutine pp_limiter()
    !! positivity preserving limiter
    use shockwave
    use fsintegration
    use polynomial
    real:: ukb(ndg,np),uukb(ndg,np),Eukb(ndg,np)
    real :: uin(3,np,ndg)
    integer :: i,j,k,h
    real :: phi_off,R
    real :: eps,rho_min,theta1,ein,pij,theta2,xb(ndg+1),plb,prb,einlb,einrb,ulb,urb
    real :: dx,flee(np,np,2),fleedge(np,2)
    real:: ug(3,np),ta(np),UL(3),UR(3),Clba,Crba,usb,xeb(np),xg(np,2)
    real::m(3),V,ukgb(np),uukgb(np),Eukgb(np),E0b(3),eing(np),Cbg(np),pbg(np)

    call gl0%gl_init_legendre(np)
    
    xg(1:np,1)=gl0%Pos
    xg(1:np,2)=gl0%weight 

    dx=1.0/ndg
    eps=1.0e-13

    do i=1,ndg,1
      do k=1,np,1
         uin(:,k,i)=ans_dg(k,:,i)
      end do
    end do
    
    do i=1,np,1
      do k=1,np,1
        flee(i,k,1)=basis((xg(k,1)+1)/2.0,i)
      end do
    end do

    do i=1,np,1
      do k=1,np,1
        flee(i,k,2)=dbasis((xg(k,1)+1)/2.0,i)/dx
      end do
    end do

    do i=1,ndg,1
       xb(i+1)=xb(i)+1.0/ndg
    end do

    do j=2,ndg,1
        UL = 0
        UR = 0
        do i=1,np
            UL = UL+ ans_dg(i,:,j-1)*basis(1.0,i)
            UR = UR+ ans_dg(i,:,j)*basis(0.0,i)
        enddo


        ulb=UL(2)/UL(1)-xb(j)*(dRb)-dRb
        urb=UR(2)/UR(1)-xb(j)*(dRb)-dRb
        einlb=UL(3)/UL(1)-0.5*(UL(2)/UL(1))**2
        einrb=UR(3)/UR(1)-0.5*(UR(2)/UR(1))**2
        plb = e2p(einlb*UL(1))
        prb = e2p(einrb*Ur(1))
        Clba=soundspeed(Plb,UL(1))
        Crba=soundspeed(Prb,UR(1))

        do h=1,np,1
            xeb(h) = xb(j)+dx*(xg(h,1)+1)/2.0
        end do

        ug=matmul(uin(:,:,j),flee(:,:,1))
        eing=ug(3,:)/ug(1,:)-0.5*(ug(2,:)/ug(1,:))**2
        pbg=ug(1,:)*eing*(gamma-1.0)-gamma*Pw
        cbg=((pbg+Pw)/ug(1,:)*gamma)
      
        m=0
        V=0
  
        do k=1,np,1
            m(:)=m(:)+(xeb(k)*Rb+Rb)**2*ug(:,k)*xg(k,2)*dx*Rb !dx
        end do
  
        do k=1,np,1
            V=V+(xeb(k)*Rb+Rb)**2*xg(k,2)*dx*Rb !dx
        end do

  
  
        if (j==2.and.(Clba < eps.or.Crba<eps.or.UL(1)<eps.or.UR(1)<eps)) then
                uin(:,2:np,1) = 0.0
                uin(:,1,1) = m(:)/V
                do k=1,np,1
                    ans_dg(k,:,1)=uin(:,k,1)
                end do
        end if

        if (j==2.and.(plb < eps.or.prb < eps)) then 
                uin(:,2:np,1) = 0.0
                uin(:,1,1) = m(:)/V
                do k=1,np,1
                    ans_dg(k,:,1)=uin(:,k,1)
                end do
        end if 

        do h=1,np,1
            if (j==2.and.(ug(1,h)<eps))then
                uin(:,2:np,1) = 0.0
                uin(:,1,1) = m(:)/V
                do k=1,np,1
                    ans_dg(k,:,1)=uin(:,k,1)
                end do
            end if
        end do
     
        do h=1,np,1
            if (j==2.and.(cbg(h)<eps.or.pbg(h)<eps))then
                uin(:,2:np,1) = 0.0
                uin(:,1,1) = m(:)/V
                do k=1,np,1
                    ans_dg(k,:,1)=uin(:,k,1)
                end do
            end if
        end do
        if (Clba < eps.or.Crba<eps.or.UL(1)<eps.or.UR(1)<eps) then
            uin(:,2:np,j) = 0.0
            uin(:,1,j) = m(:)/V
        end if
    
        if (plb < eps.or.prb < eps) then 
            uin(:,2:np,j) = 0.0
            uin(:,1,j) = m(:)/V
        end if 

        do h=1,np,1
            if (ug(1,h)<eps)then
                uin(:,2:np,j) = 0.0
                uin(:,1,j) = m(:)/V
            end if
        end do
        do h=1,np,1
            if (cbg(h)<eps.or.pbg(h)<eps)then
                uin(:,2:np,j) = 0.0
                uin(:,1,j) = m(:)/V
            end if
        end do
        do k=1,np,1
          ans_dg(k,:,j)=uin(:,k,j)
        end do
    end do

    end subroutine