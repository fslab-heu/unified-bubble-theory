      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
      !!!!! CopyRight by Liu Yunlong   Email:yunlong_liu@hrbeu.edu.cn !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
      module math
      implicit none
      real, parameter,private:: pi=3.141592653589793
      interface rotate_vect !! 旋转矢量
        module procedure rotate_vect
        module procedure rotate_vect2
      end interface   
      interface resort
        module procedure resort
        module procedure resort_i
      end interface
      
    contains
    subroutine rotate_vect_HW(input,omega)
    implicit none
    real input(3),omega(3)
    real D,Sij(3,3),Rij(3,3),Eij(3,3)
    integer i,j,k
    real L
    L=length(input)
    if(L<1e-10)return
    Eij=0;Eij(1,1)=1.0;Eij(2,2)=1.0;Eij(3,3)=1.0
    D = 2.0+0.5*sum(omega**2)
    Sij(1,:)=[0.0,-omega(3),omega(2)]
    Sij(2,:)=[omega(3),0.0,-omega(1)]
    Sij(3,:)=[-omega(2),omega(1),0.0]
    
    ! 1st order method
    Rij=Eij+Sij!(2.0*Eij+)/2.0/D
    
    input=matmul(Rij,input)
    input=input/length(input)*L
    
    return
    end subroutine
    
    
      integer function permut(input)
      ! permutations tensor of the input
      ! odd num, -1
      ! even num, 1
      ! otherwise, 0
      implicit none
      integer input(:)
      integer N,i,j,m,sw
      integer,allocatable:: temp(:)
      N=size(input)
      allocate(temp(N))
      temp=input
      m=0
      do j=1,N-1
          if(m==-1)then
              exit
          endif
          
          do i=1,N-j
            if(temp(i)>temp(i+1))then
                m=m+1
                sw=temp(i)
                temp(i)=temp(i+1)
                temp(i+1)=sw
            elseif(temp(i)==temp(i+1))then
                m=-1
                exit
            endif
          enddo
      enddo
      if(m==-1)then
          permut=0
      elseif(mod(m,2)==0)then
          permut=1
      else
          permut=-1
      endif
      deallocate(temp)
      end function
      function tensor2vect(A)
      real A(3,3)
      real tensor2vect(6)
      tensor2vect(1)=A(1,1)
      tensor2vect(2)=A(2,2)
      tensor2vect(3)=A(3,3)
      tensor2vect(4)=A(1,2)
      tensor2vect(5)=A(1,3)
      tensor2vect(6)=A(2,3)
      end function
      
      function vect2tensor(A)
      real A(6)
      real vect2tensor(3,3)
      vect2tensor(1,1)=A(1)
      vect2tensor(2,2)=A(2)
      vect2tensor(3,3)=A(3)
      vect2tensor(1,2)=A(4)
      vect2tensor(1,3)=A(5)
      vect2tensor(2,3)=A(6)
      vect2tensor(2,1)=vect2tensor(1,2)
      vect2tensor(3,1)=vect2tensor(1,3)
      vect2tensor(3,2)=vect2tensor(2,3)
      end function
      
      function stress_spin(stress,Rmat)
      real stress(6),Rmat(3,3)
      real stress_spin(6)
      real tensor(3,3)
      tensor=vect2tensor(stress)
      stress_spin=stress_spin+tensor2vect(&
          matmul(tensor,Rmat)+matmul(Rmat,tensor))
      
      end function
      
      function dist2line(node1,node2,pos)
      implicit none
      real node1(3),node2(3),pos(3)
      real dist2line
      real dir(3),sub(3)
      dir = node2-node1
      dir = dir/length(dir)
      sub = pos-node1
      dist2line = length(cross(sub,dir))
      return
      endfunction
      
      function permutation(i,j,k)
        integer i,j,k
        real permutation
        if(i==j.or.i==k.or.j==k)then
            permutation=0
        elseif(i+j+k==2.or.i+j+k==-1)then
            permutation=1
        else
            permutation=-1
        endif
        return
      end function permutation
      
      
      
      function rotate_vect(in,axis0)
        real in(:),axis0(:),rotate_vect(3)
        real ang
        integer i,j,k,N
        real e1(3),e2(3),R,dis(3),axis(3)
        rotate_vect=in
        if(size(axis0)==3)then
            axis=axis0
            ang=length(axis)
            ang=mod(ang,2*pi);
            axis=axis/ang;
        else
            axis=axis0(1:3)
            ang=axis0(4)
            ang=mod(ang,2*pi);
        endif
        
        
        if (abs(ang)<1e-8)then
            return;
        endif
        
        
        

        e1=cross(in,axis);
        e2=cross(axis,e1);

        if (length(e1)<1e-8)then
            return;
        endif
        e1=e1/length(e1);
        e2=e2/length(e2);
        R=length(cross(in,axis));
        dis=e1*R*sin(ang)+e2*(R-R*cos(ang));
        rotate_vect=-dis+in;
      end function rotate_vect
      
      function rotate_vect2(in,axis0)
        real in(:,:),axis0(:)
        real,allocatable::rotate_vect2(:,:)
        real ang,axis(3)
        integer i,j,k,N
        real e1(3),e2(3),R,dis(3)
        axis=axis0
        N=size(in(:,1))
        allocate(rotate_vect2(N,3))
        ang=length(axis)
        rotate_vect2=in
        if (ang<1e-8)then
            return;
        endif
        ang=mod(ang,2*pi);
        axis=axis/ang;
        do i=1,N
            e1=cross(in(i,:),axis);
            e2=cross(axis,e1);
            if (length(e1)<1e-8)then
                cycle
            endif
            e1=e1/length(e1);
            e2=e2/length(e2);
            R=length(cross(in(i,:),axis));
            dis=e1*R*sin(ang)+e2*(R-R*cos(ang));
            rotate_vect2(i,:)=-dis+in(i,:);
        enddo
      end function rotate_vect2
      subroutine reflect(in,out,plane)
!!DEC$ ATTRIBUTES DLLEXPORT :: reflect
        integer i,j,k,N
        real plane(2,3)
        real in(:,:),out(:,:)
        
        N=size(in(:,1))
        plane(2,:)=plane(2,:)/length(plane(2,:))
        do i=1,N
            out(i,:)=in(i,:)-2*dot_product(in(i,:)-plane(1,:),&
            plane(2,:))*plane(2,:)
        enddo
        return
      end subroutine reflect
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dis2plane(in,dis,plane)
!!DEC$ ATTRIBUTES DLLEXPORT :: dis2plane
        integer i,j,N
        real dis(:),in(:,:),plane(2,3)
        N=size(in(:,1))
        plane(2,:)=plane(2,:)/length(plane(2,:))
        do i=1,N
            dis(i)=dot_product(in(i,:)-plane(1,:),&
            plane(2,:))
        
        enddo
        return
      end subroutine dis2plane
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine resort(in,ind,N)
	  integer N
	  real in(N),temp(N)
	  integer ind(N),i
	  do i=1,N
	      temp(i)=in(ind(i))
	  enddo
	  in=temp
	  return
	end subroutine resort
	subroutine resort_i(in,ind,N)
	  integer N
	  integer in(N),temp(N)
	  integer ind(N),i
	  do i=1,N
	      temp(i)=in(ind(i))
	  enddo
	  in=temp
	  return
	end subroutine	 resort_i	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 微分
      function diff(x)
!!DEC$ ATTRIBUTES DLLEXPORT :: diff
          implicit none
          real x(:)
          integer N
          real,allocatable:: diff(:)
          N=size(x)
          allocate(diff(N-1))
          diff=x(2:N)-x(1:N-1)
      end function diff     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!c叉乘
        function cross(a,b)
!!DEC$ ATTRIBUTES DLLEXPORT :: cross
            implicit none
            real a(3),b(3),cross(3)
            cross(1)=a(2)*b(3)-b(2)*a(3)
            cross(2)=a(3)*b(1)-b(3)*a(1)
            cross(3)=a(1)*b(2)-b(1)*a(2)
            return
        end function cross
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!c模长(任意长度向量)
          function lengthv(a,m)
            integer N,i,j,m
            real a(:,:),lengthv(m)
            N=size(a(1,:))
            lengthv=0
            do j=1,m
                do i=1,N
                    lengthv(j)=lengthv(j)+a(j,i)**2d0
                enddo
            enddo
            lengthv=sqrt(lengthv)
            return
          end function lengthv
          
          function length(a)
            integer N,i
            real a(:),length
            N=size(a)
            length=0
            do i=1,N
                length=length+a(i)**2d0
            enddo
            length=sqrt(length)
            return
          end function length
          
          function coor_tranf(vect0,coor_sys)
!!DEC$ ATTRIBUTES DLLEXPORT :: coor_tranf    
          !!coor_sys 三行分别为 原点坐标 z轴方向 x轴方向
            use linear_equation,only: lie
            implicit none
            real vect0(:),coor_sys(:,:)
            real vect_len,A(3,3)
            real coor_tranf(3),vect(3)
            
            integer i,j,k
            vect=vect0
             A(:,3)=coor_sys(2,:)
            A(:,2)=cross(coor_sys(2,:),coor_sys(3,:))
            A(:,1)=cross(A(:,2),A(:,3))
            do i=1,3  !!! 单位化
                A(:,i)=A(:,i)/length(A(:,i))
            enddo
            vect(:)=vect(:)-coor_sys(1,:)
            
            call lie(A,coor_tranf,vect)
            return
          end function coor_tranf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!cc
          function coor_tranf_inv(vect0,coor_sys)
!!DEC$ ATTRIBUTES DLLEXPORT :: coor_tranf_inv    
          !!coor_sys 三行分别为 原点坐标 z轴方向 x轴方向
            use linear_equation
            implicit none
            real vect0(:),coor_sys(:,:)
            real vect_len,A(3,3)
            real coor_tranf_inv(3),vect(3)
            integer i,j,k
            vect=vect0
            A(:,3)=coor_sys(2,:)
            A(:,2)=cross(coor_sys(2,:),coor_sys(3,:))
            A(:,1)=cross(A(:,2),A(:,3))
            do i=1,3  !!! 单位化
                A(:,i)=A(:,i)/length(A(:,i))
            enddo
            
            coor_tranf_inv=matmul(A,vect)
            coor_tranf_inv=coor_tranf_inv+coor_sys(1,:)
            return
          end function coor_tranf_inv
!          subroutine coor_trans(vect_out,vect_in,coor_sys)
!!!DEC$ ATTRIBUTES DLLEXPORT :: coor_trans
!          !!coor_sys 三行分别为 原点坐标 z轴方向 x轴方向
!            use linear_equation
!            implicit none
!            real vect_out(:,:),vect_in(:,:),coor_sys(:,:)
!            real vect_len,A(3,3)
!            real vect(3)
!            
!            integer i,j,k,N
!            N=size(vect_in(:,1))
!            A(:,3)=coor_sys(2,:)
!            A(:,2)=cross(coor_sys(2,:),coor_sys(3,:))
!            A(:,1)=cross(A(:,2),A(:,3))
!            do i=1,3  !!! 单位化
!                A(:,i)=A(:,i)/length(A(:,i))
!            enddo
!            do j=1,N
!                vect=vect_in(j,:)
!                vect(:)=vect(:)-coor_sys(1,:)
!                call lie(A(:,:),vect_out(j,:),vect(:))
!            enddo
!            return
!          end subroutine coor_trans
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!cc
          subroutine coor_trans_inv(vect_out,vect_in,coor_sys)
!!DEC$ ATTRIBUTES DLLEXPORT :: coor_trans_inv
          !!coor_sys 三行分别为 原点坐标 z轴方向 x轴方向
            use linear_equation
            implicit none
            real vect_out(:,:),vect_in(:,:),coor_sys(:,:)
            real vect_len,A(3,3)
            real vect(3)
            integer i,j,k,N
            N=size(vect_in(:,1))
            A(:,3)=coor_sys(2,:)
            A(:,2)=cross(coor_sys(2,:),coor_sys(3,:))
            A(:,1)=cross(A(:,2),A(:,3))
            do i=1,3  !!! 单位化
                A(:,i)=A(:,i)/length(A(:,i))
            enddo
            do i=1,N
                vect=vect_in(i,:)
                vect_out(i,:)=matmul(A,vect)
                vect_out(i,:)=vect_out(i,:)+coor_sys(1,:)
            enddo
            return
          end subroutine coor_trans_inv

          function interp1(x,y,x0,interp_type)
!!DEC$ ATTRIBUTES DLLEXPORT :: interp_dp 
            real x(:),y(:),x0
            integer N,NN,i,j,k,itype,L,M
            character(*) interp_type
            real interp1
            N=size(x)
            if (interp_type=='linear'.or.interp_type=='Linear'.or.&
                interp_type=='LINEAR')then
                itype=1
            else
                itype=2
            endif
            
            if ((itype==1.and.N<2).or.(itype==2.and.N<3))then
                write(*,*) 'ERRO: interpotion erro for insurfficient &
                 input data'
                stop
            endif
            
            if (itype==1)then
                    if (x0<x(1).or.x0>x(N))then
                        write(*,*) 'ERRO: input x0 out of range'
                        stop
                    endif
                    do j=1,N
                        if(x0<x(j)) exit
                    enddo
                    L=j-1;M=j
                  interp1=y(L)+(x0-x(L))/(x(M)-x(L))*(y(M)-y(L))
            else
            
                write(*,*) 'Please waite for new version'
                stop
            endif
            return
        end function interp1


	  function Newton(Func,x0)
!!DEC$ ATTRIBUTES DLLEXPORT :: Newton 
			real Newton,x,err,x0,y,dy,dx
			real, External:: Func
			integer count
			err=1
			x=x0
			dx=1e-4
			count=0
			do while (err>1e-6.and.count<1000)
				y=Func(x)
				dy=(Func(x+dx)-Func(x-dx))/dx/2.0
				err=y/dy
				x=x-err
				err=abs(err)
			enddo
			Newton=x
	  end function Newton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        

        subroutine sort(A,ind,asc)
!!DEC$ ATTRIBUTES DLLEXPORT :: sort         
!!!!!!!!!!!!!!! 冒泡法
            real A(:)
            integer ind(:)
            integer,optional:: asc ! 0 for increasing , 1 for decreasing
            integer i,j,k,N,asc0
            real temp
            N=size(A)
            if(present(asc))then
                asc0=asc
            else
                asc0 = 0
            endif
            
            do i=1,N
                ind(i)=i
            enddo
            if(asc0==0)then
                do i=1,N-1
                    do j=1,N-i
                        if (A(j)>A(j+1))then
                            temp=A(j)
                            A(j)=A(j+1)
                            A(j+1)=temp
                            k=ind(j)
                            ind(j)=ind(j+1)
                            ind(j+1)=k
                        endif
                    enddo
                enddo
            else
                do i=1,N-1
                    do j=1,N-i
                        if (A(j)<A(j+1))then
                            temp=A(j)
                            A(j)=A(j+1)
                            A(j+1)=temp
                            k=ind(j)
                            ind(j)=ind(j+1)
                            ind(j+1)=k
                        endif
                    enddo
                enddo
            endif
            
        end subroutine sort
        subroutine sort_int(A,ind,asc)
!!DEC$ ATTRIBUTES DLLEXPORT :: sort_int     
!!!!!!!!!!!!!!! 冒泡法
            integer A(:)
            integer ind(:)
            integer,optional:: asc ! 0 for increasing , 1 for decreasing
            integer i,j,k,N,asc0
            real temp
            N=size(A)
            if(present(asc))then
                asc0=asc
            else
                asc0 = 0
            endif
            
            do i=1,N
                ind(i)=i
            enddo
            if(asc0==0)then
                do i=1,N-1
                    do j=1,N-i
                        if (A(j)>A(j+1))then
                            temp=A(j)
                            A(j)=A(j+1)
                            A(j+1)=temp
                            k=ind(j)
                            ind(j)=ind(j+1)
                            ind(j+1)=k
                        endif
                    enddo
                enddo
            else
                do i=1,N-1
                    do j=1,N-i
                        if (A(j)<A(j+1))then
                            temp=A(j)
                            A(j)=A(j+1)
                            A(j+1)=temp
                            k=ind(j)
                            ind(j)=ind(j+1)
                            ind(j+1)=k
                        endif
                    enddo
                enddo
            endif
        end subroutine sort_int
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
      subroutine car2cyl(pin,pout)
!!DEC$ ATTRIBUTES DLLEXPORT :: car2cyl    
        real pin(:,:),pout(:,:)
        integer N,i,j
        real R
        N=size(pin(:,1))
        do i=1,N
            pout(i,3)=pin(i,3)
            R=length(pin(i,1:2))
            pout(i,1)=R
            if (R.ne.0)then
                pout(i,2)=acos(pin(i,1)/R)
            else
                pout(i,2)=0
            endif
            if (pin(i,2)<0) pout(i,2)=2*real(pi)-pout(i,2)
        enddo
        return
      end subroutine car2cyl
      subroutine cyl2car(pin,pout)
!!DEC$ ATTRIBUTES DLLEXPORT :: cyl2car    
        real pin(:,:),pout(:,:)
        integer N,i,j
        real R,xita
        N=size(pin(:,1))
        do i=1,N
            pout(i,3)=pin(i,3)
            R=pin(i,1)
            xita=pin(i,2)
            pout(i,1)=R*cos(xita)
            pout(i,2)=R*sin(xita)
        enddo
        return
      end subroutine cyl2car
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      subroutine rotate(in,out,axis,ang)
!!!DEC$ ATTRIBUTES DLLEXPORT :: rotate   
!        real in(:,:),out(:,:),axis(2,3),ang
!        integer N,i,j
!        real trans(3,3),v1(3),v2(3)
!        
!        real ,allocatable:: mid1(:,:),mid2(:,:)
!        data v1/1,0,0/
!        data v2/0,1,0/
!        N=size(in(:,1))
!           
!        trans(1,:)=axis(1,:)
!        trans(2,:)=axis(2,:)-axis(1,:)
!        if (abs(dot_product(axis(2,:)-axis(1,:),v1))>abs(dot_product&
!         (axis(2,:)-axis(1,:),v2)))then
!            trans(3,:)=cross(axis(2,:)-axis(1,:),v2)
!        else
!            trans(3,:)=cross(axis(2,:)-axis(1,:),v1)
!        endif
!        
!        allocate(mid1(N,3),mid2(N,3))
!
!        do i=1,N
!            mid1(i,:)=coor_tranf(in(i,:),trans)
!        enddo
!        call car2cyl(mid1,mid2)
!        mid2(:,2)=mid2(:,2)+ang
!        call cyl2car(mid2,mid1)
!        do i=1,N
!            out(i,:)=coor_tranf_inv(mid1(i,:),trans)
!        enddo
!        deallocate(mid1,mid2)
!      return
!      
!      
!      end subroutine rotate

  
      end module