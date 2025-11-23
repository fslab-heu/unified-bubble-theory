subroutine advance_near_field()
   use model
   use DG
   use near_field
   use bubble, only: bubbles
   implicit none
   integer :: count
   logical :: bubble_trigered
   count = 0
   print*, '============================================'
   print*, ''
   print*, 'Starting near field simulation...'
   do while (.true.)
      call update_dt(near_DG_1e)
      call near_DG_1e%march_DG(dt, flux, num_flux, source, bc_left, bc_right)
      call update_near_field(near_DG_1e)
      count = count + 1
      t = t + dt
      if (Rs>1.1*Rb) then
         exit
      end if
   end do
   print*, 'Starting projection to multi-element DG...'
   call project_1e_to_multi(near_DG_1e, near_DG)
   call output_near_field(near_DG)
   t_bubble = 1e10;
   bubble_trigered = .false.
   print*, '============================================'
   print '(4A9)', 'Step', 'Time', 'Rs', 'Pm'
   do while (.true.)
      call update_dt(near_DG)
      call near_DG%march_DG(dt, flux, num_flux, source, bc_left, bc_right)
      call update_near_field(near_DG)
      count = count + 1
      t = t + dt
      if(mod(count,10) == 0) then
         call output_near_field(near_DG)
         if(mod(count,1000) == 0)then
            print*, '--------------------------------------------'
            print '(4A11)', 'Step', 'Time', 'Rs', 'Pm'
            print*, '--------------------------------------------'
         endif
         if(mod(count,100)==0) print"(I11,3E11.3)", count, t, Rs, Ps-P0
      end if
      
      if(Rs>10*Rb.and.mod(count,10)==0)then
         call near_shock2bubble()
      endif
      if ((Rs>10*Rb.or.t-t_arrive>tend).and..not.bubble_trigered) then
         t_bubble = t
         bubble_trigered = .true.
         call bubbles(1)%ini_bubble(Rb, dRb, 1)
      end if
      if(t>2.0*t_bubble)then
         t_far = t
         exit
      endif
   end do
   if(mod(count,100).ne.0) call output_near_field(near_DG)
   !close(101)
   !close(102)
   !close(103)
   print*, 'Near field simulation done.'
end subroutine advance_near_field

subroutine advance_far_field()
   use model
   use far_field
   implicit none
   integer :: i, count
   print*, ''
   print*, '============================================'
   print*, 'Starting far field simulation...'
   print*, '--------------------------------------------'
   print '(4A11)', 'Step', 'Time', 'Rs', 'Pm'
   print*, '--------------------------------------------'
   call ini_far_field()
   call update_dt(far_DG)
   do while (.true.)
      call press2vel()
      call far_DG%march_DG(dt, flux, num_flux, &
         source, bc_left, bc_right)
      call update_far_field(far_DG)
      t = t + dt
      call calculate_shock2bubble()
      if(mod(count,10) == 0) then
         call output_far_field(far_DG)
         if(mod(count,1000) == 0)then
            print*, '--------------------------------------------'
            print '(4A11)', 'Step', 'Time', 'Rs', 'Pm'
            print*, '--------------------------------------------'
         endif
         if(mod(count,100)==0) print"(I11,3E11.3)", count, t, Rs, Ps-P0
      endif
      count = count + 1
      if (t - t_arrive >tend.or.&
         (Rs>50*Rb.and.Rs-L>dist)) then
         exit
      end if
      call update_dt(far_DG)
   end do
   print*, 'Far field simulation done.'

end subroutine advance_far_field

subroutine advance_bubble_dynamics()
   use bubble
   use model
   implicit none
   integer :: i
   real :: pout, t_print
   print*, ''
   print*, 'Starting bubble dynamics simulation...'
   print*, '--------------------------------------------'
   print '(4A11)', 'Step', 'Time', 'Rb', 'Zb'
   print*, '--------------------------------------------'
   t = t_bubble
   inc = 0
   do while (.true.)
      inc = inc + 1
      call collect_dt()
      call bubbles(1)%advance()
      if(mod(inc,10) == 0) then
         write(201,'(3E15.6)') t, bubbles(1)%R, bubbles(1)%center(3)
         t_print = t + (dist - bubbles(1)%R)/c0 - t_arrive
         if(t_print> t_shock_end)then
            pout = rho0*(bubbles(1)%R/dist*&
                (bubbles(1)%H+0.5*bubbles(1)%dR**2.0))
            write(200,'(2E15.6)') t_print, pout + (bubbles(1)%pamb - p0)*bubbles(1)%R/dist
         endif
         if(mod(inc,100) == 0)then
            print*, '--------------------------------------------'
            print '(4A11)', 'Step', 'Time', 'Rb', 'Zb'
            print*, '--------------------------------------------'
         endif
         print"(I11,3E11.3)", inc, t, bubbles(1)%R, bubbles(1)%center(3)
      end if
      t = t + dt
      if(t- t_arrive >= tend -1e-10) then
         exit
      end if
   end do
end subroutine advance_bubble_dynamics
