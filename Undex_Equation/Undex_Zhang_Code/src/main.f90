    program main
        implicit none
        ! initialize variables
        call initialize()
        ! solve near-field shock wave
        call advance_near_field()
        ! solve far-field shock wave
        call advance_far_field()
        ! solve bubble dynamics
        call advance_bubble_dynamics()
        ! terminate program
        print*, ''
        print*, '============================================'
        print*, 'Simulation completed successfully.'
        print*, ''
        print*, 'Please find the results in the output folder.'
        print*, ''
        print*, '============================================'
    end program main