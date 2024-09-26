      program hello_world
        use AMP_INTERFACE
        use iso_c_binding, only: C_INT, C_PTR
        integer(c_int) :: argc
        type(c_ptr) :: argv

        argc = 1
        
        ! startup
        call startup( argc, argv )
        
        ! test initialized

        ! shutdown
        call shutdown()
        ! Finished
        print *, 'Finished'
      end program


