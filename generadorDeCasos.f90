module casosDinamicos
    implicit none

    contains
    !~  Generador de casos de prueba
    !~  Siempre se genera el mismo conjunto de numeros (misma semilla)

    function generaCaso(orden)
        integer(4), intent(in) :: orden
        integer(4), parameter :: longitud_semilla = 33
        integer(4), allocatable :: seed(:)
        real(8) generaCaso(orden)
        
        allocate(seed(longitud_semilla))
        seed = 314159265
        call random_seed(put=seed)
        call random_seed(get=seed)
        
        call armarVector(generaCaso, orden)
        
        deallocate(seed)
    end function generaCaso

    subroutine armarVector(vector, n)
        integer(4) i, n
        real(8) x, vector(n)
        
        do i = 1, n
            call random_number(x)
            vector(i) = x*1000
        end do
    end subroutine
    
end module casosDinamicos