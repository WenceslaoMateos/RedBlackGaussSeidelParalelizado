program tpfinal
    use moduloArchivos

!~ Generador de casos de prueba
!~ Si se quiere generar datos para parabolicas, modificar 'orden' y 'nombreArchivoParabolica'
!~ Si se quiere generar datos para elipticas, modificar 'n', 'm' y 'nombreArchivoEliptica'
!~ Notas:
!~    Siempre se genera el mismo conjunto de numeros (misma semilla)
!~    Se puede descomentar la linea 95 para evitar redundancia en elipticas
!~    Dejar comentado armarVector o armarMatriz (que solo uno se ejecute)
    INTEGER, PARAMETER :: longitud_semilla = 33
    integer, parameter :: orden = 100
    integer, parameter :: n = 4, m = 4
    integer, allocatable :: seed(:)
    integer longitud, filas, columnas
    real(8), allocatable :: vector(:), matriz(:,:)
    character(len = 20) :: nombreArchivoParabolica, nombreArchivoEliptica
    
    nombreArchivoParabolica = "parabolica1.dat"
    nombreArchivoEliptica = "eliptica1.dat"
    
    allocate(seed(longitud_semilla))
    seed = 314159265
    call random_seed(put=seed)
    call random_seed(get=seed)
    
    allocate(vector(orden))
!~     allocate(matriz(n,m))
!~     Si orden=100 -> tamanio vector = 100
!~     call armarVector(vector, orden)
!~     call armarMatriz(matriz, n, m)

    call armarVector(vector, orden)
    call escribeArchivoVector(vector, orden, 'Datos/term.dat')
    call armarVector(vector, orden)
    call escribeArchivoVector(vector, orden, 'Datos/xini.dat')
!~     call escribeArchivoMatriz(matriz, n, m, nombreArchivoEliptica)
    
!~     call leeArchivoVector(vector, longitud, nombreArchivoParabolica)
!~     call leeArchivoMatriz(matriz, filas, columnas, nombreArchivoEliptica)
    
    write(*,*) matriz
    
    if (allocated(vector)) then
        deallocate(vector)
    end if
    deallocate(seed)
end program

subroutine armarVector(vector, n)
    integer i
    integer n
    real(8) x
    real(8) vector(n)
    
    do i=1, n
        CALL random_number(x)
        vector(i) = x*1000
    end do
end subroutine

subroutine armarMatriz(matriz, n, m)
    integer n, m
    real(8) matriz(n,m)
    real(8) x
    
    matriz=0
    do i=1, m
        CALL random_number(x)
        matriz(1,i) = x
    enddo
    do i=2, n-1
        CALL random_number(x)
        matriz(i,1) = x
        CALL random_number(x)
        matriz(i,m) = x
    enddo
    do i=1, m
        CALL random_number(x)
        matriz(n,i) = x
    enddo
end subroutine
