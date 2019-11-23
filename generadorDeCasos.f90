program tpfinal

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
    real, allocatable :: vector(:), matriz(:,:)
    character(len = 20) :: nombreArchivoParabolica, nombreArchivoEliptica
    
    nombreArchivoParabolica = "parabolica1.dat"
    nombreArchivoEliptica = "eliptica1.dat"
    
    allocate(seed(longitud_semilla))
    seed = 314159265
    call random_seed(put=seed)
    call random_seed(get=seed)
    
    allocate(vector(orden))
    allocate(matriz(n,m))
!~     Si orden=100 -> tamanio vector = 100
!~     call armarVector(vector, orden)
    call armarMatriz(matriz, n, m)

!~     call escribeArchivoVector(vector, orden, nombreArchivoParabolica)
    call escribeArchivoMatriz(matriz, n, m, nombreArchivoEliptica)
    
    deallocate(vector)
    deallocate(seed)
end program

subroutine armarVector(vector, n)
    integer i
    integer n
    real x
    real vector(n)
    
    do i=1, n
        CALL random_number(x)
        vector(i) = x*1000
    end do
end subroutine

subroutine armarMatriz(matriz, n, m)
    integer n, m
    real matriz(n,m)
    real x
    
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

subroutine escribeArchivoVector(vector, n, nombreArchivo)
    real vector(n)
    integer n
    character(len = 20) :: nombreArchivo
    
    open (1, file = nombreArchivo, status='replace')
    write(1,*) n
    write(1,*) vector(:)
    close (1, status='keep')
end subroutine

subroutine escribeArchivoMatriz(matriz, n, m, nombreArchivo)
    real matriz(n, m)
    integer n
    character(len = 20) :: nombreArchivo
    
    open (1, file = nombreArchivo, status='replace')
    write(1,*) n, m
    write(1,*) matriz(1,:)
    do i=2, n-1
        write(1,*) matriz(i,:)
!~         Puedo evitar la redundancia si solo escribo los bordes
!~         write(1,*) matriz(i,1)
!~         write(1,*) matriz(i,m)
    end do
    write(1,*) matriz(n,:)
    close (1, status='keep')
end subroutine
