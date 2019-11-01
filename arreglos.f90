module arreglos

    implicit none

contains

subroutine cargarMatriz(A)
    intent(out) :: A
    integer(4) N, M, i, j
    real(8) A(:, :)

    N = size(A, dim=1)
    M = size(A, dim=2)
    do i = 1, N
        do j = 1, M
            write(*, '(A, I2, A, I2, A)', advance='NO') 'Ingrese el elemento (', i, ', ', j, ') de la matriz: '
            read(*, *) A(i, j)
        end do
    end do
end subroutine cargarMatriz

subroutine mostrarMatriz(A, formato)
    real(8), intent(in) :: A(:, :)
    character(len=*), optional :: formato
    integer(4) i, filas

    filas = size(A, dim=1)
    if (present(formato)) then
        do i = 1, filas
            write(*, formato) A(i, :)
        end do
    else
        do i = 1, filas
            write(*, *) A(i, :)
        end do
    end if
end subroutine mostrarMatriz

subroutine leerMatriz(A, archivo)
    character(len=*), intent(in) :: archivo
    real(8), intent(out) :: A(:, :)
    integer(4) N, M, i

    N = size(A, dim=1)
    M = size(A, dim=2)
    open(unit=2, file=archivo, access='SEQUENTIAL')
    do i = 1, N
        read(2, *) A(i, :M)
    end do
    close(2)
end subroutine leerMatriz

subroutine grabarMatriz(A, archivo)
    intent(in) :: A, archivo
    integer(4) i, filas
    real(8) A(:, :)
    character(len=*) :: archivo

    filas = size(A, dim=1)
    open(unit=2, file=archivo, access='SEQUENTIAL', status='REPLACE')
    do i = 1, filas
        write(2, *) A(i, :)
    end do
    close(2, status='KEEP')
end subroutine grabarMatriz

subroutine cargarVector(V)
    intent(out) :: V
    real(8) V(:)
    integer(4) N, i
    
    N = size(V)
    do i = 1, N
        write(*, '(A, I3, A)', advance='NO') 'Ingrese el elemento ', i, ' del vector: '
        read(*, *) V(i)
    end do
end subroutine cargarVector

subroutine mostrarVector(V, formato)
    real(8), intent(in) :: v(:)
    character(len=*), optional :: formato

    if (present(formato)) then
        write(*, formato) V
    else
        write(*, *) V
    end if
end subroutine mostrarVector

subroutine leerVector(V, archivo)
    character(len=*), intent(in) :: archivo
    real(8), intent(out) :: V(:)
    integer(4) N
    
    N = size(V)
    open(unit=2, file=archivo, access='SEQUENTIAL')
    read(2, *) V(:N)
    close(2)
end subroutine leerVector

subroutine grabarVector(V, archivo)
    intent(in) :: V, archivo
    real(8) V(:)
    character(len=*) :: archivo
    
    open(unit=2, file=archivo, access='SEQUENTIAL', status='REPLACE')
    write(2, *) V
    close(2, status='KEEP')
end subroutine grabarVector

end module arreglos
