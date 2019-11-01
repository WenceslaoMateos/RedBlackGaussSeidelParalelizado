!para compilar el modulo: gfortran -g -c modulo.f90
module SELs

    implicit none

    abstract interface
        function metodoDirecto(matriz, term_ind)
            real(8), dimension(:, :), intent(in) :: matriz, term_ind
            real(8) metodoDirecto(size(matriz, dim=1), size(matriz, dim=2) + size(term_ind, dim=2))
        end function
    end interface

    abstract interface
        function solucionDirecto(matriz, term_ind)
            real(8), dimension(:, :), intent(in) :: matriz, term_ind
            real(8) solucionDirecto(size(term_ind, dim=1), size(term_ind, dim=2))
        end function
    end interface
    
    abstract interface
        function mNorma(matriz)
            real(8), dimension(:, :), intent(in) :: matriz
            real(8) mNorma
        end function
    end interface
contains

function matrizAmpliada(matriz, term_ind)
    real(8), dimension(:, :), intent(in) :: matriz, term_ind
    real(8) matrizAmpliada(size(matriz, dim=1), size(matriz, dim=2) + size(term_ind, dim=2))
    integer(4) col_mat, col_amp

    col_mat = size(matriz, dim=2)
    col_amp = size(matriz, dim=2) + size(term_ind, dim=2)

    matrizAmpliada(:, :col_mat) = matriz
    matrizAmpliada(:, col_mat + 1:col_amp) = term_ind                
end function matrizAmpliada

subroutine pivotear(matriz, term_ind)
    real(8), dimension(:, :), intent(inout) :: matriz, term_ind
    real(8) auxMatriz(size(matriz, dim=2)), auxTerm(size(term_ind, dim=2))
    integer(4) j, i, filas, columnas, pivote

    filas = size(matriz, dim=1)
    columnas = size(matriz, dim=2)
    j = 1 ! j pensarlo como la diagonal
    do while (j < filas .and. j < columnas)
        pivote = j
        do i = j + 1, filas
            if (abs(matriz(pivote, j)) < abs(matriz (i, j))) then
                pivote = i
            end if
        end do
        if (pivote /= j) then
            auxMatriz = matriz(j, :)
            matriz(j, :) = matriz(pivote, :)
            matriz(pivote, :) = auxMatriz

            auxTerm = term_ind(j, :)
            term_ind(j, :) = term_ind(pivote, :)
            term_ind(pivote, :) = auxTerm
        end if
        j = j + 1
    end do
end subroutine pivotear

!!!!!!!!!!!!!!!!!!!!!!!!!!! Metodos directos !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function gauss(matriz, term_ind)
    real(8), dimension(:, :), intent(in) :: matriz, term_ind
    real(8) gauss(size(matriz, dim=1), size(matriz, dim=2) + size(term_ind, dim=2))
    integer(4) i, j, filas, columnas
    
    gauss = matrizAmpliada(matriz, term_ind)
    filas = size(gauss, dim=1)
    columnas = size(gauss, dim=2)
    do j = 1, columnas - 1
        do i = j + 1, filas
            gauss(i, :) = gauss(i, :) - gauss(j, :) * gauss(i, j)  / gauss(j, j)
            gauss(i, j) = 0.
        end do
    end do 
end function gauss

function solucionGauss(matriz, term_ind)
    real(8), dimension(:, :), intent(in) :: matriz, term_ind
    real(8) sol(size(term_ind, dim=1), size(matriz, dim=2) + size(term_ind, dim=2))
    real(8), dimension(size(term_ind, dim=1), size(term_ind, dim=2)) :: solucionGauss
    real(8), dimension(size(term_ind, dim=2)) :: sumatorias
    integer(4) i, k, filas, colterm

    sol = gauss(matriz, term_ind)
    filas = size(matriz, dim=1)
    colterm = size(matriz, dim=2) + 1
    solucionGauss(filas, :) = sol(filas, colterm:) / sol(filas, filas)
    do i = filas - 1, 1, -1
        k = i + 1
        sumatorias = matmul(sol(i, k:colterm - 1), solucionGauss(k:, :))
        solucionGauss(i, :) = (sol(i, colterm:) - sumatorias) / sol(i, i)
    end do
end function solucionGauss

function gaussJordan(matriz, term_ind)
    real(8), dimension(:, :), intent(in) :: matriz, term_ind
    real(8) gaussJordan(size(matriz, dim=1), size(matriz, dim=2) + size(term_ind, dim=2))
    integer(4) i, j, filas, columnas
    
    gaussJordan = gauss(matriz, term_ind)

    filas = size(matriz, dim=1)
    columnas = size(matriz, dim=2)
    do j = 2, columnas
        do i = 1, j - 1
            gaussJordan(i, :) = gaussJordan(i, :) - gaussJordan(j, :) * gaussJordan(i, j)  / gaussJordan(j, j)
            gaussJordan(i, j) = 0.
        end do
    end do 
end function gaussJordan

function solucionGaussJordan(matriz, term_ind)
    real(8), dimension(:, :), intent(in) :: matriz, term_ind
    real(8) sol(size(term_ind, dim=1), size(matriz, dim=2) + size(term_ind, dim=2))
    real(8) solucionGaussJordan(size(term_ind, dim=1), size(term_ind, dim=2))
    integer(4) i, filas

    sol = gaussJordan(matriz, term_ind)
    filas = size(matriz, dim=1)
    do i = 1, filas
        sol(i, :) = sol(i, :) / sol(i, i)
    end do
    solucionGaussJordan(:, :) = sol(:, size(matriz, dim=2) + 1:)
end function solucionGaussJordan

function reduccionCrout(matriz)
    intent(in) :: matriz
    real(8) matriz(:, :), reduccionCrout(size(matriz, dim=1), size(matriz, dim=2))
    integer(4) orden, i, fila, k, col
 
    reduccionCrout = matriz
    reduccionCrout(1, 2:) = reduccionCrout(1, 2:) / reduccionCrout(1, 1)
    orden = size(matriz, dim=1)
    do i = 2, orden
        ! Calculo de la columna i de L
        do fila = i, orden
            do k = 1, i - 1
                reduccionCrout(fila, i) = reduccionCrout(fila, i) - reduccionCrout(fila, k) * reduccionCrout(k, i)
            end do
        end do
 
        ! Calculo de la fila i de U
        do col = i + 1, orden
            do k = 1, i - 1
                reduccionCrout(i, col) = reduccionCrout(i, col) - reduccionCrout(i, k) * reduccionCrout(k, col)
            end do
            reduccionCrout(i, col) = reduccionCrout(i, col) / reduccionCrout(i, i)
        end do
    end do
 end function reduccionCrout
 
 function solucionCrout(matriz, term_ind)
    real(8), dimension(:, :), intent(in) :: matriz, term_ind
    real(8), dimension(size(term_ind, dim=1), size(term_ind, dim=2)) :: solucionCrout, c
    real(8) aux(size(matriz, dim=1), size(matriz, dim=2))
    integer(4) i, orden, k
 
    orden = size(matriz, dim=1)
    aux = reduccionCrout(matriz)
 
    ! Calcula c
    c(1, :) = term_ind(1, :) / aux(1, 1)
    do i = 2, orden
        c(i, :) = term_ind(i, :)
        do k = 1, i - 1
            c(i, :) = c(i, :) - aux(i, k) * c(k, :)
        end do
        c(i, :) = c(i, :) / aux(i, i)
    end do
 
    ! Calcula la Solución
    solucionCrout(orden, :) = c(orden, :)
    do i = orden -1, 1, -1
        solucionCrout(i, :) = c(i, :)
        do k = i + 1, orden
            solucionCrout(i, :) = solucionCrout(i, :) - aux(i, k) * solucionCrout(k, :)
        end do
    end do
 end function

function thomas(u_o, d_o, l_o, term_ind)
    real(8), dimension(:,:), intent(in) :: term_ind
    real(8), dimension(:), intent(in) :: u_o, d_o, l_o
    real(8), dimension(size(term_ind, dim=1), size(term_ind, dim=2)) :: thomas
    real(8), dimension(size(term_ind, dim=1)) :: u, d, l
    integer(4) filas, i

    filas = size(term_ind, DIM=1)
    u = u_o
    d = d_o
    l = l_o
    thomas = term_ind
    do i = 1, filas - 1
        u(i) = u(i) / d(i)
        thomas(i, :) = thomas(i, :) / d(i)
        d(i) = 1.0
        d(i + 1) = d(i + 1) - l(i + 1) * u(i)
        thomas(i + 1, :) = thomas(i + 1, :) - l(i + 1) * thomas(i, :)
        l(i + 1) = 0.0
    end do

    thomas(filas, :) = thomas(filas, :) / d(filas)
    do i = filas - 1, 1, -1
        thomas(i, :) = thomas(i, :) - u(i) * thomas(i + 1, :) / d(i)
    end do
end function thomas

function refinamientoIter(matriz, term_ind, tol, metodo, norma)
    real(8), dimension(:, :), intent(in) :: matriz, term_ind
    real(8), intent(in) :: tol
    procedure(solucionDirecto) :: metodo ! Tener en cuenta que solo necesita el valor de las incognitas
    procedure(mNorma) :: norma 
    real(8), dimension(size(term_ind, dim=1), size(term_ind, dim=2)) :: refinamientoIter, res
    real(8) error

    refinamientoIter = metodo(matriz, term_ind)
    res = residuo(matriz, refinamientoIter, term_ind)
    error = norma(res)
    do while (error > tol)
        refinamientoIter = refinamientoIter - metodo(matriz, res)
        res = residuo(matriz, refinamientoIter, term_ind)
        error = norma(res)
    end do
end function refinamientoIter

!!!!!!!!!!!!!!!!!!!!!!!!!! Metodos indirectos !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function jacobi(matriz, term_ind, xini, tol)
    real(8), dimension(:, :), intent(in) :: matriz, term_ind, xini
    real(8), intent(in) :: tol
    real(8), dimension(size(term_ind, dim = 1), size(term_ind, dim= 2 )) :: jacobi, xant
    real(8) e1, e2
    integer(4) i, j, n

    jacobi = xini
    n = size(jacobi, dim = 1)
    e1 = tol + 1
    e2 = tol + 1
    do while((e1 > tol) .and. (e2 > tol))
        xant = jacobi
        do i = 1, n
            jacobi(i, :) = term_ind(i, :)
            do j = 1, i - 1
                jacobi(i, :) = jacobi(i, :) - matriz(i, j) * xant(j, :)
            end do
            do j = i + 1, n
                jacobi(i, :) = jacobi(i, :) - matriz(i, j) * xant(j, :)
            end do
            jacobi(i, :) = jacobi(i, :) / matriz(i, i)
        end do
        e1 = maxval(abs(jacobi-xant))
        e2 = mNormaM(residuo(matriz, jacobi, term_ind))
    end do
end function jacobi

function gaussSeidel(matriz, term_ind, xini, tol)
    real(8), dimension(:, :), intent(in) :: matriz, term_ind, xini
    real(8), intent(in) :: tol
    real(8), dimension(size(term_ind, dim=1), size(term_ind, dim=2)) :: gaussSeidel, xant
    real(8) e1, e2
    integer(4) i, j, n, cont

    gaussSeidel = xini
    n = size(gaussSeidel, dim=1)
    e1 = tol + 1
    e2 = tol + 1
    cont = 0
    do while((e1 > tol) .and. (e2 > tol))
        xant = gaussSeidel
        do i = 1, n
            gaussSeidel(i, :) = term_ind(i, :)
            do j = 1, i - 1
                gaussSeidel(i, :) = gaussSeidel(i, :) - matriz(i, j) * gaussSeidel(j, :)
            end do
            do j = i + 1, n
                gaussSeidel(i, :) = gaussSeidel(i, :) - matriz(i, j) * gaussSeidel(j, :)
            end do
            gaussSeidel(i, :) = gaussSeidel(i, :) / matriz(i, i)
        end do
        e1 = errorRelativo(gaussSeidel, xant, mNormaM)
        e2 = mNormaM(residuo(matriz, gaussSeidel, term_ind))
        cont = cont + 1
    end do
    write(*, *) "Itereaciones: ", cont
end function gaussSeidel

function relajacion(matriz, term_ind, xini, tol)
    real(8), dimension(:, :), intent(in) :: matriz, term_ind, xini
    real(8), intent(in) :: tol
    real(8), dimension(size(term_ind, dim=1), size(term_ind, dim=2)) :: relajacion, xant, r
    real(8), parameter :: PI = 3.141592653589793238462643383279502884197169399375
    real(8) e1, e2, omega
    integer(4) i, j, n, cont, m

    relajacion = xini
    n = size(matriz, dim=1)
    m = size(matriz, dim=2)
    omega = (4./(2.+sqrt(4.-(cos(PI/(m-1.))+cos(PI/(n-1.)))**2)))
    write(*,*) 'Omega = ',omega
    e1 = tol + 1
    e2 = tol + 1
    cont = 0
    do while((e1 > tol) .and. (e2 > tol))
        xant = relajacion
        do i = 1, n
            relajacion(i, :) = term_ind(i, :)
            do j = 1, i - 1
                relajacion(i, :) = relajacion(i, :) - matriz(i, j) * relajacion(j, :)
            end do
            do j = i + 1, n
                relajacion(i, :) = relajacion(i, :) - matriz(i, j) * relajacion(j, :)
            end do
            relajacion(i, :) = relajacion(i, :) / matriz(i, i)
        end do
        r = relajacion - xant
        relajacion = xant + omega * r
        !relajacion = relajacion * omega + (1 - omega) * xant
        e1 = errorRelativo(relajacion, xant, mNormaM)
        e2 = mNormaM(residuo(matriz, relajacion, term_ind))
        cont = cont + 1
    end do
    write(*, *) "Itereaciones: ", cont
end function relajacion

function identidad(orden)
    integer(4), intent(in) :: orden
    integer(4) i
    real(8) identidad(orden, orden)

    identidad = 0.
    do i = 1, orden
        identidad(i, i) = 1.
    end do
end function identidad

function matrizInversa(matriz) ! Llamar solo con matrices cuadradas
    intent(in) :: matriz
    integer(4) orden
    real(8) matriz(:, :), matrizInversa(size(matriz, dim=1), size(matriz, dim=2))

    orden = size(matriz, dim=1)
    matrizInversa = solucionGaussJordan(matriz, identidad(orden))
end function

function residuo(mat, sol, term_ind)
    real(8), intent(in) :: mat(:, :), sol(:, :), term_ind(:, :)
    real(8) residuo(size(mat, dim=1), size(mat, dim=2))
    
    residuo = matmul(mat, sol) - term_ind
end function residuo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Medidas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function vNormaM(vector)
    intent(in) :: vector
    real(8) vector(:), vNormaM

    vNormaM = maxval(abs(vector))
end function vNormaM

function vNormaE(vector)
    intent(in) :: vector
    real(8) vector(:), vNormaE

    vNormaE = sqrt(sum(vector**2))
end function vNormaE

function mNormaM(matriz)
    intent(in) :: matriz
    real(8) matriz(:, :), mNormaM

    mNormaM = maxval(sum(abs(matriz), dim=2))
end function mNormaM

function mNormaL(matriz)
    intent(in) :: matriz
    real(8) matriz(:, :), mNormaL

    mNormaL = maxval(sum(abs(matriz), dim=1))
end function mNormaL

function mNormaF(matriz)
    intent(in) :: matriz
    real(8) matriz(:, :), mNormaF

    mNormaF = sqrt(sum(matriz**2))
end function mNormaF

function condicion(matriz, norma)
    intent(in) :: matriz
    procedure(mNorma) :: norma
    real(8) matriz(:, :), condicion

    condicion = norma(matriz) * norma(matrizInversa(matriz))
end function condicion

function errorAbsoluto(A, Aper, norma)
    real(8), dimension(:, :), intent(in) :: A, Aper
    procedure(mNorma) :: norma
    real(8) errorAbsoluto

    errorAbsoluto = norma(A - Aper)
end function errorAbsoluto

function errorRelativo(A, Aper, norma)
    real(8), dimension(:, :), intent(in) :: A, Aper
    procedure(mNorma) :: norma
    real(8) errorRelativo

    errorRelativo = errorAbsoluto(A, Aper, norma) / norma(A)
end function errorRelativo

function cotaErrorRelativo(A, Aper, b, bper, norma)
    real(8), dimension(:, :), intent(in) :: A, Aper, b, bper
    procedure(mNorma) :: norma
    real(8) cotaErrorRelativo, cond, erA, erb

    cond = condicion(A, norma)
    erA = errorRelativo(A, Aper, norma)
    erb = errorRelativo(b, bper, norma)
    cotaErrorRelativo = (cond / (1 - cond*erA)) * (erA + erb)
end function cotaErrorRelativo

end module SELs
