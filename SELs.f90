module SELs

    implicit none

contains

function thomas(u_o, d_o, l_o, term_ind)
    real(8), dimension(:), intent(in) :: u_o, d_o, l_o, term_ind
    real(8), dimension(size(term_ind, dim=1)) :: thomas, u, d, l
    integer(4) filas, i

    filas = size(term_ind, DIM=1)
    u = u_o
    d = d_o
    l = l_o
    thomas = term_ind
    do i = 1, filas - 1
        u(i) = u(i) / d(i)
        thomas(i) = thomas(i) / d(i)
        d(i) = 1.0
        d(i + 1) = d(i + 1) - l(i + 1) * u(i)
        thomas(i + 1) = thomas(i + 1) - l(i + 1) * thomas(i)
        l(i + 1) = 0.0
    end do

    thomas(filas) = thomas(filas) / d(filas)
    do i = filas - 1, 1, -1
        thomas(i) = thomas(i) - u(i) * thomas(i + 1) / d(i)
    end do
end function thomas

!!!!!!!!!!!!!!!!!!!!!!!!!! Metodos indirectos !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function gaussSeidel1D(d, ld, rd, term_ind, xini, tol)
    real(8), dimension(:), intent(in) :: d, ld, rd, term_ind, xini
    real(8), intent(in) :: tol
    real(8), dimension(size(term_ind)) :: gaussSeidel1D, xant
    real(8) e1
    integer(4) i, cont, orden

    gaussSeidel1D = xini
    orden = size(xini)
    e1 = tol + 1
    cont = 0
    do while (e1 > tol)
        xant = gaussSeidel1D
        gaussSeidel1D(1) = (term_ind(1) - rd(1) * gaussSeidel1D(2)) / d(1)
        do i = 2, orden - 1
            gaussSeidel1D(i) = (term_ind(i) - ld(i) * gaussSeidel1D(i-1) - rd(i) * gaussSeidel1D(i+1)) / d(i)
        end do
        gaussSeidel1D(orden) = (term_ind(orden) - ld(orden) * gaussSeidel1D(orden-1)) / d(orden)
        e1 = errorAbsolutoV(gaussSeidel1D, xant)
        cont = cont + 1
    end do
    !write(*, *) "Iteraciones: ", cont
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Medidas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function vNormaM(vector)
    intent(in) :: vector
    real(8) vector(:), vNormaM

    vNormaM = maxval(abs(vector))
end function vNormaM

function errorAbsolutoV(V, Vper)
    real(8), dimension(:), intent(in) :: V, Vper
    real(8) errorAbsolutoV

    errorAbsolutoV = vNormaM(V - Vper)
end function errorAbsolutoV

! Funcion de refinamiento para el metodo Thomas
function refinamiento_thomas(ld_local, d_local, rd_local, term_local, tol)
    real(8), dimension(:), intent(in) :: ld_local, d_local, rd_local, term_local
    real(8), dimension(size(term_local, dim=1)) :: refinamiento_thomas, res
    real(8), intent(in) :: tol
    real(8) error
    integer(4) n

    n = ubound(term_local, 1)
    refinamiento_thomas = thomas(ld_local, d_local, rd_local, term_local)
    res = vectorResiduo(ld_local, d_local, rd_local, refinamiento_thomas, term_local, n)
    error = vNormaM(res)
    do while (error >= tol)
        refinamiento_thomas = refinamiento_thomas - thomas(ld_local, d_local, rd_local, res)
        res = vectorResiduo(ld_local, d_local, rd_local, refinamiento_thomas, term_local, n)
        error = vNormaM(res)
    end do
end function refinamiento_thomas

! Calculo del vector residuo para Thomas
function vectorResiduo(ld, d, rd, vectorResultado, b, n)
    integer(4), intent(in) :: n
    real(8), intent(in) :: ld(n), d(n), rd(n), vectorResultado(n), b(n)
    real(8) vectorResiduo(n)
    integer(4) i
    
    vectorResiduo = 0.
    vectorResiduo(1) = d(1)*vectorResultado(1) + rd(2)*vectorResultado(2) - b(1)
    do i = 2, n - 1
        vectorResiduo(i) = ld(i)*vectorResultado(i-1) + d(i)*vectorResultado(i) + rd(i)*vectorResultado(i+1) - b(i)
    enddo
    vectorResiduo(n) = ld(n-1)*vectorResultado(n-1) + d(n)*vectorResultado(n) - b(n)
end function

end module SELs
