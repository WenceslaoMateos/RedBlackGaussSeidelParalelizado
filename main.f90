program main
    use RBGS
    use arreglos
    use SELs
    implicit none
    
    integer(4), parameter :: orden=4
    real(8), dimension(1:orden) :: d, term_ind, xini, res, ld, rd, res1
    real(8) tol

    !Elementos de prueba para Gauss-Seidel
    xini = [1., 2., 3., 4.]
    tol = 1e-5

    !Matriz y términos independientes para utilizar (ejemplo de diapositivas Crank-Nicholson)
    !term_ind = [2., 4., 4., 22.]
    !d = 4.
    !ld = [0., -1., -1., -1.]
    !rd = [-1., -1., -1., 0.]

    term_ind = [40.8, 0.8, 0.8, 200.8]
    d = 2.04
    ld = [0., -1., -1., -1.]
    rd = [-1., -1., -1., 0.]
    !res = gaussSeidel1D(d, ld, rd, term_ind, xini, tol)
    !resultados correctos [ 65.9698, 93.7784, 124.5382, 159.4795]

    !Resultados de Gauss-Seidel
    res = RBGS1D(d, ld, rd, term_ind, xini, tol)
    res1 = gaussSeidel1D(d, ld, rd, term_ind, xini, tol)

    if (this_image() == 1) then 
        write(*, *) 'Resultado de Red-Black Gauss-Seidel'
        call mostrarVector(res)
        write(*, *) 'Resultado de Gauss-Seidel Posta'
        call mostrarVector(res1)
    end if
    
end program main