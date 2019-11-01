program main
    use RBGS
    use arreglos
    use SELs
    implicit none
    
    integer(4), parameter :: orden = 4
    real(8), dimension(orden, orden) :: matriz
    real(8), dimension(orden, 1) :: term_indAUX, xini
    real(8), dimension(orden) :: term_ind
    real(8) tol

    !Elementos de prueba para Gauss-Seidel
    xini(:, 1) = [0., 0., 0., 0.]
    tol = 1e-5

    !Matriz y términos independientes para utilizar
    matriz(1, :) = [0., -3., 2., 6.] 
    matriz(2, :) = [2., 3., 2., -1.] 
    matriz(3, :) = [-3., -1., 3., 1.] 
    matriz(4, :) = [1., 2., -3., -1.] 
    term_ind = [-8., -8., 0., 0.]

    !Resultados de Gauss-Seidel
    term_indAUX = gaussSeidel(matriz, term_indAUX, xini, tol)
    write(*, *) 'Resultado de Gauss-Seidel'
    call mostrarVector(term_indAUX(:, 1))

    
end program main