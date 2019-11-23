error = normaVectorResiduoConVector(matrizInicial, vectorSumatoriaFila, n, m)

real function normaVectorResiduoConVector(matrizInicial, vectorResultado, n, m)
    integer n, m
    real matrizInicial(n, m), vectorResultado(n)
    integer i, j
    real vectorR(n)
    
    vectorR = 0
    
    do i=1, n
      do j=1, n
        vectorR(i) = vectorR(i) + matrizInicial(i,j)*vectorResultado(j)
      enddo
      vectorR(i) = vectorR(i) - matrizInicial(i,m)
    enddo
    
    normaVectorResiduoConVector = normaVector(vectorR, n)
    
  endfunction

  real function normaVector(vector, n)
    integer n
    real vector(n)
    
    normaVector = maxval(abs(vector))
  endfunction
