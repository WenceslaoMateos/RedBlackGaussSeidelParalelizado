real function normaVectorResiduoConVector(ld, d, rd, vectorResultado, b, n)
    integer n
    real ld(n), d(n), rd(n), vectorResultado(n), b(n)
    integer i, j
    real vectorR(n)
    
    vectorR = 0
  
    do i=1, n
      vectorR(i) = ld(i)*vectorResultado(i) + d(i)*vectorResultado(i) + rd(i)*vectorResultado(i) - b(i)
    enddo
    
    normaVectorResiduoConVector = normaVector(vectorR, n)
    
  endfunction

  real function normaVector(vector, n)
    integer n
    real vector(n)
    
    normaVector = maxval(abs(vector))
  endfunction
