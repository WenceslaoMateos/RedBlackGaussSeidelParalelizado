real function normaVectorResiduoConVector(ld, d, rd, vectorResultado, b, n)
    integer n
    real ld(n), d(n), rd(n), vectorResultado(n), b(n)
    integer i, j
    real vectorR(n)
    
    vectorR = 0
    
    vectorR(i) = d(1)*vectorResultado(1) + rd(2)*vectorResultado(2) - b(1)
    do i=2, n-1
      vectorR(i) = ld(i)*vectorResultado(i-1) + d(i)*vectorResultado(i) + rd(i)*vectorResultado(i+1) - b(i)
    enddo
    vectorR(i) = ld(n-1)*vectorResultado(n-1) + d(n)*vectorResultado(n) - b(n)
    
    normaVectorResiduoConVector = normaVector(vectorR, n)
    
  end function

  real function normaVector(vector, n)
    integer n
    real vector(n)
    
    normaVector = maxval(abs(vector))
  end function
