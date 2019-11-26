program graficos
    implicit none
    
    call plot('benchmarkWen.dat', 'Intel Core i7-6500U 2.50GHz')
    call plot('benchmarkBraulio.dat', 'Intel Celeron N2840 2.16GHz')
    
    contains
    subroutine plot(archivo, titulo)
        character (LEN=*), intent(in) :: archivo, titulo

        open(unit=2, file="temporal.p", access='SEQUENTIAL', status='REPLACE')
        write(2, *) "set autoscale"
        write(2, *) "unset log"
        write(2, *) "unset label"
        write(2, *) "set grid"
        write(2, *) "set logscale y 10"
        write(2, *) "set logscale x 10"
        write(2, *) "set title '", titulo,"'"
        write(2, *) "set xlabel 'Orden'"
        write(2, *) "set ylabel 'Tiempo'"
        write(2, *) "plot '", archivo, "' using 1:2 title 'Gauss-Seidel' with lines,\"
        write(2, *) " '", archivo, "' using 1:3 title 'Thomas' with lines,\"
        write(2, *) " '", archivo, "' using 1:4 title 'RBGS' with lines,\"
        write(2, *) " '", archivo, "' using 1:5 title 'RBGS + Burocracia' with lines"
        call system('gnuplot -persist temporal.p')
        close(2, STATUS='DELETE')
    end subroutine plot
end program graficos