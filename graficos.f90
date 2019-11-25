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
        write(2, *) "set xlabel 'Tiempo'"
        write(2, *) "set ylabel 'Orden'"
        write(2, *) "plot '", archivo, "' using 2:1 title 'Gauss-Seidel' with lines,\"
        write(2, *) " '", archivo, "' using 3:1 title 'Thomas' with lines,\"
        write(2, *) " '", archivo, "' using 4:1 title 'RBGS' with lines,\"
        write(2, *) " '", archivo, "' using 5:1 title 'RBGS + Burocracia' with lines"
        call system('gnuplot -persist temporal.p')
        close(2, STATUS='DELETE')
    end subroutine plot
end program graficos