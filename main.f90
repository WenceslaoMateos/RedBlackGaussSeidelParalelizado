program main
    use RBGS
    use arreglos
    use SELs
    use moduloArchivos
    implicit none
    
    real(8), dimension(:), codimension[:], allocatable :: d, term_ind, xini, res, ld, rd
    real(8), dimension(:), allocatable :: d_local, term_local, xini_local, ld_local, rd_local, res1
    real(8) tol, t_ini, t_fin, t_normal, t_concurrente, t_resultado
    integer(4) tam_divisiones[*], im_act, im_tot, i, inicio, fin, orden, remanente, n

    im_act = this_image() ! Mi imagen
    im_tot = num_images() ! Cantidad total de imagenes
    tol = 1e-5

    ! La imagen 1 se encarga de la carga de datos
    if (im_act == 1) then
        ! Estos son los unicos datos que se deben tocar para modificar el programa
        call leeArchivoVector(xini_local, n, 'Datos/xini.dat')
        call leeArchivoVector(term_local, n, 'Datos/term.dat')
        allocate(d_local(n), ld_local(n), rd_local(n))
        d_local = 5.
        ld_local = -1.
        ld_local(1) = 0.
        rd_local = -1.
        rd_local(n) = 0.
        
        !call leeArchivoVector(d_local, n, 'Datos/d.dat')
        !call leeArchivoVector(ld_local, n, 'Datos/ld.dat')
        !call leeArchivoVector(rd_local, n, 'Datos/rd.dat')

        orden = size(term_local)
        tam_divisiones = ceiling(real(orden) / real(im_tot))
        ! Hay que propagar la cantidad de divisiones a todas las imagenes
        do i = 2, im_tot
            tam_divisiones[i] = tam_divisiones
        end do

        call CPU_TIME(t_ini)
        res1 = gaussSeidel1D(d_local, ld_local, rd_local, term_local, xini_local, tol)
        call CPU_TIME(t_fin)
        !write(*, *) 'Resultado de Gauss-Seidel Posta'
        !call mostrarVector(res1)
        t_normal = t_fin-t_ini
        write(*, *) 'Tiempo = ', t_normal

    end if
    
    
    ! No funciona si no todas las imagenes tienen el mismo tam_divisiones
    sync all
    allocate(d(tam_divisiones)[*], term_ind(tam_divisiones)[*], xini(tam_divisiones)[*], &
        res(tam_divisiones)[*], ld(tam_divisiones)[*], rd(tam_divisiones)[*])
    sync all
    ! La imagen 1 se encarga de distribuir los datos
    if (im_act == 1) then
        inicio = 1
        fin = tam_divisiones
        do i = 1, im_tot - 1
            d(:)[i] = d_local(inicio:fin)
            ld(:)[i] = ld_local(inicio:fin)
            rd(:)[i] = rd_local(inicio:fin)
            xini(:)[i] = xini_local(inicio:fin)
            term_ind(:)[i] = term_local(inicio:fin)
            inicio = fin + 1
            fin = fin + tam_divisiones
        end do

        ! La última imagen caso puede quedar con remanente
        remanente = mod(orden, tam_divisiones)
        if (remanente /= 0) then
            inicio = orden - remanente + 1
            d(1:remanente)[im_tot] = d_local(inicio:)
            ld(1:remanente)[im_tot] = ld_local(inicio:)
            rd(1:remanente)[im_tot] = rd_local(inicio:)
            xini(1:remanente)[im_tot] = xini_local(inicio:)
            term_ind(1:remanente)[im_tot] = term_local(inicio:)

            ! Estos son los valores para que las componentes sobrantes no interfieran
            ! en la resolucion del metodo
            do i = remanente + 1, tam_divisiones
                d(i)[im_tot] = 1.
                ld(i)[im_tot] = 0.
                rd(i)[im_tot] = 0.
                xini(i)[im_tot] = 0.
                term_ind(i)[im_tot] = 0.
            end do
        else
            d(:)[i] = d_local(inicio:fin)
            ld(:)[i] = ld_local(inicio:fin)
            rd(:)[i] = rd_local(inicio:fin)
            xini(:)[i] = xini_local(inicio:fin)
            term_ind(:)[i] = term_local(inicio:fin)
        end if
        
        deallocate(d_local, term_local, xini_local, ld_local, rd_local)
        write(*, *)
    end if
    

    sync all
    call CPU_TIME(t_ini)
    call RBGSlineal(res, d, ld, rd, term_ind, xini, tol)
    sync all
    call CPU_TIME(t_fin)

    if (im_act == 1) then 
        !write(*, *) 'Resultado de Red-Black Gauss-Seidel'
        do i = 1, im_tot
            !call mostrarVector(res(:)[i])
        end do
        t_concurrente = t_fin-t_ini
        write(*, *) 'Tiempo = ', t_concurrente
        write(*, *) 'Optimización = ', t_normal/t_concurrente *100,'%'
    end if
    
    deallocate(res, d, ld, rd, term_ind, xini)
end program main