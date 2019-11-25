program main
    use RBGS
    use SELs
    use casosDinamicos
    implicit none
    
    real(8), dimension(:), codimension[:], allocatable :: d, term_ind, xini, res, ld, rd
    real(8), dimension(:), allocatable :: d_local, term_local, xini_local, ld_local, rd_local, res1
    real(8) tol[*], t_ini, t_fin, t_normal, t_concurrente, t_resultado, t_thomas, error, t_trans, t_concurrente_tot
    integer(4) tam_divisiones[*], im_act, im_tot, i, inicio, fin, orden[*], remanente, cant

    im_act = this_image() ! Mi imagen
    im_tot = num_images() ! Cantidad total de imagenes

    cant = 50    

    ! La imagen 1 se encarga de la carga de datos
    if (im_act == 1) then
		write (*, *) "Ingrese el orden:"
        read(*, *) orden
        write(*, *) "Ingrese la tolerancia:"
        read(*, *) tol
        ! Estos son los unicos datos que se deben tocar para modificar el programa
        xini_local = generaCaso(orden)
        term_local = generaCaso(orden)
        allocate(d_local(orden), ld_local(orden), rd_local(orden))
        d_local = 4.
        ld_local = -1.
        ld_local(1) = 0.
        rd_local = -1.
        rd_local(orden) = 0.
        
        tam_divisiones = ceiling(real(orden) / real(im_tot))
        ! Hay que propagar la cantidad de divisiones a todas las imagenes
        ! Ademas se propaga el orden y la tolerancia
        do i = 2, im_tot
            tam_divisiones[i] = tam_divisiones
            orden[i] = orden
            tol[i]= tol
        end do
        t_normal = 0
        do i = 1, cant
            call cpu_time(t_ini)
            res1 = gaussSeidel1D(d_local, ld_local, rd_local, term_local, xini_local, tol)
            call cpu_time(t_fin)
            t_normal = t_normal + t_fin-t_ini
        end do
        t_normal = t_normal/cant
        ! write(*, *) 'Resultado de Gauss-Seidel'
        ! call mostrarVector(res1)
        write(*, '(A,F10.7)') 'Tiempo Gauss-Seidel = ', t_normal

        t_thomas = 0
        do i = 1, cant
            call cpu_time(t_ini)
            res1 = refinamiento_thomas(ld_local, d_local, rd_local, term_local, tol)
            call cpu_time(t_fin)
            t_thomas = t_thomas + t_fin-t_ini
        end do
        t_thomas = t_thomas/cant
        ! write(*, *) 'Resultado de Thomas'
        ! call mostrarVector(res1)
        write(*, '(A,F10.7)') 'Tiempo Thomas = ', t_thomas
        write(*, *)

    end if
    
    ! No funciona si no todas las imagenes tienen el mismo tam_divisiones
    sync all
    allocate(d(tam_divisiones)[*], term_ind(tam_divisiones)[*], xini(tam_divisiones)[*], &
        res(tam_divisiones)[*], ld(tam_divisiones)[*], rd(tam_divisiones)[*])
    sync all
    ! La imagen 1 se encarga de distribuir los datos
    if (im_act == 1) then
        call cpu_time(t_ini)
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
        
        call cpu_time(t_fin)
        t_trans = t_fin - t_ini
        write(*, '(A,F12.7)') 'Tiempo de transmision: ',t_trans

        deallocate(d_local, term_local, xini_local, ld_local, rd_local)
    end if

    
    t_concurrente = 0
    do i = 1, cant
        call cpu_time(t_ini)
        sync all
        call RBGSlineal(res, d, ld, rd, term_ind, xini, tol)
        sync all
        call cpu_time(t_fin)
        t_concurrente = t_concurrente + t_fin-t_ini
    end do

    if (im_act == 1) then 
        ! write(*, *) 'Resultado de Red-Black Gauss-Seidel'
        ! do i = 1, im_tot
        !     call mostrarVector(res(:)[i])
        ! end do
        t_concurrente = t_concurrente /cant
        t_concurrente_tot = t_concurrente + t_trans
        write(*, '(A,F10.7)') 'Tiempo RBGS = ', t_concurrente
        write(*, '(A,F10.7)') 'Tiempo RBGS + Burocracia = ', t_concurrente_tot
        write(*, *)
        write(*, '(A,F12.7,A)') 'Optimización RBGS vs GS = ', (t_normal/t_concurrente),' veces mejor'
        write(*, '(A,F12.7,A)') 'Optimización RBGS vs Thomas = ', (t_thomas/t_concurrente),' veces mejor'
        write(*, *)
        write(*, '(A,F12.7,A)') 'Optimización RBGS + Burocracia vs GS = ', (t_normal/t_concurrente_tot),' veces mejor'
        write(*, '(A,F12.7,A)') 'Optimización RBGS + Burocracia vs Thomas = ', (t_thomas/t_concurrente_tot),' veces mejor'
    end if
    
    deallocate(res, d, ld, rd, term_ind, xini)
end program main
