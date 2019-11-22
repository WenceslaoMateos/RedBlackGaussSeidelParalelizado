program main
    use RBGS
    use arreglos
    use SELs
    implicit none
    
    real(8), dimension(:), codimension[:], allocatable :: d, term_ind, xini, res, ld, rd
    real(8), dimension(:), allocatable :: d_local, term_local, xini_local, ld_local, rd_local
    real(8) tol
    integer(4) tam_divisiones[*], im_act, im_tot, i, inicio, fin, orden, remanente

    im_act = this_image()
    im_tot = num_images()
    tol = 1e-5

    if (im_act == 1) then
        xini_local = [1., 2., 3., 4.]
        term_local = [40.8, 0.8, 0.8, 200.8]
        d_local = [2.04, 2.04, 2.04, 2.04]
        ld_local = [0., -1., -1., -1.]
        rd_local = [-1., -1., -1., 0.]
        orden = size(term_local)
        tam_divisiones = ceiling(real(orden) / real(im_tot))
        call co_broadcast(tam_divisiones, 1)
    end if

    allocate(d(tam_divisiones)[*], term_ind(tam_divisiones)[*], xini(tam_divisiones)[*], &
        res(tam_divisiones)[*], ld(tam_divisiones)[*], rd(tam_divisiones)[*])
    
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

        remanente = mod(orden, tam_divisiones)
        if (remanente /= 0) then
            inicio = orden - remanente + 1
            d(1:remanente)[im_tot] = d_local(inicio:)
            ld(1:remanente)[im_tot] = ld_local(inicio:)
            rd(1:remanente)[im_tot] = rd_local(inicio:)
            xini(1:remanente)[im_tot] = xini_local(inicio:)
            term_ind(1:remanente)[im_tot] = term_local(inicio:)
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
    end if

    sync all
    write(*, *) xini
    call RBGSlineal(res, d, ld, rd, term_ind, xini, tol)

    sync all
    if (this_image() == 1) then 
        write(*, *) 'Resultado de Red-Black Gauss-Seidel'
        do i = 1, im_tot
            call mostrarVector(res(:)[i])
        end do
    end if
    
    deallocate(res, d, ld, rd, term_ind, xini)
end program main