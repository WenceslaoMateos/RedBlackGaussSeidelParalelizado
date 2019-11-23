module RBGS
    use SELs
    implicit none
    
contains

    subroutine RBGSlineal(res, d, ld, rd, term_ind, xini, tol)
        real(8), dimension(:), codimension[*], intent(in) :: d, ld, rd, term_ind, xini
        real(8), dimension(:), codimension[*], intent(inout) :: res
        real(8), intent(in) :: tol
        real(8), dimension(:), codimension[:], allocatable :: xant
        real(8) e
        real(8), codimension[:], allocatable :: error_local
        integer(4) i, cont, tam_divisiones, im_act, arranca_rojo, arranca_negro

        im_act = this_image()
        tam_divisiones = size(xini)
        res = xini
        allocate(xant(tam_divisiones)[*], error_local[*])
        e = tol + 1
        cont = 0
        do while (e >= tol)
            xant = res

            ! Calculo de rojos

            ! Calculo de los nodos del principio
            if (this_image() == 1) then
                !nodo primero de todos
                res(1) = (term_ind(1) - rd(1) * res(2)) / d(1)
                arranca_rojo = 3
                arranca_negro = 2
            else
                ! El array local empieza en rojo si cada array tiene una cantidad par de elementos
                ! o estamos en una imagen impar
                if (mod(tam_divisiones, 2) == 0 .or. mod(im_act, 2) == 1) then
                    !nodo impar, osea rojo al principio, y calculo con la imagen anterior
                    res(1) = (term_ind(1) - ld(1) * res(tam_divisiones)[im_act-1] - rd(1) * res(2)) / d(1)
                    arranca_rojo = 3
                    arranca_negro = 2
                else
                    ! Nodo par al principio, asi que vamos al siguiente primero rojo
                    arranca_rojo = 2
                    arranca_negro = 3
                end if
            end if

            ! Calculo de los nodos del medio
            do i = arranca_rojo, tam_divisiones - 1, 2
                res(i) = (term_ind(i) - ld(i) * res(i-1) - rd(i) * res(i+1)) / d(i)
            end do

            ! El ultimo nodo es rojo
            if (i == tam_divisiones) then
                if (im_act == num_images()) then
                    res(tam_divisiones) = (term_ind(tam_divisiones) - ld(tam_divisiones) * res(tam_divisiones-1)) &
                        / d(tam_divisiones)
                else
                    res(tam_divisiones) = (term_ind(tam_divisiones) - ld(tam_divisiones) * res(tam_divisiones-1) &
                        - rd(tam_divisiones) * res(1)[im_act+1]) / d(tam_divisiones)
                end if
            end if

            sync all

            ! Calculo de negros 

            ! Negro en posicion 1 (nunca puede darse en la imagen 1)
            if (arranca_negro == 3) then
                res(1) = (term_ind(1) - ld(1) * res(tam_divisiones)[im_act-1] - rd(1) * res(2)) / d(1)
            end if

            ! Calculo de los nodos del medio
            do i = arranca_negro, tam_divisiones - 1, 2
                res(i) = (term_ind(i) - ld(i) * res(i-1) - rd(i) * res(i+1)) / d(i)
            end do

            ! El ultimo nodo es negro
            if (i == tam_divisiones) then
                if (im_act == num_images()) then
                    res(tam_divisiones) = (term_ind(tam_divisiones) - ld(tam_divisiones) * res(tam_divisiones-1)) &
                        / d(tam_divisiones)
                else
                    res(tam_divisiones) = (term_ind(tam_divisiones) - ld(tam_divisiones) * res(tam_divisiones-1) &
                        - rd(tam_divisiones) * res(1)[im_act+1]) / d(tam_divisiones)
                end if
            end if

            error_local = errorAbsolutoV(res, xant, vNormaM)
            sync all
            e = maxCoarrayEscalar(error_local)
            cont = cont + 1
        end do
        
        ! if (im_act == 1) then
        !     write(*, *) "Iteraciones: ", cont
        ! end if
        deallocate(xant, error_local)
    end subroutine RBGSlineal

    function maxCoarrayEscalar(coarray)
        real(8), intent(in) :: coarray[*]
        real(8) maxCoarrayEscalar
        integer(4) i

        maxCoarrayEscalar = coarray[1]
        do i = 2, num_images()
            if (coarray[i] > maxCoarrayEscalar) then
                maxCoarrayEscalar = coarray[i]
            end if
        end do
    end function maxCoarrayEscalar

end module RBGS