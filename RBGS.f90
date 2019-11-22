module RBGS
    use SELs
    implicit none
    
contains

    function RBGS1D(d, ld, rd, term_ind, xini, tol)
        implicit none
        real(8), dimension(:), intent(in) :: d, ld, rd, term_ind, xini
        real(8), intent(in) :: tol
        real(8), dimension(size(term_ind)) :: RBGS1D
        real(8), dimension(:), codimension[:], allocatable :: res, xant
        real(8) e1
        integer(4) i, cont, orden, inicio, fin, tam_divisiones, im_act, arranca, termina, cota
        
        !definimos los inicios y finales del coarray correspondiente al problema real
        orden = size(xini)
        tam_divisiones = ceiling(real(orden) / real(num_images()))
        inicio = tam_divisiones * this_image() - tam_divisiones + 1
        fin = this_image() * tam_divisiones
        allocate(res(tam_divisiones)[*])
        if (this_image() < num_images()) then
            res = xini(inicio:fin)
        else
            res(:ubound(xini,1)-inicio) = xini(inicio:)
        end if
        e1 = tol + 1
        cont = 0
        cota = orden - tam_divisiones * (num_images()-1)
        !cambiar la cota de iteraciones
        do while (cont < 36)
            xant = res

            !calculamos los rojos (impares)
            
            !calculo de los nodos del principio
            if (this_image() == 1) then
                !nodo primero de todos
                res(1) = (term_ind(1) - rd(1) * res(2)) / d(1)
                arranca = 3
            else
                !nodos al principio de cada imagen, que no son el primero
                if (mod(inicio, 2) == 1) then
                    !nodo impar, osea rojo al principio, y calculo con la imagen anterior
                    res(1) = (term_ind(inicio) - ld(inicio) * res(tam_divisiones)[this_image()-1] - rd(inicio) * res(2)) / d(inicio)
                    arranca = 3
                else
                    !nodo par al principio, asi que vamos al siguiente primero rojo
                    arranca = 2
                end if
            end if

            !Calculo de nodos del final
            if (this_image() == num_images()) then
                !calculo del nodo del final de todos
                if (mod(fin, 2) == 1) then
                    !final rojo
                    res(cota) = (term_ind(orden) - ld(orden) * res(cota-1)) / d(orden)
                end if
                !sea o no el final un nodo rojo, lo salteamos, por que si es, ya lo calculamos, y si no es, no lo usamos
                termina = cota-1
            else
                if (mod(fin, 2) == 1) then
                    !nodo al final es rojo
                    res(tam_divisiones) = (term_ind(fin) - ld(fin) * res(tam_divisiones-1) &
                        - rd(fin) * res(1)[this_image()+1]) / d(fin)
                end if
                !si el nodo al final no es rojo lo salteamos, y si es, ya esta calculado
                termina = tam_divisiones - 1
            end if

            !calculo de los nodos del medio
            do i = arranca, termina, 2
                res(i) = (term_ind(i+inicio) - ld(i+inicio) * res(i-1) - rd(i+inicio) * res(i+1)) / d(i+inicio)
            end do

            sync all

            !calculamos los negros (pares)

            !calculo de los nodos del principio (nunca va a ser negro el primero del vector original)
            if (mod(inicio, 2) == 0) then
                res(1) = (term_ind(inicio) - ld(inicio) * res(tam_divisiones)[this_image()-1] - rd(inicio) * res(2)) / d(inicio)
                arranca = 3
            else
                arranca = 2
            end if

            !calculo de los nodos del final del todo
            if (this_image() == num_images()) then
                !estoy en la ultima imagen
                if (mod(fin, 2) == 0) then
                    !final negro
                    res(cota) = (term_ind(orden) - ld(orden) * res(cota-1)) / d(orden)
                end if
                !si el final es negro, ya lo calcule, y si es rojo, no preciso calcularlo
                termina = cota-1
            else
                if (mod(fin, 2) == 0) then
                    !el ultimo es negro y lo calcculo
                    res(tam_divisiones) = (term_ind(fin) - ld(fin) * res(tam_divisiones-1) &
                        - rd(fin) * res(1)[this_image()+1]) / d(fin)
                end if
                !si es negro, ya lo calcule y lo evito, sino, no lo toco
                termina = tam_divisiones - 1
            end if

            !calculo de los nodos del medio
            do i = arranca, termina, 2
                res(i) = (term_ind(i+inicio) - ld(i+inicio) * res(i-1) - rd(i+inicio) * res(i+1)) / d(i+inicio)
            end do
            sync all

            e1 = errorAbsolutoV(res, xant, vNormaM)
            cont = cont + 1
        end do
        !write(*, *) "Iteraciones: ", cont
        sync all
        if (this_image() == 1) then
            do i = 1, num_images() - 1 
                RBGS1D((i-1)*tam_divisiones+1:i*tam_divisiones) = res(:)[i]
            end do
            RBGS1D(1+tam_divisiones*(num_images()-1):) = res(:ubound(RBGS1D,1)-1+tam_divisiones*(num_images()-1))[num_images()]
        end if
    end function RBGS1D

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

            !calculo de los nodos del medio
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
        
        write(*, *) "Iteraciones: ", cont
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