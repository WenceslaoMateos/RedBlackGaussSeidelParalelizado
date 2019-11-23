module moduloArchivos
    implicit none
  
    contains
  
    subroutine escribeArchivoVector(vector, n, nombreArchivo)
        real vector(n)
        integer n
        character(len = 20) :: nombreArchivo
        
        open (1, file = nombreArchivo, status='replace')
        write(1,*) n
        write(1,*) vector(:)
        close (1, status='keep')
    end subroutine

    subroutine escribeArchivoMatriz(matriz, n, m, nombreArchivo)
        integer n, m, i
        real matriz(n, m)
        character(len = 20) :: nombreArchivo
        
        open (1, file = nombreArchivo, status='replace')
        write(1,*) n, m
        write(1,*) matriz(1,:)
        do i=2, n-1
            write(1,*) matriz(i,:)
!~             Puedo evitar la redundancia si solo escribo los bordes
!~             write(1,*) matriz(i,1)
!~             write(1,*) matriz(i,m)
        end do
        write(1,*) matriz(n,:)
        close (1, status='keep')
    end subroutine

    subroutine leeArchivoVector(vector, n, nombreArchivo)
        integer n
        real, allocatable :: vector(:)
        character(len = 20) :: nombreArchivo
        logical :: existeArchivo
        
        if (.not. allocated(vector)) then
            inquire(file = nombreArchivo, exist=existeArchivo)
            if (existeArchivo) then
                open (1, file = nombreArchivo)
                read (1, *) n
                allocate(vector(n))
                read (1, *) vector(:)
                close (1, status='keep')
            else
                write(*,*) "No existe el archivo: ", nombreArchivo, "Por favor, crearlo."
                error stop "No existe el archivo."
            end if
        else
            error stop "El vector ya esta inicializado."
        end if
    end subroutine

    subroutine leeArchivoMatriz(matriz, n, m, nombreArchivo)
        integer n, m, i
        real, allocatable :: matriz(:,:)
        character(len = 20) :: nombreArchivo
        logical :: existeArchivo

        if (.not. allocated(matriz)) then
            inquire(file = nombreArchivo, exist=existeArchivo)
            if (existeArchivo) then
                open (1, file = nombreArchivo)
                read (1, *) n, m
                allocate(matriz(n, m))
                do i=1, n
                    read (1, *) matriz(i, :)
                end do
                close (1, status='keep')
            else
                write(*,*) "No existe el archivo: ", nombreArchivo, "Por favor, crearlo."
                error stop "No existe el archivo."
            end if
        else
            error stop "La matriz ya esta inicializada."
        end if
    end subroutine
endmodule
