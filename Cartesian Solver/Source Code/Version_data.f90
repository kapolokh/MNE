subroutine Version_data(out)
        implicit none
        integer, intent(in) :: out
        character(len=32) :: execution_date, execution_time
        call date_and_time(date=execution_date,time=execution_time)
        write(out,'(a)') "Transport Solver"
        write(out,'(a)') "Version 4.0"
        write(out,'(a)') "Kirill Polokhalo"
        write(out,'(a)') "Date:", execution_date
        write(out,'(a)') "Time:", execution_time
        write(out,'(a)') "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    end subroutine Version_data