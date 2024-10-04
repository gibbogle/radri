!-----------------------------------------------------------------------------------------
! Main program
!-----------------------------------------------------------------------------------------
PROGRAM radri_main
use radri_mod
use global
implicit none
integer :: res
real(8) :: summarydata(100)
character*(128) :: infile, outfile, runfile
integer :: status, nlen, cnt, i, inbuflen, outbuflen
integer :: jstep,  irun     
character*(128) :: b, c, progname
real(8) :: t1, t2
integer :: idist, ndist = 40
real(8) :: PE, dist(40), ddist = 50   

real(8) :: progress(30)
integer :: phase(30), Nph
logical :: use_single

call get_command (b, nlen, status)
if (status .ne. 0) then
    write (*,*) 'get_command failed with status = ', status
    stop
end if
call get_command_argument (0, c, nlen, status)
if (status .ne. 0) then
    write (*,*) 'Getting command name failed with status = ', status
    stop
end if
progname = c(1:nlen)
cnt = command_argument_count ()
write (*,*) 'number of command arguments = ', cnt
if (cnt < 1) then
    write(*,*) 'Use: ',trim(progname),' input_file'
    stop
endif

do i = 1, cnt
    call get_command_argument (i, c, nlen, status)
    if (status .ne. 0) then
        write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        stop
    end if
    if (i == 1) then
        infile = c(1:nlen)
        write(*,*) 'Input file: ',infile
    endif
    use_PEST = .false.
end do
if (cnt == 2) then
    call get_command_argument (2, c, nlen, status)
    PEST_outputfile = c(1:nlen)     ! PEST is the parameter estimation program
    use_PEST = .true.
elseif (cnt /= 1) then
	write(*,*) 'Error: wrong number of arguments'
	stop
endif

if (use_PEST) then
    outfile = PEST_outputfile
else
    outfile = 'radri_main.out' 
endif

! Synchronisation of cell IR
use_synchronise = .false.
use_single = .false.    ! to simulate a cell (or cells) at specified synch_phase and synch_progress
synch_phase = S_phase
synch_fraction = 0.5
nph = 1
if (use_synchronise) then
    if (use_single) then
        nph = 1
        progress(1) = synch_fraction
    else
        nph = 10
        do i = 1,nph
            progress(i) = (i-1)*1.0/nph
        enddo
    endif
endif

do irun = 1,nph
    if (use_synchronise) then
        synch_fraction = progress(irun)
        write(*,*)
    	write(*,'(a,2i4,f6.3)') 'radri_main: irun, synch_phase, synch_fraction: ',irun,synch_phase,synch_fraction
    endif
	inbuflen = len(infile)
	outbuflen = len(outfile)
	write(*,*) 'call execute'
    res = irun
	call execute(infile,inbuflen,outfile,outbuflen,res)
	if (res /= 0) stop
	t1 = wtime()
	write(*,*) 'did execute: nsteps: ',nsteps
	do jstep = 1,Nsteps+1
		call simulate_step(res)
		if (res == 1) then
		    write(*,*) 'Success'
		    exit
		elseif (res == 2) then
			write(*,*) 'All cells have died'
			stop
		elseif (res /= 0) then
			write(*,*) 'Error exit: ',res
			stop
		endif
	enddo
    write(*,*) 'res: ',res
    if (res == 0) write(*,*) 'Exceeded nsteps, not all cells reached mitosis, increase ndays'
	call terminate_run(res)
	t2 = wtime()
	write(*,*) 'time: ',t2-t1
enddo
end

