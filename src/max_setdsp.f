        subroutine setdsp(nmdisp)
c-----
c       Changes
c       21 JUL 2002 - ensure proper format for GAMMA
c-----
c       interactively set up dispersion file
c-----
        parameter (LIN=5, LOT=6, LER=0)
        character nmdisp*(*)
        character*1 lorr(2), porg(3)
        character ostr*12
        data lorr/'L','R'/,porg/'C','U','G'/

        write(LOT,*)
     1 ' Interactively setting up dispersion file:',nmdisp
c-----
c       open file
c-----
        open(2,file=nmdisp,status='unknown',form='formatted',
     1           access='sequential')
        rewind 2
c-----
c       NEW UNITS ARE ALWAYS km AND DATA IS ALWAYS PERIOD
c-----
c       prompt for dispersion information
c-----
        write(LOT,*)
     1 ' Enter ilvry,iporg,imode,per,val,dval'
        write(LOT,*)
     1  ' ilvry=1(Love)'
        write(LOT,*)
     1  '      =2(Rayleigh)'
        write(LOT,*)
     1  ' iporg=1 (phase velocity km/s)'
        write(LOT,*)
     1  '      =2 (group velocity km/s)'
        write(LOT,*)
     1  '      =3 (gamma 1/km)'
        write(LOT,*)
     1  ' imode (mode number) e.g., 0=fundamental, 1=first'
        write(LOT,*)
     1  ' per=the period '
        write(LOT,*)
     1 '    val=dispersion value, velocity or gamma'
        write(LOT,*)
     1 '   dval=error in dispersion value'
        write(LOT,*)
     1  '       (Enter 1.0 if stderr from residuals)'
        write(LOT,*)
     1 '   NOTE: Enter all zeros or negative to terminate input'
 1000 continue
        read(LIN,*,end=1001,err=1001)ilvry,iporg,imode,per,val,dval
        if(ilvry.le.0 .or. imode.lt.0)goto 1001
            ostr(1:7) = 'SURF96 '
            ostr(8:8) = lorr(ilvry)
            ostr(9:9) = ' '
            ostr(10:10) = porg(iporg)
            ostr(11:12) = ' X'
            obserr = 0.0
            if(iporg.eq.1 .or. iporg.eq.2)then
            write(2,'(a12,1x,i3,1x,3f11.4)') ostr,imode,per,val,dval
            else if(iporg.eq.3)then
            write(2,'(a12,1x,i3,1x,f11.4,2e11.3)') 
     1         ostr,imode,per,val,dval
            endif
        goto 1000
 1001 continue
        close (2)
        return
        end
