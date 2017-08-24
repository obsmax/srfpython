        subroutine setmod(nmmodl,mmax)
c-----
c       COMMENTS
c       07 MAY 2002 - eliminated a write(2 to an unopened file
c-----
c       interactively set up model file
c-----
        parameter (LIN=5, LOT=6, LER=0)
        common/param/qaqb,itype,dlam,invdep
        character nmmodl*(*)

        common/modtit/title
        character title*80
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/depref/refdep

        write(LOT,*)
     1 ' Interactively setting up initial model file:', nmmodl
c-----
c       get model format
c-----
        WRITE(LOT,*)'Is model flat (0) or spherical (1)'
        read(LIN,*)iflsph
        WRITE(LOT,*)'Enter descriptive title for this model'
        read(LIN,'(a)')title

        write(LOT,*)' Enter d,a,b,rho,qa,qb'
        write(LOT,*)'  d=0.0 or EOF  indicates halfspace ',
     1     'and end of input'
c-----
c       get model data
c-----
        mmax = 0
 1000   continue
            read(LIN,*,err=1000,end=1001)v1,v2,v3,v4,v5,v6
            write(LOT,*)v1,v2,v3,v4,v5,v6
c-----
c           safety check
c-----
            if(v2.lt.v3)then
                WRITE(LOT,*)' Error: P Velocity < S velocity:',
     1              v2, '<', v3
                WRITE(LOT,*)'        Reenter this layer'
                go to 1000
            endif
            if(v4 .le. 0.0)then
                WRITE(LOT,*)'Error:  Density <= 0:',v4
                WRITE(LOT,*)'        Reenter this layer'
                go to 1000
            endif
            if(v5 .lt. 0.0 .or. v6.lt.0.0)then
                WRITE(LOT,*)'Error: qa or qb not >= 0',v5,v5
                WRITE(LOT,*)'        Reenter this layer'
                go to 1000
            endif
            mmax = mmax + 1
            d(mmax)      = v1
            a(mmax)      = v2
            b(mmax)      = v3
            rho(mmax)    = v4
            if(v5.gt.1.0)v5 = 1.0/v5
            if(v6.gt.1.0)v6 = 1.0/v6
            qa(mmax)     = v5
            qb(mmax)     = v6
            etap(mmax)   = 0.0
            etas(mmax)   = 0.0
            frefp(mmax)  = 1.0
            frefs(mmax)  = 1.0
            
            if(v1.eq.0.0)goto 1001
            goto 1000
 1001   continue
        d(mmax) = 0.0
        
c-----
c       lun I*4 - logical unit for writing model file. This
c                 unit is released after the use of this routine
c       nmmodl  C*(*)   - model name
c       mmax    I*4 - number of layers in the model, last layer is
c                    halfspace
c       title   C*(*)   - title of the model file
c       iunit   I*4 - 0 Kilometer, Gram, Sec
c       iiso    I*4 - 0 isotropic 
c                 1 transversely anisotropic 
c                 2 general anisotropic 
c       iflsph  I*4 - 0 flat earth model
c                 1 spherical earth model
c       idimen  I*4 - 1 1-D
c               - 2 2-D
c               - 3 3-D
c       icnvel  I*4 - 0 constant velocity
c                 1 variable velocity
c------
        iunit = 0
        iiso = 0
        idimen = 1
        icnvel = 0
        lt = lgstr(title)
        call putmod(2,nmmodl,mmax,title(1:lt),iunit,iiso,iflsph,
     1      idimen,icnvel,.true.)
        return
        end
