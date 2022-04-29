c     
c     Program ldbnds.f, using subroutines GLAN and LANDEM
c
c     Performs computations related to group sequential boundaries 
c     using the spending function approach.  Interactive input 
c     routine has four options which compute:
c     (1) bounds given analysis times and a spending function;
c     (2) drift parameters corresponding to given power and bounds;
c     (3) probabilities given times, bounds and drift parameters;
c     (4) a confidence interval given times, bounds and final
c         test statistic value.
c     
c     The essential computations are done in two subroutines:
c     GLAN and LANDEM.  GLAN computes probabilities from bounds,
c     and LANDEM computes bounds given a spending fuction.  The 
c     long input routine is designed to be run interactively, but 
c     can be run non-interactively.  
c     
c     TO DO:
c      bisect routine won't compute power correctly if drift < 0;
c      use structured progamming to simplify input routine code.
c     
c     References:
c      Armitage, McPherson & Rowe 1969 JRSSA
c      McPherson & Armitage 1971 JRSSA
c      Pocock 1977 Biometrika
c      Lan & DeMets 1983 Biometrika
c      Kim & DeMets 1987 Biometrics
c      Lan, Reboussin & DeMets 1994 Comm Stat A
c      Reboussin, DeMets, Kim & Lan 1992 UW Biostat TR 60
c      Reboussin, DeMets, Kim & Lan 1992 UW Biostat TR 95
c      http://www.biostat.wisc.edu/landemets/
c     
c     Partial variable list:
c      n is the number of interim analyses.
c      t is the vector of analysis times on (0,1] scale.
c      t2 is the second time scale for covariances.
c      alpha is the type I error level.
c      conf is the confidence level (or target power).
c      side indicates one or two sided bounds.
c      za is the standardized lower boundary.
c      zb is the standardized upper boundary.
c      drift(i) is the true mean of process
c      pstp(i) is the prob of reaching ith analysis and stopping.
c      qpos(i) is the prob of reaching ith and exceeding upper.
c      qneg(i) is the prob of reaching ith and exceeding lower.
c      pr is the total prob of rejecting, sum(1,n)(qpos+qneg)
c      pd(i) is the difference in exit probability betw i, i-1
c      pe(i) is the cumulative exit probabilty
c      power(i) is the prob of rejecting under drift(i)
c      limit(i) contains the confidence limits for drift
c      est is the expected stopping time (not implemented)
c     
      program ldbnds
c
      integer ndim
      parameter(ndim=25)
c     
      double precision t(ndim),t2(ndim),za(ndim),zb(ndim),drift(ndim)
      double precision pstp(ndim),qpos(ndim),qneg(ndim)
      double precision pd(ndim),pe(ndim)
      double precision power(ndim)
      double precision alpha,conf,side,pr,qrcum
      integer n
c     
      double precision znorm
      double precision limit(2),est
      integer nopt,ndrifs,iuse,iinf
      character char
c
 10   continue
c     
      call input(nopt,n,ndim,t,iinf,t2,alpha,side,iuse,za,zb,
     .     ndrifs,drift,conf,pd,pe)
c     
      if (nopt.eq.1) then
c        Output bounds from spending function.
         call outld(n,alpha,iuse,side,iinf,t,t2,pd,za,zb)
      else if (nopt.eq.2) then
c        Return drift for given bounds and power
c        Find drift resulting in power "conf"
c        Starting value: what's a good one?
         drift(1) = (zb(n)+znorm(conf))/dsqrt(t(n))
c     
         call bisect(n,t,za,zb,conf, drift(1))
         call glan(n,t,za,zb,drift(1), pstp,qpos,qneg,est,pr)
         call outgl(n,pr,drift(1),est,t,za,zb,pstp,qpos,qneg,qrcum)
      else if (nopt.eq.3) then
c        Probabilities from bounds, possibly with non-zero drift.
         do 50 i=1,ndrifs
            call glan(n,t,za,zb,drift(i), pstp,qpos,qneg,est,pr)
            call outgl(n,pr,drift(i),est,t,za,zb,pstp,qpos,qneg,qrcum)
            power(i) = qrcum
 50      continue
         if (ndrifs .gt. 1) then
c           Prompt for graphics.
            write(6,*)
            write(6,*) ' Do you want to see a graph? (1=yes,0=no)'
            read(5,100) char
 100        format(10a)
            if ((char.eq.'1').or.(char.eq.'y')) then
               call graphp(ndrifs,drift,power)
            end if
         end if
      else if (nopt.eq.4) then
c        Confidence limits from bounds and final statistics value.
         call ci(conf,n,t,za,zb, limit)
         call outci(conf,limit)
      end if
c     
c     Prompt for return to start
      write(6,*)
      write(6,*) ' Do you wish to start again? (1=yes,0=no)'
      read(5,200) char
 200  format(10a)
      if ((char.eq.'1').or.(char.eq.'y')) then
         go to 10
      end if
c     
      write(6,*) 
      write(6,*) ' Done.'
      write(6,*) 
c                                                                            
      stop                                                                  
      end                                                                    
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccc  I N P U T  R O U T I N E  cccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine input(nopt,n,ndim,t,iinf,t2,alpha,side,iuse,za,zb,
     .     ndrifs,drift,conf,pd,pe)
c
c     Interactive or noninteractive input of options, analysis times,
c     and other information.  Invalid input in noninteractive mode 
c     terminates the program unless a sensible default exists.  In
c     interactive mode, user is usually re-prompted.  Some types of
c     incorrect input cause the program to crash (eg character input 
c     for integer), so it's not completely robust.
c     
c     User inputs option, analysis times, spending function or bounds,
c     truncation value for bounds, optional drift parameters, power, or 
c     confidence level. 
c     
      double precision t(*),t2(*),za(*),zb(*),drift(*),pd(*),pe(*)
      double precision alpha, side, conf
      integer inter, nopt, n, ndim, iinf, iuse, ndrifs, idrif
c     
      double precision zninf, ztrun, value
      integer iequl,iside,ispnd,itrun,isymm,ifile,ibad
      character char
c     
      data inter/-1/
c     
c     Control for input from screen vs. input from a file.
c     Skip if the program is re-executing.
      if ((inter.ne.0).and.(inter.ne.1)) then
         write(6,*) ' Is this an interactive session? (1=yes,0=no)'
         read(5,50) char
 50      format(10a)
         if ((char.eq.'1').or.(char.eq.'y')) then
            inter = 1
         else if ((char.eq.'0').or.(char.eq.'n')) then
            inter = 0
         else
c           Invalid response terminates program!
            write(6,*) ' Invalid response.  Program terminating.'
            stop
         end if
         write(6,*) ' interactive = ', inter
      end if
c
c     Choose whether to compute bounds, probabilities or confidence
c     limits:
 100  if (inter.eq.1) then
         write(6,*) ' Enter number for your option:'
         write(6,*) ' (1) Compute bounds for given spending function.'
         write(6,*) ' (2) Compute drift for given power and bounds'
         write(6,*) ' (3) Compute probabilities for given bounds.'
         write(6,*) ' (4) Compute confidence interval.'
      end if
      read(5,*) nopt
      if (nopt.eq.1) then
         if (inter.eq.1) write(6,*) ' Option 1: You will be prompted',
     .        ' for a spending function.'
      else if (nopt.eq.2) then
         if (inter.eq.1) write(6,*) ' Option 2: You will be prompted',
     .        ' for bounds and a power level.'
      else if (nopt.eq.3) then 
         if (inter.eq.1) write(6,*) ' Option 3: You will be prompted',
     .        ' for bounds or a spending',' function to compute them.'
      else if (nopt.eq.4) then 
         if (inter.eq.1) write(6,*) ' Option 4: You will be prompted',
     .        ' for bounds and a confidence level.'
      else
         write(6,*) ' Invalid option.'
c        If running interactively, re-prompt.
         if (inter.eq.1) goto 100
c        If running noninteractively, program stops.
         stop
      end if
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccc INTERIM ANALYSIS TIMES cccccccccccccccccccccccccc
 200  if (inter.eq.1) write(6,*) ' Number of interim analyses?'
      read(5,*) n
      write(6,*) n,' interim analyses.'
c     Upper range check would be appropriate . . .
c     # of looks should be < ndim in calling routine.
      if (n.le.0) then
         write(6,*) ' Number must be positive.'
         if (inter.eq.1) goto 200
         stop
      else if (n.gt.ndim) then
         write(6,250) ndim
 250     format(' Number must be <=',i3,' or ndim must ',
     .        'be increased in main program.')
         if (inter.eq.1) goto 200
         stop
      end if
c     
c     Prompt for equally spaced analyses
 500  if (inter.eq.1) write(6,*) ' Equally spaced times between',
     .     ' 0 and 1? (1=yes,0=no)'
      read(5,510) char
 510  format(10a)
      if ((char.eq.'1').or.(char.eq.'y')) then
         iequl = 1
      else if ((char.eq.'0').or.(char.eq.'n')) then
         iequl = 0
      else
         write(6,*) ' Invalid response.'
         if (inter.eq.1) goto 500
c     Invalid responses cause manual entry
         iequl = 0
      end if
c     Compute equal spacing
 1000 if (iequl.eq.1) then
         do 1010 i=1,n
            t(i) = real(i)/real(n)
 1010    continue
         iequl = 0
      else
c        One time for each look.
         if (inter.eq.1) write(6,*) ' Times of interim analyses:',
     .        ' (>0 & <=1)'
         read(5,*)  (t(i),i=1,n)
      end if
c     
      write(6,1020) (t(i),i=1,n)
 1020 format(' Analysis times: ',25(f8.3,1x))
c     
c     Verify range and ordering of times
      ibad = 0
      if ((t(1).le.0.0).or.(t(n).gt.1.0)) ibad = 1
      do 1030 i=2,n
         if (t(i).le.t(i-1)) ibad = 1
 1030 continue
c     
      if (ibad.eq.1) then
c     Times out of range or order is bad.
         write(6,*) ' Times must be strictly increasing and <=1.'
         if (inter.eq.1) goto 1000
         stop
      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
cccc OPTIONS 2, 3 AND 4 DECIDE WHETHER TO USE SPENDING FCT ccccc
c     Users can enter either a spending function to compute 
c     bounds or enter a set of bounds manually.
      if (nopt.eq.1) then
         ispnd = 1
      else
 1100    if (inter.eq.1) write(6,*) ' Are you using a spending',
     .        ' function to determine bounds? (1=yes,0=no)'
         read(5,1110) char
 1110    format(10a)
         if ((char.eq.'1').or.(char.eq.'y')) then
            ispnd = 1
            write(6,*) ' Spending function will determine bounds.'
         else if ((char.eq.'0').or.(char.eq.'n')) then
            ispnd = 0
            write(6,*) ' You must enter a set of bounds.'
         else
            write(6,*) ' Invalid response.'
            if (inter.eq.1) goto 1100
c           Invalid responses are assumed to mean spending function!
            ispnd = 1
         end if
      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
cccccc DETERMINE BOUNDS USING SPENDING FUNCTION cccccccccccccccc
      if (ispnd.eq.1.) then
c     
         do 1200 i=1,n
            t2(i) = t(i)
 1200    continue
c     SECOND TIME SCALE
c        Will exact or estimated information be used?
         if (nopt.eq.1) then
            if (inter.eq.1)
     .           write(6,*) ' Do you wish to specify',
     .           ' a second time/information scale? (e.g.',
     .           ' number of patients or number of events,',
     .           ' as in Lan & DeMets 89?) ',
     .           ' (1=yes, 0=no)'
            read(5,1250) char
 1250       format(10a)
            if ((char.eq.'1').or.(char.eq.'y')) then
               iinf = 1
            else if ((char.eq.'0').or.(char.eq.'n')) then
               iinf = 0
            else
c           Invalid responses are assumed to mean time only!
               write(6,*) ' Invalid response.'
               iinf = 0
            end if
            if ((inter.eq.1).and.(iinf.eq.1)) 
     .           write(6,*) ' Second scale will estimate covariances.'
c     
c           If information is being used, enter.
            if (iinf.eq.1) then
               if (inter.eq.1) write(6,*) ' Information:'
               read(5,*) (t2(i),i=1,n)
               write(6,1260) (t2(i),i=1,n)
 1260          format(' Information ',25(f7.3,1x))
c     
c              Verify range and ordering of information
               ibad = 0
               if (t2(1).le.0.0) ibad = 1
               do 1270 i=2,n
                  if (t2(i).le.t2(i-1)) ibad = 1
 1270          continue
               if (ibad.eq.1) then
c                 Information <=0 or out of order is bad.
                  write(6,*) ' Information must be strictly',
     .                 ' increasing and > 0.  Input is ignored,',
     .                 ' computation uses first time scale only.'
                  do 1290 i=1,n
                     t2(i) = t(i)
 1290             continue
               end if
            end if
         end if
c     
c        Total type I error probability.
 2000    if (inter.eq.1) write(6,*) ' Overall significance level?',
     .        ' (>0 and <=1)'
         read(5,*) alpha
         if ((alpha.le.0.0).or.(alpha.gt.1.0)) then
            write(6,*) ' Invalid response.'
            write(6,*) ' Alpha must be > 0 and <= 1.'
            if (inter.eq.1) goto 2000
c           Out of range sets alpha = 0.05!
            alpha = 0.05
            write(6,*) ' Alpha has been set to 0.05.'
         end if
         write(6,2010) alpha
 2010    format('  alpha = ',f5.3)
c     
c        Upper and lower bounds, or upper only?
 2100    if (inter.eq.1) write(6,*) ' One(1) or two(2)-sided symmetric?'
         read(5,2110) char
 2110    format(10a)
         if ((char.eq.'1').or.(char.eq.'o')) then
            side = 1.d0
         else if ((char.eq.'2').or.(char.eq.'t')) then
            side = 2.d0
         else
            write(6,*) ' Invalid response.'
            if (inter.eq.1) goto 2100
c           Invalid responses are assumed to be two-sided!
            side = 2.d0
         end if
         write(6,2120) side
 2120    format(1x,f2.0,'-sided test')
c        
c        Alpha-star use functions can be added & edited in 
c        subroutine alphas below.
         if (inter.eq.1) then
            write(6,*) ' Use function? (1-5)'
            write(6,*) ' (1) OBrien-Fleming type'
            write(6,*) ' (2) Pocock type'
            write(6,*) ' (3) alpha * t'
            write(6,*) ' (4) alpha * t^1.5'
            write(6,*) ' (5) alpha * t^2'
         end if
         read(5,*)  iuse
         write(6,*) ' Use function alpha-star',iuse
c     
c     Does the design include a truncation point?
         if (inter.eq.1) write(6,*) ' Do you wish to truncate the',
     .        ' standardized bounds? (1=yes, 0=no)'
         read(5,2400) char
 2400    format(10a)
         if ((char.eq.'1').or.(char.eq.'y')) then
            itrun = 1
         else if ((char.eq.'0').or.(char.eq.'n')) then
            itrun = 0
         else
c     Invalid responses are assumed to mean no truncation!
            write(6,*) ' Invalid response.'
            itrun = 0
         end if
c     Note: ztrun must be > 0 or it is ignored.
         if (itrun.eq.1) then
            if (inter.eq.1) write(6,*) ' Enter truncation point: '
            read(5,*) ztrun
            write(6,2450) ztrun
 2450       format(' Bounds will be truncated at ',f5.2)
         else
            ztrun = -1.d0
c           write(6,*) ' Bounds will not be truncated.'
         end if
c     
c     Evaluate bounds by calling landem
         call landem(n,0.d0,alpha,t,t2,side,iuse,ztrun, pe,pd,za,zb)
c     
c     OPTION (1) RETURNS HERE
         if (nopt.eq.1) return
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc USERS ENTER THEIR OWN BOUNDS MANUALLY ccccccccccccccc
      else
c        User enters choice for one or two sided bounds.
 3000    if (inter.eq.1) write(6,*) ' One(1)- or two(2)-sided?'
         read(5,3010) char
 3010    format(10a)
         if ((char.eq.'1').or.(char.eq.'o')) then
            iside = 1
         else if ((char.eq.'2').or.(char.eq.'t')) then
            iside = 2
         else
            write(6,*) ' Invalid response.'
            if (inter.eq.1) goto 3000
c           Invalid responses are assumed to be two-sided!
            iside = 2
         end if
         write(6,*) ' ',iside,'-sided test'
c        
c        Two sided bounds may be symmetric or asymmetric.
 3100    if (iside.eq.2) then
            if (inter.eq.1) write(6,*) ' Symmetric bounds? (1=yes,0=no)'
            read(5,3120) char
 3120       format(10a)
            if ((char.eq.'1').or.(char.eq.'y')) then
               isymm = 1
            else if ((char.eq.'0').or.(char.eq.'n')) then
               isymm = 0
            else
               write(6,*) ' Invalid response.'
               if (inter.eq.1) goto 3100
c              Invalid responses are assumed to be symmetric
               isymm = 1
            end if
            if (isymm.eq.0) then
               write(6,*) ' Two sided asymmetric bounds.'
            else
               write(6,*) ' Two sided symmetric bounds.'
            end if
         else
            isymm = 1
            zninf = -8.d0
            write(6,*) ' One sided bounds.'
         end if
c        
c     Users may input bounds to a file, here 'userin', instead
c     of inputting from the screen or standard input.  This 
c     will usually require modification based on local system.
c***  open(unit=20,file='userin')
c      if (inter.eq.1) write(6,*) ' Will bounds be input from a',
c     .     ' file? (1=yes,0=no)'
c      read(5,3500) char
c 3500 format(10a)
c      if ((char.eq.'1').or.(char.eq.'y')) then
c         ifile = 1
c      else if ((char.eq.'0').or.(char.eq.'n')) then
c         ifile = 0
c      else
cc     Invalid responses are assumed to be no.
c         write(6,*) ' Invalid response.'
c         ifile = 0
c      end if
c        
         ifile = 0
         if (ifile.eq.0) then
            if (inter.eq.1) write(6,*) ' Enter upper bounds',
     .           ' (standardized): '
            read(5,*) (zb(i),i=1,n)
            if (iside.eq.2) then
               if (isymm.eq.0) then
                  if (inter.eq.1) write(6,*) ' Enter lower bounds',
     .                 ' (standardized): '
                  read(5,*) (za(i),i=1,n)
               else
                  do 3600 i=1,n
                     za(i) = -zb(i)
 3600             continue
               end if
            else
               do 3700 i=1,n
                  za(i) = zninf
 3700          continue
            end if
         else
            write(6,*) ' Bounds being entered from a file.'
            read(20,*) (zb(i),i=1,n)
            if (iside.eq.2) then
               if (isymm.eq.0) then
                  read(20,*) (za(i),i=1,n)
               else
                  do 3800 i=1,n
                     za(i) = -zb(i)
 3800             continue
               end if
            else
               do 3900 i=1,n
                  za(i) = zninf
 3900          continue
            end if
         end if
         write(6,*) ' Bounds entered.'
      end if
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c     Write out bounds
      write(6,*) '       Time                Bounds '
      do 4050 i=1,n
         write(6,4060) t(i),za(i),zb(i)
 4050 continue
 4060 format(4x,f5.2,12x,f7.4,4x,f7.4)
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
cccccccccccc FINDING DRIFT FOR GIVEN POWER ccccccccccccccccccccc
c     Default is 1 drift equal to zero
      ndrifs = 1
      drift(1) = 0.d0
c     
      if (nopt.eq.2) then
c        Enter desired power.
 5100    if (inter.eq.1) write(6,*) ' Desired power?',
     .        ' (>0 and <=1)'
         read(5,*) conf
         if ((conf.le.0.0).or.(conf.gt.1.0)) then
            write(6,*) ' Invalid response.'
            write(6,*) ' Power must be > 0 and <= 1.'
            if (inter.eq.1) goto 5100
c           Out of range sets conf = 0.80!
            conf = 0.80
            write(6,*) ' Power has been set to 0.80.'
         end if
         write(6,5110) conf
 5110    format('  Power is ',f5.3)
      end if
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccc COMPUTING POWER CURVES cccccccccccccccccccccccccccc
      if (nopt.eq.3) then
c     Are drift parameters to be used?
         if (inter.eq.1) write(6,*) ' Do you wish to use drift',
     .        '  parameters? (1=yes, 0=no)' 
         read(5,5700) char
 5700    format(10a)
         if ((char.eq.'1').or.(char.eq.'y')) then
            idrif = 1
         else if ((char.eq.'0').or.(char.eq.'n')) then
            idrif = 0
         else
c     Invalid responses are assumed to mean no drift parameters!
            write(6,*) ' Invalid response.'
            idrif = 0
         end if
         if (idrif.eq.1) then
 5720       if (inter.eq.1) write(6,*) ' How many drift parameters',
     .           ' do you wish to enter?'
            read(5,*) ndrifs
c     Upper range should be verified . . . 
            if (ndrifs.le.0) then
               write(6,*) ' Number of drifts must be positive'
               if (inter.eq.1) goto 5720
c     Program terminates.
               stop
            end if
            write(6,*) ndrifs,' drift parameters.'
c     
            if (inter.eq.1) write(6,*) ' Enter drift parameters: '
            read(5,*) (drift(i),i=1,ndrifs)
            write(6,5750) (drift(i),i=1,ndrifs)
 5750       format(' Drift parameters: ',25(f8.3,1x))
            write(6,*) ' Drift is equal to the standard',
     .           ' treatment difference times the square',
     .           ' root of total information per arm.'
         end if
      end if
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
cccccccc LEVEL FOR CONFIDENCE INTERVAL ccccccccccccccccccccccccc
      if (nopt.eq.4) then
c        Confidence interval replaces last bound with value
         if (inter.eq.1) write(6,*) ' Enter the standardized',
     .        ' statistic at the last analysis:'
         read(5,*) value
c     
         write(6,6770) value
 6770    format(4x,'Last value: ',15x,f7.4)
c        
c        Replace last bound with value
         zb(n) = value
c        
 6900    if (inter.eq.1) write(6,*) 
     .        ' Enter confidence level (>0 and <1):'
         read(5,*) conf
c        Range check
         if ((conf.le.0.0).or.(conf.gt.1.0)) then
            write(6,*) ' Invalid response.'
            write(6,*) ' Confidence level must be > 0 and < 1.'
            if (inter.eq.1) goto 6900
c           Out of range sets conf = 0.95!
            conf = 0.95
            write(6,*) ' Confidence level has been set to 0.95.'
         end if
         write(6,6910) 100.0*conf
 6910    format(f5.0,' percent confidence interval'/)
      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccc  O U T P U T  R O U T I N E S  ccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine outgl(n,pr,drift,est,t,za,zb,pstp,qpos,qneg,qrcum)
c     
c     Print out results from glan
c     Note: pr, est and pstp are included but not used.
c     
      double precision pr,drift,est,t(*),za(*),zb(*),
     .     pstp(*),qpos(*),qneg(*)
      integer n
c     
      double precision qplsr,qrcum
      integer i
c     
      write(6,10) n,drift
 10   format(/5x,'n =',i3,', drift =',f8.4)
c     
      write(6,20)
 20   format(/' look',5x,'time',6x,'lower',4x,'upper',4x,
     .     ' exit  probability ',3x,'cum exit pr'/)
c     .     'alpha(i)-alpha(i-1)',3x,'cum alpha'/)
c     .     'qplsr',5x,'qrcum'/)
c     .     'pdif',5x,'qplsr',5x,'pstp',5x,'qrcum'/)
c     
      qrcum = 0.d0
      do 100 i=1,n
         qplsr = qpos(i)+qneg(i)
         qrcum = qrcum + qplsr
c         if (i.eq.1) then
c            pdif = pstp(i)
c         else 
c            pdif = pstp(i)-pstp(i-1)
c         end if
c     
         write(6,110) i,t(i),za(i),zb(i),
     .        qplsr,qrcum
c     .        qpos(i),qrcum
c     .        pdif,qplsr,pstp(i),qrcum
 100  continue
 110  format(i3,4x,f7.2,2x, 2(2x,f7.4), 8x,f7.5,11x,f7.5)
      return
      end                                                                    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine outld(n,alpha,iuse,side,iinf,ti,t2,pd,za,zb)
c     
c     landem output for bounds computed from spending function
c     
      double precision alpha,side,ti(*),t2(*),pd(*),za(*),zb(*)
      integer n,iuse,iinf
c     
      double precision pr
      integer i
      character char
c     
c     Print out results.
      if (side .eq. 1.d0) then
         write(6,10)
         write(6,20) n,alpha,iuse
 10      format(//' This program generates one-sided boundaries.')
 20      format(' n = ',i2/' alpha =',f6.3/' use function = ',i1)
c     
      else
         write(6,30)
         write(6,40) n,alpha,iuse,iuse
 30      format(//' This program generates two-sided ',
     .        'symmetric boundaries.')
 40      format(' n = ',i2/' alpha =',f6.3
     .        /' use function for the lower boundary = ',i1
     .        /' use function for the upper boundary = ',i1)
      end if
c
c     Improve labelling?
      if (iinf.eq.0) then
         write(6,*) '       Time                Bounds ',
     .        '        alpha(i)-alpha(i-1)   cum alpha'
         pr = 0.d0
         do 100 i=1,n
            pr = pr + pd(i)
            write(6,150) ti(i),za(i),zb(i),pd(i),pr
 100     continue
 150     format(4x,f5.2,12x,f7.4,4x,f7.4,9x,f7.5,9x,f7.5)
c     
      else
         write(6,*) '    Time    Information     Bounds ',
     .        '     alpha(i)-alpha(i-1)  cum alpha'
         pr = 0.d0
         do 200 i=1,n
            pr = pr + pd(i)
            write(6,250) ti(i),t2(i),za(i),zb(i),pd(i),pr
 200     continue
 250     format(2x,f5.2,5x,f7.2,4x,f7.4,2x,f7.4,7x,f7.5,7x,f7.5)
      end if
c     
c     Prompt for graphics
      if (n .gt. 1) then
         write(6,*)
         write(6,*) ' Do you want to see a graph? (1=yes,0=no)'
         read(5,500) char
 500     format(10a)
         if ((char.eq.'1').or.(char.eq.'y')) then
            call graphb(n,ti,za,zb,side)
         end if
      end if
c     
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine outci(conf,limit)
      double precision conf,limit(*)
c     
      write(6,100) 100*conf,limit(1),limit(2)
 100  format(//2x,f3.0,' percent confidence interval:',
     .     ' (',f7.4,',' f7.4,')')
c     
      write(6,*) ' Drift is equal to the standard',
     .     ' treatment difference times the square',
     .     ' root of total information per arm.'
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc C O N F I D E N C E   I N T E R V A L   R O U T I N E  cccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ci(conf,n,t,za,zb, limit)
c     
c     Computes confidence intervals, per Kim & DeMets 87.
c     Given n interim analyses, bounds for the first n-1 and
c     the standardized statistic for the nth, the program finds
c     drift parameters satisfying confidence limit equations.
c
c     Search algorithm uses naive confidence limits
c     as starting values.  It brackets the limit and 
c     applies a bisection algorithm.
c     
c     Input variables:
c     conf is confidence level
c     n is the number of analyses.
c     t  is the vector of analysis time point on (0,1).
c     zb is the vector of upper bounds (standardized).
c     za is the vector of lower bounds (standardized).
c     limit(1) is lower confidence limit
c     limit(2) is upper confidence limit
c     
      double precision conf,t(*),za(*),zb(*),limit(*)
      integer n
c     
      double precision znorm
      double precision zcrit,target
      integer i
c     
c     Use naive limits as starting values
      zcrit = znorm(1.d0-(1.d0-conf)/2.d0)
      limit(1) = (zb(n) - zcrit)/dsqrt(t(n))
      limit(2) = (zb(n) + zcrit)/dsqrt(t(n))
c     
c     
      write(6,*) ' Starting computation for lower limit . . .'
      do 5000 i=1,2
         if (i.eq.2) write(6,*)
     .        ' Lower limit computed, starting on upper limit . . .'
c     
         target = real(i-1)*conf + (1.d0-conf)/2.d0
         call bisect(n,t,za,zb,target,limit(i))
c     
 5000 continue
c     
c     
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bisect(n,t,za,zb,target, drift)
c     Bisection routine to find drift resulting in exit probability.
c     drift should contain a decent nonzero starting value
c     
      double precision t(*),za(*),zb(*),target,drift
      integer n
c     
      double precision pstp(25),qpos(25),qneg(25),est,pr
      double precision lo,hi,prev
      integer gotlo, gothi
      integer j
c     
      double precision tol,del
      parameter(tol=1.d-6,del=0.25d0)
c     
c     Compute upper exit prob for drift, should be target
      gotlo = 0
      gothi = 0
      prev  = 0.d0
      pr = 0.d0
 1000 call glan(n,t,za,zb,drift, pstp,qpos,qneg,est,pr)      
      pr = 0.d0
      do 1300 j=1,n
c     
c     THIS OUGHT TO BE pr + pstp FOR POWER!  It works for
c     positive drifts since qneg is tiny but for drift < 0
c     this ought to fail miserably.
         pr = pr + qpos(j)
c     
 1300 continue
      
      if (dabs(pr-target).le.tol) then
c     pr is within tol of target
         go to 2000
      else if (dabs(drift-prev).le.tol/10.d0) then
c     If drift changes by less than tol/10, stop.
         write(6,*) ' Warning: Convergence problem.'
         go to 2000
c     
      else if (pr.gt.target+tol) then
c     Make sure hi gives pr .le. target + tol
         hi = drift
         drift = drift - del
         gothi = 1
      else 
c     Make sure lo gives pr .le. target + tol
         lo = drift
         drift = drift + del
         gotlo = 1
      end if
c     If bracketing values have been found, bisect
      if (gotlo*gothi .eq. 1) then
         prev  = drift
         drift = (lo+hi)/2.d0
      end if
      go to 1000
c     
 2000 continue
      return
      end
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccc S U B R O U T I N E  G L A N cccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine glan(nn,t,za,zb,drift, pstp,qpos,qneg,est,pr)
c     
c     Computes alpha level and other probabilities for the number of
c     tests of hypotheses.  Method used is a modification of Armitage,
c     McPherson and Rowe.  Consistent with Lan & DeMets 83, does not
c     assume equally spaced information times.
c
c     Input variables:
c     nn is the number of analyses.
c     t  is the vector of analysis time point on (0,1).
c     zb is the vector of upper bounds (standardized).
c     za is the vector of lower bounds (standardized).
c     drift is the noncentrality; at t=1, process has mean drift.
c
c     Returned variables:
c     pstp(i) is the prob of reaching ith analysis and stopping.
c     qpos(i) is the prob of reaching ith and exceeding upper.
c     qneg(i) is the prob of reaching ith and exceeding lower.
c     est is expected stopping time (not implemented).
c     pr is the total type I error used.
c
c     Local variables:
c     fn is the joint density for the current analysis.
c     last is the joint density at the last analysis.
c     yb is the vector of upper integration limits.
c     ya is the vector of lower integration limits.
c     stdv is the standard deviation of the increment
c           square root of t(i)-t(i-1).
c     sdproc is the standard deviation of the process.
c
c     Program parameters:
c     ndim is the maximum allowable number of analyses
c     h is the grid size for trapezoidal rule numerical integration
c     lnn dimensions local vectors in glan, NB if this is changed,
c        user must also change parameters in cprob, fcab and other!!!
c     
c     Subroutines and functions called directly:
c     pnorm is the standard normal cumulative distribution function.  
c     cprob computes exit probabilities for current analysis.
c     first computes joint density at first analysis.
c     other computes joint density at second and later analyses.
c     
      double precision t(*),za(*),zb(*),drift,
     .     pstp(*),qpos(*),qneg(*),est,pr
      integer nn
c     
      double precision h
      integer ndim,lnn
      parameter(ndim=25,lnn=5000,h=.05d0)
c     
      double precision ya(ndim),yb(ndim),stdv(ndim),sdproc(ndim),
     .     last(lnn)     
      integer nints(ndim)
      double precision pnorm
c
c     Check limit for number of interim analyses (looks).
      if(nn.gt.ndim) then
         write(6,10) nn,ndim
 10      format(/5x,' glan: No. of analyses too large, n =',i5,'>',i5)
         return
      end if
c                                                                           
c     Compute standard deviations of process and increments.
      call sd(nn,t, stdv,sdproc)
c     
c     Set upper (yb) and lower (ya) integration limits.  
c     Check limit for number of grid points.
      do 100 i=1,nn
         yb(i) = zb(i)*sdproc(i)-drift*t(i)
         ya(i) = za(i)*sdproc(i)-drift*t(i)
         nints(i) = idint((yb(i)-ya(i))/(h*stdv(i))) + 1
         if (nints(i).gt.lnn) then
            write(6,20) i,lnn
 20         format(5x,' glan: Bounds too wide, analysis',i3,' 
     .           (b-a)/h > ',i6/)
            return
         end if
 100  continue                                                              
c                                                                           
c     Begin loop calculating probabilities for each look.
c                                                                           
c     For first look some direct calculation is possible
c     using normal cdf.  Bounds in standard form are
c     adjusted for mean of process at t(1).
      pstp(1) = 1.d0 - (pnorm(zb(1)-drift*t(1)/stdv(1))
     .                - pnorm(za(1)-drift*t(1)/stdv(1)))
      qpos(1) = 1.d0 -  pnorm(zb(1)-drift*t(1)/stdv(1))
      qneg(1) =         pnorm(za(1)-drift*t(1)/stdv(1))
c     
c     After the first look, numerical integration is necessary.
      do 1000 i=2,nn                                                         
c     
c     Compute density of process at first look.
         if (i.eq.2) then
            call first(ya(1),yb(1),h,stdv(1), last,nints)
         end if
c     
c     Compute pstp, qpos and qneg in subroutine cprob.
         call cprob(last,nints,ya,yb,i,h,stdv(i),
     .        pstp(i),qpos(i),qneg(i))
         pstp(i) = 1.d0 - pstp(i)
c     
c     If density will be needed for next step, call other.
         if (i.ne.nn) then
            call other(ya,yb,i,stdv(i),h, last,nints)
         end if
c     
c     Accumulate boundary crossing probability.
         pr  = pr  + (qpos(i)+qneg(i))
c     
c     Expected stopping time.  (Not implemented.)
c        est = est + (qpos(i)+qneg(i))*(1.d0-t(i))
c     
 1000 continue
c     
c     Compute expected stopping time.  (Not implemented.)
c     est = 1.d0 - est
c     
      return
      end
c     
cccccccc SUBROUTINE SD cccccccccccccccccccccccccccccccccccccccccc
      subroutine sd(nn,t, stdv,sdproc)
c     
c     Standard deviations are calculated.  
c     
      double precision t(*),stdv(*),sdproc(*)
      integer nn
c     
      integer i
c     
c     Loop for standard deviations
c     Standard deviation for first analysis:
      stdv(1)   = dsqrt(t(1))
      sdproc(1) = stdv(1)
c     Standard deviation for second and later analyses:
      do 100 i=1,nn
         stdv(i) = dsqrt((t(i)-t(i-1)))
         sdproc(i) = dsqrt(t(i))
 100  continue
c     
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      subroutine cprob(last,nints,ya,yb,i,h,stdv, pstp,qpos,qneg)
c
c     Compute probabilities for current analysis.  Applying
c     Fubini's Theorem to Armitage, McPherson and Rowe's
c     formulae, the inner integral is a normal cdf and the outer
c     is computed by a trapezoidal rule numerical integration.
c
c     Input variables:
c     last is the vector of density values from the previous step.
c     ya is the vector of lower integration limits.
c     yb is the vector of upper integration limits.
c     i is the index of the current step.
c     h is the grid size (not used).
c     stdv is the standard deviation of the process increment.
c
c     Returned variables:
c     pstp is the value of the integral from ya to yb.
c     qpos is the value of the integral from yb to infinity.
c     qneg is the value of the integral from ya to minus infinity.
c     
c     Local variables:
c     pmid is vector of integrand values for pstp.
c     pupr is vector of integrand values for qpos.
c     plow is vector of integrand values for qneg.
c     grid is the argument for pmid, pupr and plow.
c     nints is the number of steps of size h between ya and yb.
c     hlast is the grid size from last step.
c     
c     Program parameters:
c     nwork dimensions pmid, pupr, plow.
c     
c     Subroutines and functions called directly:
c     pnorm is the standard normal cdf.
c     trap is a trapezoidal rule numerical integration routine.
c     
      double precision last(*),ya(*),yb(*),h,stdv,pstp,qpos,qneg
      integer i,nints(*)
c
      integer nwork
      parameter(nwork=5000)
c
      double precision pupr(nwork),plow(nwork),
     .     grid,tpstp,tqpos,tqneg,hlast
      integer j
      double precision pnorm
c
c     Previous grid size.
      hlast = (yb(i-1)-ya(i-1))/nints(i-1)
c     
c     Function values to be passed to numerical integration
c     routine are calculated.
      do 100 j=1,nints(i-1)+1
         grid = ya(i-1) + (hlast*(j-1))
c         pmid(j)  =        (pnorm((yb(i)-grid)/stdv)
c     .                    - pnorm((ya(i)-grid)/stdv)) * last(j)
         pupr(j)  = (1.d0 - pnorm((yb(i)-grid)/stdv)) * last(j)
         plow(j)  = (       pnorm((ya(i)-grid)/stdv)) * last(j)
 100  continue
c     
c     Calls to numerical integration routine.
c     pstp is not computed
c      call trap(pmid,nints(i-1),ya(i-1),yb(i-1),hlast, tpstp)
      tpstp = 0.d0
      call trap(pupr,nints(i-1),ya(i-1),yb(i-1),hlast, tqpos)
      call trap(plow,nints(i-1),ya(i-1),yb(i-1),hlast, tqneg)
      pstp = tpstp
      qpos = tqpos
      qneg = tqneg
c     
c     
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccc  S U B R O U T I N E   L A N D E M  ccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine landem(nn,drift,alpha,t,t2,side,iuse,ztrun, 
     .     pe,pd,za,zb)
c
c     Subroutine to calculate boundaries corresponding to
c     particular type I error spending rates and times of
c     interim analyses, per Lan & DeMets 83.
c
c     Input Variables:
c     nn is the number of interim analyses.
c     drift is the true mean of the process, drift of br.mot.
c       drift is set to zero, and so does not affect program.
c       NB drift MAY NOT be changed without adding code!
c     alpha is the total type I error
c     t is the information times
c     t2 are the numbers of events
c     side indicates one or two sided boundaries.
c     iuse indicates type I error spending function.
c     ztrun is the user selected truncation on integration limits.
c
c     Returned variables:
c     time is the expected first exit time.
c     pe is a vector of exit probabilities.
c     pd(i) = pe(i)-pe(i-1)
c     za is the vector of lower standardized boundaries
c     zb is the vector of upper standardized boundaries
c
c     Local variables:
c     stdv is the standard deviation of the process increment.
c     stproc is the standard deviation of the process.
c     last is the grid of the joint density at the last analysis.
c     nints is the number of intervals for numerical integration.
c     ya is the vector of lower integration limits.
c     yb is the vector of upper integration limits.
c     zninf is "negative infinity" for the program
c
c     Program parameters:
c     h is the grid size for numerical integration by trapezoidal rule.
c     lnn is the maximum dimension for last().
c     maxnn is the maximum number of interim analyses.
c
c     Subroutines and functions called directly:
c     znorm is the inverse standard normal cdf.
c     pnorm is the standard normal cdf.
c     alphas computes pe() according to specified use function.
c     search finds second and later integration bounds
c     first computes joint density of process at first analysis.
c     other computes joint density of process at second and later.
c     qpos computes tail probabilities.
c     
      double precision za(*),zb(*),alpha,t(*),t2(*),pe(*),pd(*),
     .     side,drift,ztrun,tol
      integer nn,iuse
c
      double precision h
      integer lnn,maxnn
      parameter(maxnn=50,lnn=5000,h=0.05d0)
c
      double precision stdv(maxnn),sdproc(maxnn),last(lnn),
     .     ya(maxnn),yb(maxnn),zninf
      integer i,nints(maxnn)
      double precision znorm, pnorm, qpos
c
c     
      zninf = -8.d0
c      write(6,10) zninf
c 10   format(' Negative infinity set to ',f7.3,/)
c     
c     Note: if ztrun <= 0 it is set to +infinity, i.e. ignored.
      if (ztrun.le.0.d0) ztrun = -zninf
c     
c     Tolerance for type I error to spend after 1st analysis.
      tol = 1.d-07
c     
c     Calculate probabilities according to use function.
      call alphas(iuse,nn,alpha,side,t,pe,pd)
c     
c     Calculate standard deviations of increments and process.
      call sd(nn,t2, stdv,sdproc)
c
c     Begin loop calculating boundaries.
c     
c     Direct calculations can be made for first analysis.
c     
c     Check type I error to spend
      if ((pd(1) .lt. 0.d0) .or. (pd(1) .gt. 1.d0)) then
         write(6,*) ' Error in spending function.'
         pd(1) = min(1.d0,pd(1))
         pd(1) = max(0.d0,pd(1))
      end if
c
c     Spending probability is zero (or less)
      if (pd(1) .eq. 0.d0) then
         zb(1) = -zninf
         if (zb(1) .gt. ztrun) then
            zb(1) = ztrun
            pd(1) = side * (1 - pnorm(zb(1)))
            pe(1) = pd(1)
            if (nn .gt. 1) then
               pd(2) = pe(2) - pe(1)
            end if
         end if
         yb(1) = zb(1)*stdv(1)
c     
c     Spending probability is one (or more)
      else if (pd(1) .eq. 1.d0) then
         zb(1) = 0.d0
         yb(1) = zb(1)*stdv(1)
      else
c
c     First bound based on normal distribution.
         zb(1) = znorm(1.d0-(pd(1)/side))
         if (zb(1) .gt. ztrun) then
            zb(1) = ztrun
            pd(1) = side * (1.d0 - pnorm(zb(1)))
            pe(1) = pd(1)
            if (nn .gt. 1) then
               pd(2) = pe(2) - pe(1)
            end if
         end if
         yb(1) = zb(1)*stdv(1)
      end if
c     
c     Lower bound is either "negative infinity" or -yb.
      if (side.eq.1.d0) then 
         za(1) = zninf
         ya(1) = za(1)*stdv(1)
      else
         za(1) = -zb(1)
         ya(1) = -yb(1)
      end if
c     
c     Number of intervals for numerical integration.
      nints(1) = idint((yb(1)-ya(1))/(h*stdv(1))) + 1
c     
c     Calculations for second and later analyses.
      do 100 i=2,nn
c
c     Calculate joint density for use in next step.
         if (i.eq.2) call first(ya(1),yb(1),h,stdv(1), last,nints)
c     
c     Check type I error to spend
      if ((pd(i) .lt. 0.d0) .or. (pd(i) .gt. 1.d0)) then
c        write(6,*) ' Error in spending function.'
         pd(i) = min(1.d0,pd(i))
         pd(i) = max(0.d0,pd(i))
      end if
c     
c     Spending probability is zero (or less than tol)
      if (pd(i) .lt. tol) then
         zb(i) = -zninf
         if (zb(i) .gt. ztrun) then
            zb(i) = ztrun
            pd(i) = side * qpos(zb(i)*sdproc(i),last,nints(i-1),ya(i-1), 
     .           yb(i-1),h,stdv(i)) 
            pe(i) = pd(i) + pe(i-1)
            if (nn .gt. i) then
               pd(i+1) = pe(i+1) - pe(i)
            end if
         end if
         yb(i) = zb(i)*sdproc(i)
c     
c     Spending probability is one (or more)
      else if (pd(i) .eq. 1.d0) then
         zb(i) = 0.d0
         yb(i) = zb(i)*stdv(i)
      else
c
c     Bounds are found using a search starting at the
c     bound from the previous analysis.
         call search(last,nints,i,h,pd(i)/side,stdv(i),ya, yb)
         zb(i) = yb(i)/sdproc(i)
         if (zb(i) .gt. ztrun) then
            zb(i) = ztrun
            pd(i) = side * qpos(zb(i)*sdproc(i),last,nints(i-1),ya(i-1),
     .           yb(i-1), h, stdv(i))   
            pe(i) = pd(i) + pe(i-1)
            if (nn .gt. i) then
               pd(i+1) = pe(i+1) - pe(i)
            end if
         end if
         yb(i) = zb(i)*sdproc(i)
      end if
c     
c     Lower bound is either "negative infinity" or -yb.
         if (side.eq.1.d0) then
            ya(i) = zninf*sdproc(i)
            za(i) = zninf
         else
            ya(i) = -yb(i)
            za(i) = -zb(i)
         end if
c     Number of intervals for numerical integration.
         nints(i) = idint((yb(i)-ya(i))/(h*stdv(i))) + 1
c     
c     Calculate joint density for use in next step.
         if (i.ne.nn) call other(ya,yb,i,stdv(i),h, last,nints)
 100  continue
c     
      return
      end
c     
      subroutine alphas(iuse,nn,alpha,side,t,pe,pd)
c
c     This use function, or type I error spending rate,
c     corresponds to the alpha 1 star in Lan & DeMets(1983).
c     Alpha 2 star in Lan & DeMets (1983).
c     Alpha 3 star in Lan & DeMets (1983).
c     Alpha 4 star in Kim & DeMets (1987).
c     Alpha 5 star in Kim & DeMets (1987).
c
c     Input variables:
c     iuse indicates type I error spending rate function.
c     nn is the number of interim analyses so far.
c     alpha is the desired overall size.
c     side is the number of sides for the test.
c     t is the vector of information times.
c     
c     Returned variables:
c     pe is the vector of type I error spent at each analysis.
c     pd is the differences intype I error spent between analyses.
c     
c     Local variables:
c     e is the base for natural logarithms.
c     
c     Subroutines and functions called directly:
c     pnorm is the standard normal cdf.
c     znorm is the inverse standard normal cdf.
c     
      double precision alpha,side,t(*),pe(*),pd(*)
      integer iuse,nn
c     
      integer i
      double precision pnorm,znorm
c     
      double precision e,tol
      e = 2.71828 18284 59045 23536d0
      tol = 1.d-13
c     
c     Calculate probabilities according to use function.
      do 50 i=1,nn
         if (iuse .eq. 1) then
            pe(i)=2.d0*
     .       (1.d0-pnorm(znorm(1.d0-(alpha/side)/2.d0)/dsqrt(t(i))))
         else if (iuse .eq. 2) then
            pe(i)=(alpha/side)*dlog(1.d0 + (e-1.d0)*t(i))
         else if (iuse .eq. 3) then
            pe(i)=(alpha/side)*t(i)
         else if (iuse .eq. 4) then
            pe(i)=(alpha/side)*(t(i) ** 1.5d0)
         else if (iuse .eq. 5) then
            pe(i)=(alpha/side)*(t(i) ** 2.0d0)
c     Add other spending function options here: e.g.
c        else if (iuse.eq.6) then . . .
         else
            write(6,*) ' Warning: invalid use function.'
         end if
c     
         pe(i)=side*pe(i)
c     
c     pd is the change in type I error spent
         if (i.eq.1) then
            pd(i) = pe(1)
         else
            pd(i) = pe(i) - pe(i-1)
         end if
c     Check type I error to spend
         if ((pd(i) .lt. 0.d0) .or. (pd(i) .gt. 1.d0)) then
            write(6,*) ' Error in spending function.'
            pd(i) = min(1.d0,pd(i))
            pd(i) = max(0.d0,pd(i))
         end if
c
         if (pd(i) .lt. tol) then
            write(6,*) ' Type I error spent too small, analysis', i
            write(6,*) ' Zero used as approximation for ',pd(i)
            pd(i) = 0.d0
         end if
c     
 50   continue
c     
      return
      end

      subroutine search(last,nints,i,h,pd,stdv,ya, yb)
c
c     A naive searching algorithm.  This starts at the previous
c     boundary, moves in the direction of the current boundary,
c     changes direction and reduces step size when target is
c     overstepped.  Stops when probability of estimate is within
c     eps of target probability.
c     
c     Users may want to substitute a better routine: some unused
c     variables are included to make this easier.
c
c     Input variables:
c     last is joint density from the previous analysis.
c     nints is the number of intervals for numerical integration.
c     i is the number of the current analysis.
c     h is the grid size (not used by current code).
c     pd is the target probability.
c     stdv is the standard deviation of the process increment.
c
c     Returned variables:
c     ya is the vector of lower integration limits.
c     yb is the vector of upper integration limits.
c
c     Local variables:
c     uppr is the estimate for yb(i).
c     q is the probability associated with uppr.
c     del is the searching step size.
c     eps is the tolerance for stopping.
c
c     Subroutines and functions called directly:
c     qpos computes tail probabilities.
c     
      double precision last(*),ya(*),yb(*),h,pd,stdv
      integer i,nints(*)
c     
      double precision uppr,del,eps,q
      double precision qpos
c     
c     Initialize tolerance for abs(q-pd).
      eps = 1.d-7
c
c     Initialize step size.
      del = 10.0d0
c
c     Initialize estimates at previous integration limit.
      uppr = yb(i-1)
      q = qpos(uppr,last,nints(i-1),ya(i-1),yb(i-1),h,stdv)
c
c     Begin searching.
 10   continue
c
c     If q and pd are nearly equal, set yb(i) and return.
      if (dabs(q-pd).le.eps) then
         yb(i) = uppr
         return
c
c     If q is too large, start increasing uppr by steps.
      else if (q.gt.(pd+eps)) then
c
c     Reduce step size by a factor of 10.
         del = del/10.d0
c
c     Increase uppr by del and check whether q is near pd.
         do 30 j=1,50
            uppr = uppr + del
            q = qpos(uppr,last,nints(i-1),ya(i-1),yb(i-1),h,stdv)
            if (q.le.(pd+eps)) go to 10
c     
c     If many iterations do not converge, print warning.
            if (mod(j,10).eq.0) write(6,90) i
 30      continue
c
c     After many steps, abort.
         write(6,100) i
         return
c
c     If q is too small, follow analogous procedure.
      else if (q.lt.(pd-eps)) then
         del = del/10.d0
         do 80 j=1,50
            uppr = uppr - del
            q = qpos(uppr,last,nints(i-1),ya(i-1),yb(i-1),h,stdv)
            if (q.ge.(pd-eps)) go to 10
            if (mod(j,10).eq.0) write(6,90) i
 80      continue
         write(6,100) i
         return
      end if
c     
 90   format('Large change in bounds, possible error analysis',i4)
 100  format('Error in search: not converging.  Abort analysis',i4)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     landem probability routines
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function qpos(xq,last,nints,yam1,ybm1,h,stdv)
c
c     Routine calculates tail probability by numerical integration.
c
c     Input variables:
c     xq is the argument for evaluation
c     last is the vector of density values from the previous step.
c     nints is the number of intervals for numerical integration.
c     yam1 is the lower integration limit from the previous step.
c     ybm1 is the upper integration limit from the previous step.
c     h is the grid size (not used).
c     stdv is the standard deviation of the process increment.
c     
c     Local variables:
c     fun is a temporary vector of function values.
c     grid is the step, yam1 + h(j-1)
c     
c     Program parameters:
c     nwork dimensions fun.
c     
c     Subroutines and functions called directly:
c     pnorm is the standard normal cdf.
c     trap is a trapezoidal rule numerical integration.
c     
      double precision xq,last(*),yam1,ybm1,h,stdv
      integer nints
c     
      integer nwork
      parameter(nwork=5000)
c     
      double precision fun(nwork),grid,tqpos,hlast
      integer j
      double precision pnorm
c     
c     Previous grid size.
      hlast = (ybm1-yam1)/nints
c     
c     Compute function at grid points.
      do 10 j=1,nints+1
         grid = yam1 + (hlast*(j-1))
         fun(j) = last(j)*(1.d0-pnorm((xq-grid)/stdv))
 10   continue
c     
c     Numerical integration.
      call trap(fun,nints,yam1,ybm1,hlast,tqpos)
      qpos = tqpos
c     
      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     numerical integration and normal probability routines
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine first(ya,yb,h,stdv, last,nints)
c
c     Joint density at first step is the normal.
c
c     Input variables:
c     ya is the lower integration limit for the first analysis.
c     yb is the upper integration limit for the first analysis.
c     h is the grid size, determines function argument (not used).
c     stdv is the standard deviation of the process at first analysis.
c     
c     Returned variables:
c     last is the vector of function values over the grid.
c     nints is the number of whole steps of size h between ya and yb.
c     
c     Local variables:
c     hh is the grid size used.
c     
c     Subroutines and functions called directly:
c     g is the standard normal density function.
c     
      double precision ya,yb,last(*),stdv,h
      integer nints(*)
c     
      double precision hh
      integer j
      double precision g
c     
c     Compute grid size to be used.
      hh = (yb-ya)/nints(1)
c     
c     Evaluate function (normal density) at grid points.
      do 10 j=1,nints(1)+1
         last(j) = g((ya+(hh*(j-1)))/stdv)/stdv
 10   continue
c     
      return
      end

      subroutine other(ya,yb,i,stdv,h, last,nints)
c
c     Computes joint density for current step using fcab.
c
c     Input variables:
c     ya is the vector of lower integration limits.
c     yb is the vector of upper integration limits.
c     h is the grid size, determines function argument (not used).
c     stdv is the standard deviation of the process increment.
c     last on entry is the density values from the previous step.
c     
c     Returned variables:
c     last is the vector of function values over the grid.
c     nints is the number of whole steps of size h between ya and yb.
c     
c     Local variables:
c     grid is the argument for function evaluation.
c     fn is a temporary vector of function values--last from
c      previous step is needed to compute each current value.
c     hh is the actual grid size used.
c     
c     Program parameters:
c     nlim dimensions fn.  Must be compatible with dimension of last.
c     
c     Subroutines and functions called directly:
c     fcab computes current values of the joint density.
c     
      double precision ya(*),yb(*),stdv,h,last(*)
      integer i,nints(*)
c     
      integer nlim
      parameter(nlim=5000)
c     
      double precision grid,fn(nlim),hh,hlast
      integer j
      double precision fcab
c     
c     Current grid size to be used.
      hh = (yb(i)-ya(i))/nints(i)
c     
c     Previous grid size.
      hlast = (yb(i-1)-ya(i-1))/nints(i-1)
c     
c     Evaluate function over grid ya + (j-1)*hh, j=1,nints+1
      do 10 j=1,nints(i)+1
         grid = ya(i) + (hh*(j-1))
         fn(j) = fcab(last,nints(i-1),ya(i-1),yb(i-1),hlast, grid,stdv)
 10   continue
c     
c     Copy fn into last.
      do 20 j=1,nints(i)+1
         last(j) = fn(j)
 20   continue
c     
      return
      end

      double precision function fcab(last,nints,yam1,ybm1,h, x,stdv)
c
c     This function computes the density function of (W(t(i)),t(i)) 
c     where t(i) are the stopping times and W(t(i)) is the standard
c     Brownian motion process on (0,1).
c
c     The joint density for the previous step is multiplied by
c     the density for the increment to the current step (just the
c     normal density with proper mean and variance), and the
c     result is numerically integrated.  
c
c     Input variables:
c     last is the vector of values of the previous density.
c     nints is the number of steps of size h between yam1 and ybm1.
c     yam1 is the lower integration limit for the previous step.
c     ybm1 is the upper integration limit (not used).
c     h is the grid size.
c     x is the argument at which density function is evaluated.
c     stdv is the square root of the diffusion of a Brownian motion.
c     
c     Local variables:
c     f is a temporary vector of function values to pass to trap.
c     
c     Program parameters:
c     nlim dimensions f.  Must be compatible with dimension of last.
c     
c     Subroutines and functions called directly:
c     g is the standard normal density function.
c     trap is a trapezoidal rule numerical integration.
c     
      double precision last(*),yam1,ybm1,h,x,stdv
      integer nints
c     
      integer nlim
      parameter(nlim=5000)
c     
      double precision f(nlim),grid,tfcab
      integer j
      double precision g
c     
c     Evaluate function over grid
      do 10 j=1,nints+1
         grid = yam1 + (h*(j-1))
         f(j) = last(j) * (g((x-grid)/stdv)/stdv)
 10   continue
c     
c     Numerical integration.
      call trap(f,nints,yam1,ybm1,h, tfcab)
      fcab = tfcab
c     
      return
      end
c     

      subroutine trap(f,n,a,b,h, area)
c
c     Numerical integration by composite trapezoidal rule.
c
c     Input variables:
c     f is a vector of function values.
c     n is the number of trapezoids with base h
c     a is the lower integration limit, the first grid point (not used).
c     b is the upper integration limit, the last grid point (not used).
c     h is the grid size
c     
c     Returned variables:
c     area is the sum of the specified function values
c     
c     Local variables:
c     sum accumulates the area of the trapezoids
c     
      double precision f(*),a,b,h,area
      integer n
c     
      double precision sum
c     
c     Initialize sum with the first function value.
      sum = f(1)
c     
c     All but the first and last function values appear
c      twice in the summation.
      do 100 j=2,n
         sum = sum + 2.d0*f(j)
 100  continue
c
c     Add the last function value.
      sum = sum + f(n+1)
c
c     Multiply the sum by h/2
      area = (h/2.d0)*sum 
c     
      return
      end
c     
      double precision function g(x)
c     g(x) is the normal (0,1) density function
      double precision x,pi
      parameter(pi=3.1415926536d0)
      g=dexp(-0.5d0*x*x)/dsqrt(2.d0*pi)
      return                                                                
      end
c     
c     
      double precision function znorm(p)
c     
c     znorm is the inverse normal cdf.
c     From statlib@temper.stat.cmu.edu ("send 241 from apstat") with 
c     only minor modification (e.g. IFAULT not returned).
c     
C       ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
C
C       Produces the normal deviate Z corresponding to a given lower
C       tail area of P; Z is accurate to about 1 part in 10**16.
C
C       The hash sums below are the sums of the mantissas of the
C       coefficients.   They are included for use in checking
C       transcription.
C
        DOUBLE PRECISION ZERO, ONE, HALF, SPLIT1, SPLIT2, CONST1,
     *          CONST2, A0, A1, A2, A3, A4, A5, A6, A7, B1, B2, B3,
     *          C0, C1, C2, C3, C4, C5, C6, C7, D1, D2, D3, D4, D5,
     *          D6, D7, E0, E1, E2, E3, E4, E5, E6, E7, F1, F2, F3,
     *          F4, F5, F6, F7, P, Q, R
        PARAMETER (ZERO = 0.D0, ONE = 1.D0, HALF = 0.5D0,
     *          SPLIT1 = 0.425D0, SPLIT2 = 5.D0,
     *          CONST1 = 0.180625D0, CONST2 = 1.6D0)
C
C       Coefficients for P close to 0.5
C
        PARAMETER (A0 = 3.38713 28727 96366 6080D0,
     *             A1 = 1.33141 66789 17843 7745D+2,
     *             A2 = 1.97159 09503 06551 4427D+3,
     *             A3 = 1.37316 93765 50946 1125D+4,
     *             A4 = 4.59219 53931 54987 1457D+4,
     *             A5 = 6.72657 70927 00870 0853D+4,
     *             A6 = 3.34305 75583 58812 8105D+4,
     *             A7 = 2.50908 09287 30122 6727D+3,
     *             B1 = 4.23133 30701 60091 1252D+1,
     *             B2 = 6.87187 00749 20579 0830D+2,
     *             B3 = 5.39419 60214 24751 1077D+3,
     *             B4 = 2.12137 94301 58659 5867D+4,
     *             B5 = 3.93078 95800 09271 0610D+4,
     *             B6 = 2.87290 85735 72194 2674D+4,
     *             B7 = 5.22649 52788 52854 5610D+3)
C       HASH SUM AB    55.88319 28806 14901 4439
C
C       Coefficients for P not close to 0, 0.5 or 1.
C
        PARAMETER (C0 = 1.42343 71107 49683 57734D0,
     *             C1 = 4.63033 78461 56545 29590D0,
     *             C2 = 5.76949 72214 60691 40550D0,
     *             C3 = 3.64784 83247 63204 60504D0,
     *             C4 = 1.27045 82524 52368 38258D0,
     *             C5 = 2.41780 72517 74506 11770D-1,
     *             C6 = 2.27238 44989 26918 45833D-2,
     *             C7 = 7.74545 01427 83414 07640D-4,
     *             D1 = 2.05319 16266 37758 82187D0,
     *             D2 = 1.67638 48301 83803 84940D0,
     *             D3 = 6.89767 33498 51000 04550D-1,
     *             D4 = 1.48103 97642 74800 74590D-1,
     *             D5 = 1.51986 66563 61645 71966D-2,
     *             D6 = 5.47593 80849 95344 94600D-4,
     *             D7 = 1.05075 00716 44416 84324D-9)
C       HASH SUM CD    49.33206 50330 16102 89036
C
C       Coefficients for P near 0 or 1.
C
        PARAMETER (E0 = 6.65790 46435 01103 77720D0,
     *             E1 = 5.46378 49111 64114 36990D0,
     *             E2 = 1.78482 65399 17291 33580D0,
     *             E3 = 2.96560 57182 85048 91230D-1,
     *             E4 = 2.65321 89526 57612 30930D-2,
     *             E5 = 1.24266 09473 88078 43860D-3,
     *             E6 = 2.71155 55687 43487 57815D-5,
     *             E7 = 2.01033 43992 92288 13265D-7,
     *             F1 = 5.99832 20655 58879 37690D-1,
     *             F2 = 1.36929 88092 27358 05310D-1,
     *             F3 = 1.48753 61290 85061 48525D-2,
     *             F4 = 7.86869 13114 56132 59100D-4,
     *             F5 = 1.84631 83175 10054 68180D-5,
     *             F6 = 1.42151 17583 16445 88870D-7,
     *             F7 = 2.04426 31033 89939 78564D-15)
C       HASH SUM EF    47.52583 31754 92896 71629
C
        IFAULT = 0
        Q = P - HALF
        IF (ABS(Q) .LE. SPLIT1) THEN
          R = CONST1 - Q * Q
          znorm = Q * (((((((A7 * R + A6) * R + A5) * R + A4) * R + A3)
     *                  * R + A2) * R + A1) * R + A0) /
     *                (((((((B7 * R + B6) * R + B5) * R + B4) * R + B3)
     *                  * R + B2) * R + B1) * R + ONE)
          RETURN
        ELSE
          IF (Q .LT. ZERO) THEN
            R = P
          ELSE
            R = ONE - P
          END IF
          IF (R .LE. ZERO) THEN
            IFAULT = 1
            znorm = ZERO
            RETURN
          END IF
          R = SQRT(-LOG(R))
          IF (R .LE. SPLIT2) THEN
            R = R - CONST2
            znorm = (((((((C7 * R + C6) * R + C5) * R + C4) * R + C3)
     *                  * R + C2) * R + C1) * R + C0) /
     *               (((((((D7 * R + D6) * R + D5) * R + D4) * R + D3)
     *                  * R + D2) * R + D1) * R + ONE)
          ELSE
            R = R - SPLIT2
            znorm = (((((((E7 * R + E6) * R + E5) * R + E4) * R + E3)
     *                  * R + E2) * R + E1) * R + E0) /
     *               (((((((F7 * R + F6) * R + F5) * R + F4) * R + F3)
     *                  * R + F2) * R + F1) * R + ONE)
          END IF
          IF (Q .LT. ZERO) znorm = - znorm
          RETURN
        END IF
        END
c     
c
        double precision function pnorm(Z)
c     
c     pnorm is the normal cdf.
c     From statlib@temper.stat.cmu.edu ("send 66 from apstat") with
c     changes: subroutine to function.
c     
C       REFERENCE: ADAMS,A.G. AREAS UNDER THE NORMAL CURVE,
C       ALGORITHM 39, COMPUTER J., VOL. 12, 197-8, 1969.
C
C       P, Q = PROBABILITIES TO THE LEFT AND RIGHT OF Z
C       FOR THE STANDARD NORMAL DISTRIBUTION.
C       PDF  = THE PROBABILITY DENSITY FUNCTION
C
C       LATEST REVISION - 23 JANUARY 1981
C
C********************************************************************
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DATA A0,A1,A2,A3,A4,A5,A6,A7/0.5D0, 0.398942280444D0,
     1  0.399903438504D0, 5.75885480458D0, 29.8213557808D0,
     2  2.62433121679D0, 48.6959930692D0, 5.92885724438D0/,
     3  B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11/0.398942280385D0,
     4  3.8052D-8, 1.00000615302D0, 3.98064794D-4, 1.98615381364D0,
     5  0.151679116635D0, 5.29330324926D0, 4.8385912808D0,
     6  15.1508972451D0, 0.742380924027D0, 30.789933034D0,
     7  3.99019417011D0/
C
        ZABS = ABS(Z)
        IF(ZABS.GT.12.7D0) GO TO 20
        Y = A0*Z*Z
        PDF = EXP(-Y)*B0
        IF(ZABS.GT.1.28D0) GO TO 10
C
C       Z BETWEEN -1.28 AND +1.28
C
        Q = A0-ZABS*(A1-A2*Y/(Y+A3-A4/(Y+A5+A6/(Y+A7))))
        IF(Z.LT.0.D0) GO TO 30
        pnorm = 1.D0-Q
        RETURN
C
C       ZABS BETWEEN 1.28 AND 12.7
C
   10   Q = PDF/(ZABS-B1+B2/(ZABS+B3+B4/(ZABS-B5+B6/(ZABS+B7-B8/
     1  (ZABS+B9+B10/(ZABS+B11))))))
        IF(Z.LT.0.D0) GO TO 30
        pnorm = 1.D0-Q
        RETURN
C
C       Z FAR OUT IN TAIL
C
   20   Q = 0.D0
        PDF = 0.D0
        IF(Z.LT.0.D0) GO TO 30
        pnorm = 1.D0
        RETURN
C
C       NEGATIVE Z, INTERCHANGE P AND Q
C
 30     pnorm = Q
c        Q = 1.D0-P
        RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine graphp(n,x,y)
c     
c     Specialized graphics routine to plot power function.  Assumes
c     scale for x-axis.
c     
c     Input variables:  
c     n is the length of x and y
c     x is locations of points on the x axis
c     y is location of points on the y axis
c     
c     Local variables:
c     hspace is # of x-axis spaces
c     vspace is # of y-axis spaces, multiple of 5
c     symbol is the plotting character on the graph                     
c     xmax,xmin are the maximum and minimum X axis values 
c     ymax,ymin are the maximum and minimum Y axis values 
c     dely is the value of one space on Y axis
c     window is a two-dimensional character array
c     xloc,yloc are the locations of a point in the window               
c     count indicates units on Y axis for labeling
c     
      double precision x(*),y(*)
      integer n
c     
      double precision ymax,ymin,xmax,xmin
      double precision dely,ylab
      integer i,xloc,yloc
      integer hspace,vspace
      parameter(hspace=42,vspace=25)
      character*1 window(0:hspace,0:vspace),xaxis(0:hspace),symbol
      parameter(symbol='*')
c     
      character*1 cints(0:9)
      data cints/'0','1','2','3','4','5','6','7','8','9'/
c     
c     Initialize window and xlabels to blanks
      do 150 xloc=0,hspace
         do 120 yloc=0,vspace
            window(xloc,yloc)=' '
 120     continue
         xaxis(xloc)=' '
 150  continue
c
c     Maximum and minimum for X axis and Y axis             
      ymax=y(1)
      ymin=y(1)
      do 200 i=2,n
         if(y(i) .gt. ymax) ymax=y(i)
         if(y(i) .lt. ymin) ymin=y(i)
 200  continue
c     
      xmax = x(1)
      xmin = x(1)
      do 210 i=2,n
         if(x(i) .gt. xmax) xmax = x(i)
         if(x(i) .lt. xmin) xmin = x(i)
 210  continue
c
c     Determine limits for x axis: use closest integers
      if (xmin .lt. 0.0) then
         xmin = aint(xmin) - 1.0
      else
         xmin = aint(xmin)
      end if 
      xmax = aint(xmax) + 1
c     
c     Determine limits for y axis: 
c       either 0 to .5, .5 to 1, or 0 to 1
      if (ymin .lt. 0.5) then
         ymin = 0.0
      else
         ymin = 0.5
      end if
      if (ymax .gt. 0.5) then 
         ymax = 1.0
      else
         ymax = 0.5
      end if
c     
c     Determine and mark window locations for each (x,y)
      do 500 i=1,n
c        No rocket science here.
         xloc = anint(((x(i)-xmin)/(xmax-xmin))*hspace)
         yloc = anint(((y(i)-ymin)/(ymax-ymin))*vspace)

         window(xloc,yloc)=symbol
 500  continue
c     
c     Compute numeric value of each y space in window
      dely=dabs((ymax-ymin)/real(vspace))
c     
c     Output graph: 2 blank lines, one with axis extension
      write(6,*)
      write(6,650)
c     Now start writing out lines, labelling 
      ylab = ymax
      do 610 yloc = vspace,0,-1
         write(6,640) ylab, (window(xloc,yloc),xloc=0,hspace)
         ylab = ylab - dely
 610  continue
 640  format(1x,f8.2,':',1x,100a)
 650  format(9x,':',1x,100a)
c     
c     Write out x-axis
      write(6,700) ('.',xloc=0,hspace+4)
 700  format(9x,100a)
c     
c     Lable x axis: this is a little tricky
c     We only label integer values of x, and we do it by
c     creating a character vector with the integers in it.
c     
c     First label, last label
      if (xmin .lt. 0) then
         xaxis(0) = '-'
         i = -xmin
         xaxis(1) = cints(mod(i,10))
      else
         i = xmin
         xaxis(0) = cints(mod(i,10))
      end if
      if (xmax .lt. 0) then
         xaxis(hspace-1) = '-'
         i = int(-xmax)
         xaxis(hspace) = cints(mod(i,10))
      else
         i = int(xmax)
         xaxis(hspace) = cints(mod(i,10))
      end if
c     
c     Insert other labels
      do 720 i=xmin+1,xmax-1
         xloc = anint(((real(i)-xmin)/(xmax-xmin))*hspace)
c     
         if (i .lt. 0) then
            xaxis(xloc-1) = '-'
            xaxis(xloc) = cints(mod(-i,10))
         else
            xaxis(xloc) = cints(mod(i,10))
         end if
 720  continue
c     
      write(6,730) (xaxis(xloc),xloc=0,hspace)
 730  format(11x,100a)
c     
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine graphb(n,x,y1,y2,side)
c     
c     Specialized graphics routine to plot bounds.
c     
c     Input variables:  
c     n is the length of x and y
c     x is points on the x axis
c     y1 is lower boundary
c     y2 is upper boundary
c     side indicates 1 or 2 sided bounds
c     
c     Local variables:
c     hspace is # of x-axis spaces
c     vspace is # of y-axis spaces
c     symbol is the plotting character on the graph                     
c     xmax,xmin are the maximum and minimum X axis values 
c     ymax,ymin are the maximum and minimum Y axis values 
c     dely is the value of one space on Y axis
c     window is a two-dimensional character array
c     xloc,yloc are the locations of a point in the window               
c     count indicates units on Y axis for labeling
c     
c     ylim is the upper limit for the y-axis labels
c     flag indicates points outside plotting range
c     
      double precision x(*),y1(*),y2(*),side
      integer n
c     
      double precision ymax,ymin,xmax,xmin
      double precision dely,ylab,ylim
      integer i,xloc,yloc1,yloc2,flag
      integer hspace,vspace
      parameter(hspace=42,vspace=25)
      character*1 window(0:hspace,0:vspace),xaxis(0:hspace),symbol
      parameter(symbol='*')
c     
      character*1 cints(0:9)
      data cints/'0','1','2','3','4','5','6','7','8','9'/
c     
c     Initialize y-axis label limit (always > 0) and flag
      ylim = 4.0
      flag = 0
c     
c     Initialize window and xlabels to blanks
      do 150 xloc=0,hspace
         do 120 yloc=0,vspace
            window(xloc,yloc)=' '
 120     continue
         xaxis(xloc)=' '
 150  continue
c     
      if (side .eq. 1.d0) then
         do 180 i=1,n
            y1(i) = 0.01
 180     continue
      end if
c     
c     Determin maximum and minimum values in y1 and y2
c     NOTE: y2 >= y1!
      ymax=y2(1)
      ymin=y1(1)
      do 200 i=2,n
         if(y2(i) .gt. ymax) ymax=y2(i)
         if(y1(i) .lt. ymin) ymin=y1(i)
 200  continue
c     
c     Determine limits for y axis: use closest integers
      if (ymin .lt. 0.0) then
         ymin = aint(ymin) - 1.0
      else
         ymin = aint(ymin)
      end if 
      ymax = aint(ymax) + 1.0
c     
c     If limits are beyond ylim, some points will be off plot
      if (ymin .lt. -ylim) ymin = -ylim
      if (ymax .gt.  ylim) ymax =  ylim
c     
c     Determine limits for x axis: 
c       either 0 to .5, .5 to 1, or 0 to 1
      if (x(1) .lt. 0.5) then
         xmin = 0.0
      else
         xmin = 0.5
      end if
      if (x(n) .gt. 0.5) then 
         xmax = 1.0
      else
         xmax = 0.5
      end if
c     
c     Determine and mark window locations for each (x,y)
      do 500 i=1,n
         xloc = anint(((x(i)-xmin)/(xmax-xmin))*hspace)
c     
c        If y value is outside window, place it at +/- ylim,
c        which is 0 or vspace for yloc, and designate it
c        with a special symbol, like "^" or "v".  And set flag.
         if (y1(i) .lt. ymin) then 
            yloc1 = 0
            if (side .eq. 2.d0) window(xloc,yloc1)='v'
            flag = 1
         else
            yloc1 = anint(((y1(i)-ymin)/(ymax-ymin))*vspace)
            if (side .eq. 2.d0) window(xloc,yloc1)=symbol
         end if
c     
         if (y2(i) .gt. ymax) then
            yloc2 = vspace
            window(xloc,yloc2)='^'
            flag = 1
         else
            yloc2 = anint(((y2(i)-ymin)/(ymax-ymin))*vspace)
            window(xloc,yloc2)=symbol
         end if
 500  continue
c     
c     Compute numeric value of each y space in window
      dely=dabs((ymax-ymin)/real(vspace))
c     
c     Output graph: Start with one blank line
      write(6,*)
c
c     Add a title
      write(6,*) '                  GROUP SEQUENTIAL BOUNDARIES'
      write(6,*)
c
c     Extend the top of the y-axis with one line, add y title
      write(6,610)
 610  format(3x,'Z',5x,':',1x,100a)
c     
c     Now start writing out lines, labelling 
      ylab = ymax
      do 620 yloc = vspace,0,-1
         write(6,630) ylab, (window(xloc,yloc),xloc=0,hspace)
         ylab = ylab - dely
 620  continue
 630  format(1x,f8.2,':',1x,100a)
c     
c     Write out x-axis
      write(6,700) ('.',xloc=0,hspace+4)
 700  format(9x,100a)
c     
c     Lable x axis: this is a little tricky.
c     cints is a character with digits 0 to 9 which
c     will be the labels after a decimal point is added.
c     
c     First label, last label
      if (xmin .lt. 0.5) then
c        Should this be i = int(xmin)?
         i = xmin
         xaxis(0) = cints(mod(i,10))
      else
         i = 5
         xaxis(0) = '.'
         xaxis(1) = cints(i)
      end if
      if (xmax .lt. 1.0) then
         xaxis(hspace-1) = '.'
         i = 5
         xaxis(hspace) = cints(mod(i,10))
      else
         i = int(xmax)
         xaxis(hspace) = cints(mod(i,10))
      end if
c     
c     Insert other labels
      xmin = 10*xmin
      xmax = 10*xmax
      do 720 i=xmin+1,xmax-1
         xloc = anint(((real(i)-xmin)/(xmax-xmin))*hspace)
c     
         xaxis(xloc-1) = '.'
         xaxis(xloc) = cints(mod(i,10))

 720  continue
c     
      write(6,730) (xaxis(xloc),xloc=0,hspace)
 730  format(11x,100a)
c     
c     Title for x-axis
      write(6,*) '                     Information Fraction'
      write(6,*)
c     
c     Warning for special symbols
      if (flag .eq. 1) write(6,800) ylim
 800  format(10x, 'Symbols ^ and v indicate bounds beyond +/- ',f4.1)
c     
      return
      end

