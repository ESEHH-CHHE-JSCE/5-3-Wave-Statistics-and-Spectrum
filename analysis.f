c *** Execution of zero up-cross method and Fourier transform **********
c                                2022/10/08 Produced by Katsuya Hirayama
c **********************************************************************

      parameter( kmax=7200, ichmax=1 )
      character*80 title
      character*40 fname
      integer ke( kmax )
c
      real eta( kmax ), spc( 0:kmax/2 )
      real fq( 0:kmax/2 )
      real reph( 4 ), rept( 4 ), resh( 4 )
c
c
      fname = '901.txt'  !時系列データ名
      ichr = 7           !データ名の文字数（半角英数字）
      n  = 2400          !データ総数
      dt = 0.5           !データサンプリング時間[s]
      nf = 25            !放物線フィルターの幅
c
      open(20,file=fname(1:ichr),status='old')
      read(20,'(10x,i5)') kch
      do i = 1, n
        read(20,'(10x,f10.4)') eta(i)
      end do
      close(20)
c
c --- 平均水深の算定と水面変動データの作成 ---
      sume = 0.
      do i = 1, n
        sume = sume + eta(i)
      end do
      eavg = sume/real(n)
      do i = 1, n
        eta(i) = eta(i) - eavg
      end do
c ----------------------------------------------------------
c
c --- スペクトル波高の算定 ------
      sigma = 0.
      do i = 1, n
        sigma = sigma + eta(i)**2
      end do
      sigma = sigma/real(n)
      
      resh(1) = 5.090*sqrt(sigma)
      resh(2) = 4.004*sqrt(sigma)
      resh(3) = resh(2)/1.597
      resh(4) = sqrt(8.*sigma)
c -------------------------------
c
          call repwv( eta,n,dt,reph,rept,nowv,fname,ichr )
      nb = n
	  call spcwv( eta,n,nb,spc,sigma,dt,df,fq,nf,fname,ichr )
	  call parspc( spc,n,df,fq,f1,f2,rnp,ep )
c
      write( *,'(a40,f10.4)') fname,dt
      write( *,'(a4,a10,a5,6a10)')  '  CH','    MSL[m]','  WAV',
     &           '  H1/10[m]','   H1/3[m]','   Hbar[m]',
     &           '  T1/10[s]','   T1/3[s]','   Tbar[s]'
      write( *,'(i4,f10.4,i5,6f10.4)') kch,eavg,nowv,(reph(l),l=2,4),
     &                                               (rept(l),l=2,4)
c
      open(30,file=fname(1:ichr-4)//'_tok'//'.txt',status='unknown')
      write(30,'(a40,f10.4)') fname,dt
      write(30,'(a5,a10,a5,16a10)') '   CH','    MSL[m]','  WAV',
     &            '   Hmax[m]','  H1/10[m]','   H1/3[m]','   Hbar[m]',
     &            '   Tmax[s]','  T1/10[s]','   T1/3[s]','   Tbar[s]',
     &            ' SH1/10[m]','  SH1/3[m]','  SHbar[m]','  SHrms[m]',
     &            '    T01[s]','    T02[s]','       NYU','       EPS'
      write(30,'(i5,f10.4,i5,16f10.4)') kch,eavg,nowv,(reph(l),l=1,4),
     &             (rept(l),l=1,4),(resh(l),l=1,4),1./f1,1./f2,rnp,ep
      close(30)
c
      open(50,file=fname(1:ichr-4)//'_spc'//'.txt',status='unknown')
      write(50,'(a40)') fname
      write(50,'(2a10)') '    f [Hz]','  Sp [m2s]'
      do k = 0, n/2
        write(50,'(2f10.4)') fq(k),spc(k)
      end do 
      close(50)
c
      write(*,*) 'normal terminated'
      pause
      end
c
c-----------------------------------------------------------------------
      subroutine  repwv( eta,n,dt,reph,rept,nowv,fname,ichr )

c --- 時系列データから代表波の波高，周期を求める．

      parameter( nmax = 1000 )
      real  eta(n)
      real  reph(4), rept(4)
      real  wh(nmax), wp(nmax)
      integer  istrt(nmax)

      character*40 fname

      rmwl = 0.0
      do  i = 1, n
         rmwl = rmwl+eta( i )
      end do
      rmwl = rmwl/real( n )

c --- zero up cross method
      ic = 1
      do  i = 1, n-1
         if( ( eta(i).lt.rmwl ) .and. ( eta(i+1).ge.rmwl ) ) then
            istrt( ic ) = i+1
            ic = ic+1
         endif
      end do
      nowv = ic-2  

c --- nowv : number of waves
      if( nowv .lt. 2 ) go to 950

c --- wave height and wave period of each wave
      do  i = 1, nowv
         a = rmwl-eta( istrt(i)-1 )
         b = eta( istrt(i) )-rmwl
         ws = istrt( i )-b/( a+b )
         a = rmwl-eta( istrt(i+1)-1 )
         b = eta( istrt(i+1) )-rmwl
         we = istrt( i+1 )-b/( a+b )
         wp( i ) = ( we-ws )*dt   !  ゼロアップクロス波の周期
         emax = -1.0E6
         emin = 1.0E6
         do  j = istrt(i), istrt(i+1)-1
            if( eta( j ) .ge. emax ) then
               emax = eta( j )
               imax = j
            endif
            if( eta( j ) .le. emin ) then
               emin = eta( j )
               imin = j
            endif
         end do
         a1 = 0.5*( eta(imax-1)-2*eta(imax)+eta(imax+1) )
         b1 = 0.5*( eta(imax+1)-eta(imax-1) )
         if( abs( a1 ) .gt. 1.0e-6 ) emax = emax-b1*b1/( 2*a1 )
         a2 = 0.5*( eta(imin-1)-2*eta(imin)+eta(imin+1) )
         b2 = 0.5*( eta(imin+1)-eta(imin-1) )
         if( abs( a2 ) .gt. 1.0e-6 ) emin = emin-b2*b2/( 2*a2 )
         wh( i ) = emax-emin  !  ゼロアップクロス波の波高(2次近似)
      end do

c ----------------------------------------------------------------------
      open(40,file=fname(1:ichr-4)//'_wav'//'.txt',status='unknown')
      write(40,'(a40)') fname
      write(40,'(a5,2a10)') '   No','     H [m]','     T [s]'
      do i = 1, nowv
        write(40,'(i5,2f10.4)') i,wh(i),wp(i)
      end do 
      close(40)
c ----------------------------------------------------------------------

c --- 波高による並び替え
      do  i = 1, nowv-1
         rmax = wh( i )
         tmax = wp( i )
         imax = i
         do  j = i+1, nowv
            if( wh(j) .gt. rmax ) then
               rmax = wh( j )
               tmax = wp( j )
               imax = j
            endif
         end do
         wh( imax ) = wh( i )
         wh( i ) = rmax
         wp( imax ) = wp( i )
         wp( i ) = tmax
      end do

c --- 各代表波の計算
      reph(1) = wh(1)
      rept(1) = wp(1)

      nwic = nowv/10

      h13 = 0.0
      t13 = 0.0
      do  i = 1, nwic
         h13 = h13+wh( i )
         t13 = t13+wp( i )
      end do
      reph(2) = h13/real( nwic )
      rept(2) = t13/real( nwic )

      nwic = nowv/3

      h13 = 0.0
      t13 = 0.0
      do  i = 1, nwic
         h13 = h13+wh( i )
         t13 = t13+wp( i )
      end do
      reph(3) = h13/real( nwic )
      rept(3) = t13/real( nwic )

      hme = 0.0
      tme = 0.0
      do  i = 1, nowv
         hme = hme+wh( i )
         tme = tme+wp( i )
      end do
      reph(4) = hme/real( nowv )
      rept(4) = tme/real( nowv )
 950  continue

      end

c-----------------------------------------------------------------------
	  subroutine spcwv( eta,n,nb,spc,sigma,dt,df,fq,nf,fname,ichr )
	  PARAMETER( PI = 3.141592E0 )

	  parameter( max = 8192 )	! 作業領域用の配列サイズ
	  real eta( n ), spc( 0:n/2 ), spc0( 0:max/2 )
	  real fq( 0:n/2 )
	  
	  real*8 aa( 0:max )
	  real*8 bb( 0:max )
	  real*8 b2( 0:max )
	  
          character*40 fname
	  
	  if( max < n ) then
		 write( 6,* ) ' Error in spcwv: max < n '
		 stop
	  end if
	
c --- データウィンドー ---
	  ll = n/10
	  do i = 1, n
	    if( (i.ge.1).and.(i.lt.ll) )then
	      b2(i) = 0.5*( 1.-cos( PI*real(i)/real(ll) ) )
	    else if( (i.ge.ll).and.(i.le.n-ll) )then
	      b2(i) = 1.
	    else if( (i.gt.n-ll).and.(i.le.n) )then
	      b2(i) = 0.5*( 1.-cos( PI*(real(n-i))/real(ll) ) )
	    end if
	  end do
	
c --- フーリエ係数の計算 ---
	  do k = 0, n/2
	    sum_a = 0.
	    sum_b = 0.
	    do i = 1, n
	      sum_a = sum_a + eta(i)*cos( PI*2.*real(k)/real(n)*real(i) )
	      sum_b = sum_b + eta(i)*sin( PI*2.*real(k)/real(n)*real(i) )
	    end do
	    aa(k) = 2./real(n)*sum_a
	    if( (k.ge.1).and.(k.le.n/2-1) )then
	      bb(k) = 2./real(n)*sum_b
	    else
	      bb(k) = 0.
	    end if
	  end do
	  
c --- ピリオドグラムの計算 ---
	  uu = 0.
	  do i = 1, n
	    uu = uu + b2(i)**2
	  end do
	  
	  uu = uu/real(n)
	  alpha = real(nb)**2/real(n)/uu
	  
	  do k = 0, n/2
	    spc0(k) = alpha * ( aa(k)**2 + bb(k)**2 )
	  end do
	  	  
c --- ピリオドグラムの平滑化 ---
	  call paraflr( spc0(0),spc(0),n/2+1,nf )
	  
c --- エネルギーレベルの最終調整 ---
	  df = 1./dt/real(n)
	  
	  sumspc = 0.
	  do k = 0, n/2
	    fq(k) = real(k)*df
	    if( k.eq.0 )then
!	      fq(k) = 1./4.*df
	      sumspc = sumspc + spc(k)*df/2.
	    else if( k.eq.n/2 )then
!	      fq(k) = (real(n/2)-1./4.)*df
	      sumspc = sumspc + spc(k)*df/2.
	    else
!	      fq(k) = real(k)*df
	      sumspc = sumspc + spc(k)*df
	    end if
	  end do
	  
	  do k = 0, n/2
	    if( sumspc.gt.0. )then
	      spc(k) = spc(k) * sigma/sumspc
	    else
	      spc(k) = 0.
	    end if
	  end do	  

	  end

	  SUBROUTINE  PARAFLR( SPIN,SPOUT,NDATA,NF )

	  REAL  SPIN( NDATA ),SPOUT( * )
	  if( nf .eq. 0 ) then
		 do i = 1, ndata
			spout( i ) = spin( i )
		 end do
		 return
	  end if
	  RNF = REAL( NF )
	  cof = (1.333333E0*rnf-0.3333333E0/rnf)
	  
	  
	  do i = 1, nf
		 spout( i ) = spin( i )
	  end do

	  DO I = NF+1, NDATA-NF
		 spout( i ) = 0.0
		 DO J = 1, NF-1
			FLTR = (1.0-( J/RNF )**2 ) / cof
			SPOUT( I ) = SPOUT( I )+FLTR*( SPIN(I-J)+SPIN(I+J) )
		 end do
		 spout( i ) = spout( i ) + spin( i ) / cof
	  end do
 	
	  do i = ndata-nf+1, ndata
		 spout( i ) = spin( i )
	  end do

	  end

c-----------------------------------------------------------------------
	  subroutine parspc( spc,n,df,fq,f1,f2,rnp,ep )
	  
	  real spc( 0:n/2 )
	  real fq( 0:n/2 )
	  real rm( 0:4 )
	  
	  do j = 0, 4
		 rm(j) = 0.0
	  end do
	  
	  do i = 0, n/2
	    do j = 0, 4
	      rm(j) = rm(j)+fq(i)**real(j) *spc(i)*df
	    end do
	  end do
	  
	  f1 = rm(1)/rm(0)
	  f2 = sqrt( rm(2)/rm(0) )
	  rnp = sqrt( rm(0)*rm(2)/rm(1)**2-1. )
	  ep  = sqrt( 1.-rm(2)**2/rm(0)/rm(4) )
	
	  end
