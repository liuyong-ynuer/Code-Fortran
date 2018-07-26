
	implicit none
      integer,Parameter:: nx=288,ny=145,nz1=17,nz2=8,nz3=12,nx_out=128,ny_out=64,nt=12
      Real u(nx,ny,nz1,nt),v(nx,ny,nz1,nt),t(nx,ny,nz1,nt),p(nx,ny,1,nt),q(nx,ny,nz2,nt),z(nx,ny,nz1,nt),rh(nx,ny,nz2,nt),omg(nx,ny,nz3,nt)
      real u_out(nx_out,ny_out,nz1,nt),v_out(nx_out,ny_out,nz1,nt),t_out(nx_out,ny_out,nz1,nt),p_out(nx_out,ny_out,1,nt),q_out(nx_out,ny_out,nz2,nt),z_out(nx_out,ny_out,nz1,nt),rh_out(nx_out,ny_out,nz2,nt),omg_out(nx_out,ny_out,nz3,nt)
      real lon_in(nx),lat_in(ny),lon_out(nx_out),lat_out(ny_out)
      integer i,j,k,l


      open(11,file='E:\study\王旻\LBM\分段\气候平均\Zl17.clim.1989-2009.bin',form='binary')
      open(12,file='E:\study\王旻\LBM\分段\气候平均\RHl8.clim.1989-2009.bin',form='binary')
      open(13,file='E:\study\王旻\LBM\分段\气候平均\Ql8.clim.1989-2009.bin',form='binary')
      open(14,file='E:\study\王旻\LBM\分段\气候平均\Tl17.clim.1989-2009.bin',form='binary')
      open(15,file='E:\study\王旻\LBM\分段\气候平均\Ul17.clim.1989-2009.bin',form='binary')
      open(16,file='E:\study\王旻\LBM\分段\气候平均\Vl17.clim.1989-2009.bin',form='binary')
      open(17,file='E:\study\王旻\LBM\分段\气候平均\Omgl12.clim.1989-2009.bin',form='binary')
      open(18,file='E:\study\王旻\LBM\分段\气候平均\Ps.clim.1989-2009.bin',form='binary')
      open(CONVERT= 'BIG_ENDIAN',unit=51,file='E:\study\王旻\LBM\分段\插值后数据\jra.clim.P2.t42.grd',form='unformatted')
      open(CONVERT= 'BIG_ENDIAN',unit=52,file='E:\study\王旻\LBM\分段\插值后数据\jra.clim.P2.ps.t42.grd',form='unformatted')

         read(11)z
         read(12)rh
         read(13)q
         read(14)t
         read(15)u
         read(16)v
         read(17)omg
         read(18)p
         

         do i=1,nx
            lon_in(i)=(i-1)*1.25
         end do
         do j=1,ny
            lat_in(j)=(j-1)*1.25-90
         end do
         do i=1,nx_out
            lon_out(i)=(i-1)*2.8125
         end do
         data lat_out /-87.864, -85.097, -82.313, -79.526, -76.737, -73.948, -71.158, -68.368, -65.578, -62.787, -59.997, -57.207, -54.416, -51.626, -48.835, -46.045, -43.254, -40.464, -37.673, -34.883, -32.092, -29.301, -26.511, -23.720, -20.930, -18.139, -15.348, -12.558, -9.767, -6.976, -4.186, -1.395, 1.395, 4.186, 6.976, 9.767, 12.558, 15.348, 18.139, 20.930, 23.720, 26.511, 29.301, 32.092, 34.883, 37.673, 40.464, 43.254, 46.045, 48.835, 51.626, 54.416, 57.207, 59.997, 62.787, 65.578, 68.368, 71.158, 73.948, 76.737, 79.526, 82.313, 85.097, 87.864/
       
    do k=1,nt
       call INTP(nx,ny,nz1,nx_out,ny_out,lon_in,lat_in,lon_out,lat_out,-9.9e8,z(:,:,:,k),z_out(:,:,:,k))
       call INTP(nx,ny,nz1,nx_out,ny_out,lon_in,lat_in,lon_out,lat_out,-9.9e8,t(:,:,:,k),t_out(:,:,:,k))
       call INTP(nx,ny,nz1,nx_out,ny_out,lon_in,lat_in,lon_out,lat_out,-9.9e8,u(:,:,:,k),u_out(:,:,:,k))
       call INTP(nx,ny,nz1,nx_out,ny_out,lon_in,lat_in,lon_out,lat_out,-9.9e8,v(:,:,:,k),v_out(:,:,:,k))
       call INTP(nx,ny,nz2,nx_out,ny_out,lon_in,lat_in,lon_out,lat_out,-9.9e8,q(:,:,:,k),q_out(:,:,:,k))
       call INTP(nx,ny,nz2,nx_out,ny_out,lon_in,lat_in,lon_out,lat_out,-9.9e8,rh(:,:,:,k),rh_out(:,:,:,k))
       call INTP(nx,ny,nz3,nx_out,ny_out,lon_in,lat_in,lon_out,lat_out,-9.9e8,omg(:,:,:,k),omg_out(:,:,:,k))
       call INTP(nx,ny,1,nx_out,ny_out,lon_in,lat_in,lon_out,lat_out,-9.9e8,p(:,:,:,k),p_out(:,:,:,k))
    end do

  !输出结果***************************************** 
    do k=1,nt
       do l=1,nz1
          write(51)z_out(:,:,l,k)
       end do
       do l=1,nz2
          write(51)rh_out(:,:,l,k)
       end do
       do l=1,nz2
          write(51)q_out(:,:,l,k)
       end do
       do l=1,nz1
          write(51)t_out(:,:,l,k)
       end do
       do l=1,nz1
          write(51)u_out(:,:,l,k)
       end do
       do l=1,nz1
          write(51)v_out(:,:,l,k)
       end do
       do l=1,nz3
          write(51)omg_out(:,:,l,k)
       end do
       write(52)p_out(:,:,1,k)*0.01
   end do
  ! print*,t_out
  
  end program
  
   
   include "intp.f90"