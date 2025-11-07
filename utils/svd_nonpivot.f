ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c 
c 
c       This is the end of the debugging code, and the beginning of the
c       SVD code proper
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains two user-callable subroutines: svd_nonpivot
c       and svd_nonpivot2. Below is a short description of these two
c       subroutines.
c 
c  svd_nonpivot - produces the Singular Value Decomposition (SVD) 
c       of the user-supplied matrix a, representing A in the form
c 
c                 a=u(n,n) * D(n,n) * (v(m,n))^*                    (1)
c
c       where columns columns of the matrices u,v are orthonormal, 
c       and the matrrix D is diagonal. A is expected to be (more 
c       or less) full rank, so the subroutine does NOT use any 
c       preliminary compression.
c
c  svd_nonpivot2 - calculates the Singular Values (but no singular
c       vectors of the user-supplied matrix a.
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c 
c 
c 
        subroutine svd_nonpivot2(a,n,m,rlams,w)
        implicit real *8 (a-h,o-z)
        dimension a(n,m),rlams(1),w(1)
c
c       This subroutine calculates the Singular Values of the 
c       user-supplied matrix a.
c
c                     Input parameters:
c
c  a - the n \times m matrix to be decomposed. Please note that it is NOT
c       damaged by the subroutine in any way
c  n,m - the dimensionalities of the matrix
c  
c                     Output parameters:
c
c  rlams - singular values of a
c  
c                     Work arrays:
c
c  w - must be at least m*n+20n+2*m+1000
c
c      . . . if n > m, transpose; otherwise, decompose and get out
c
c
        if(n .gt. m) goto 2000
c
        call svd_nonpivot02(a,n,m,rlams,w)
        return
c
 2000 continue
c
c       n is greater than m. svd the transposed matrix
c
        call svd_nonpivot_transp(a,n,m,w)
c
c       . . . now, a is an m \times n matrix; decompose
c
        call svd_nonpivot02(a,m,n,rlams,w)
c
c       transpose a back
c
        call svd_nonpivot_transp(a,m,n,w)
c
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_nonpivot(a,n,m,u,v,rlams,w)
        implicit real *8 (a-h,o-z)
        dimension a(n,m),v(m,m),u(n,n),rlams(1),w(1)
c
c       This subroutine produces the Singular Value Decomposition
c       (SVD) of the user-supplied matrix a, representing A in the 
c       form
c 
c                 a=u(n,n) * D(n,n) * (v(m,n))^*                    (1)
c
c       where columns columns of the matrices u,v are orthonormal,
c       and the matrrix D is diagonal. A is expected to be (more 
c       or less) full rank, so the subroutine does NOT use any 
c       preliminary compression.
c
c                     Input parameters:
c
c  a - the n \times m matrix to be decomposed. Please note that it is NOT
c       damaged by the subroutine in any way
c  n,m - the dimensionalities of the matrix
c  
c                     Output parameters:
c
c  u,v - matrices in (1)
c  rlams - diagonal elements of the matrix D in (1); will be 
c       non-negative and in non-increasing order
c  
c                     Work arrays:
c
c  w - must be at least 5*m*m+12n+1000
c
c      . . . if n > m, transpose; otherwise, decompose and get out
c
c
        if(n .gt. m) goto 2000
c
        call svd_nonpivot0(a,n,m,u,v,rlams,w)
        return
c
 2000 continue
c
c       n is greater than m. svd the transposed matrix
c
        call svd_nonpivot_transp(a,n,m,w)
c
c       . . . now, a is an m \times n matrix; decompose
c
        call svd_nonpivot0(a,m,n,v,u,rlams,w)
c
c       transpose a back
c
        call svd_nonpivot_transp(a,m,n,w)
c
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_nonpivot02(a,n,m,rlams,w)
        implicit real *8 (a-h,o-z)
        dimension a(n,m),rlams(1),w(1)
c
c       allocate memory
c
        ida=1
        lda=n+2
c
        isuba=ida+lda
        lsuba=n+2
c
c       attention!
c
        iuut=1
        luut=m*2+10
c
        iw=iuut+luut 
c
        lw=n*m+18*n+100
c
        ltot=lw+iw
        call prinf('and ltot=*',ltot,1)
c
        call svd_nonpivot002(a,n,m,rlams,
     1      w(ida),w(isuba),w(iw),w(iuut))
c
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_nonpivot002(a,n,m,rlams,
     1      da,suba,v_as,uut)
        implicit real *8 (a-h,o-z)
        dimension a(n,m),v_as(n,m),da(1),suba(1),
     1      rlams(1),uut(2,1)
c
c       This subroutine constructs a partial SVD of the
c       user-supplied n \times m matrix a. 
c
        call svd_nonpivot_copy(a,v_as,n*m)
c
c       bidiagonalize the matrix a
c
        call svd_nonpivot_bidiag2(v_as,n,m,uut)
c
c       SVD the obtained bidiagonal matrix
c
        do 2200 i=1,n-1
c
        da(i)=v_as(i,i)
        suba(i)=v_as(i+1,i)
 2200 continue
c
        da(n)=v_as(n,n)
c
        call qleigen_bidiag_vals(ier,n,da,suba,
     1      rlams,v_as)
c
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_nonpivot_bidiag2(a,n,m,uut)
        implicit real *8 (a-h,o-z)
        dimension a(n,m),bb(3),ub(2,2),uut(2,1)
c
c       This subroutine bidiagonalizes the the user-supplied 
c       matrix a(n,m). It also returns the orthogonal matrices 
c       v(m,m), u(n,n), such that 
c 
c                 b=u(n,n) * a(n,m) * v(m,m),
c
c       where a is the original matrix, and b is the didiagonal
c       one. Please note that the b is LOWER bidiagonal, i.e. 
c       it consists of a diagonal and a SUBdiagonal.
c
c                 Input paramaters:
c
c  a - the n \times m - matrix to be partially diagonalized.
c       Please note that on output, it is replaced with the 
c       partially bodiagonalized matrix

c
c        set v, ul to identity
c
        ich=0
        do 2600 ii=1,n
c
        ich=ich+1
        if(ich .eq. 10) then
            call prinf('ii=*',ii,1)
            ich=0
        endif
c
c       reduce the ii-th row to its ii-th element
c
        do 2200 i=ii+1,m
c
        bb(1)=a(ii,ii)
        bb(2)=a(ii,i)
        call svd_nonpivot_rotfnd(bb,ub)
c
        i0=ii
        call svd_nonpivot_cols_rotate(ub,a,n,ii,i,i0)
c
 2200 continue
c
c       reduce the ii-th column to its ii+1-st element
c
        icol=ii
        jj0=ii+2
        njj=n-1-ii
c
        if(njj .le. 0) goto 2600
c
        call svd_nonpivot_rows_block2(a,n,m,icol,ii,jj0,njj,uut)
c
 2600 continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_nonpivot_rows_block2(a,n,m,icol,ii,jj0,njj,uut)
        implicit real *8 (a-h,o-z)
        dimension a(n,n),bb(3),ub(2,2),uut(2,1)
c 
c       construct the vector of values to be eliminated
c
        val0=a(ii+1,icol)
c
c       construct the array of eliminating matrices
c
        do 1600 i=1,njj
c
        val1=a(jj0+i-1,icol)


        bb(1)=val0
        bb(2)=val1
        call svd_nonpivot_rotfnd(bb,ub)
c
        d1=ub(1,1)*val0+ub(1,2)*val1
        d2=ub(2,1)*val0+ub(2,2)*val1
c
        val0=d1
c
        uut(1,i)=ub(1,1)
        uut(2,i)=ub(2,1)
c
 1600 continue
c
c       now, apply the obtained matrices to the rows of the matrix a
c
        do 3000 j=ii,m,4
c
        if(j .lt. m-5) goto 2000

        do 1800 i=1,njj
c
        jj=jj0+i-1
c
        d1=uut(1,i)*a(ii+1,j)+uut(2,i)*a(jj,j)
        a(jj,j)=uut(2,i)*a(ii+1,j)-uut(1,i)*a(jj,j)
        a(ii+1,j)=d1
c
        if(j+1 .gt. m) goto 1650
        d1=uut(1,i)*a(ii+1,j+1)+uut(2,i)*a(jj,j+1)
        a(jj,j+1)=uut(2,i)*a(ii+1,j+1)-uut(1,i)*a(jj,j+1)
        a(ii+1,j+1)=d1
 1650 continue
c
        if(j+2 .gt. m) goto 1700
c
        d1=uut(1,i)*a(ii+1,j+2)+uut(2,i)*a(jj,j+2)
        a(jj,j+2)=uut(2,i)*a(ii+1,j+2)-uut(1,i)*a(jj,j+2)
        a(ii+1,j+2)=d1
 1700 continue
c
        if(j+3 .gt. m) goto 1800
c
        d1=uut(1,i)*a(ii+1,j+3)+uut(2,i)*a(jj,j+3)
        a(jj,j+3)=uut(2,i)*a(ii+1,j+3)-uut(1,i)*a(jj,j+3)
        a(ii+1,j+3)=d1
 1800 continue
c
        goto 3000
c
 2000 continue
c
c       we are far away from the right border of a. Drop 
c       all precautions
c

        do 2600 i=1,njj
c
        jj=jj0+i-1
c
        d1=uut(1,i)*a(ii+1,j)+uut(2,i)*a(jj,j)
        a(jj,j)=uut(2,i)*a(ii+1,j)-uut(1,i)*a(jj,j)
        a(ii+1,j)=d1
c
        d1=uut(1,i)*a(ii+1,j+1)+uut(2,i)*a(jj,j+1)
        a(jj,j+1)=uut(2,i)*a(ii+1,j+1)-uut(1,i)*a(jj,j+1)
        a(ii+1,j+1)=d1
c
        d1=uut(1,i)*a(ii+1,j+2)+uut(2,i)*a(jj,j+2)
        a(jj,j+2)=uut(2,i)*a(ii+1,j+2)-uut(1,i)*a(jj,j+2)
        a(ii+1,j+2)=d1
c
        d1=uut(1,i)*a(ii+1,j+3)+uut(2,i)*a(jj,j+3)
        a(jj,j+3)=uut(2,i)*a(ii+1,j+3)-uut(1,i)*a(jj,j+3)
        a(ii+1,j+3)=d1
c
 2600 continue
c
 3000 continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_nonpivot0(a,n,m,u,v,rlams,w)
        implicit real *8 (a-h,o-z)
        dimension a(n,m),v(m,m),u(n,n),rlams(1),w(1)
c
c       allocate memory
c
        iuus=1
        luus=(n**2-3*n+2)+10
c
        ivvs=iuus+luus
        lvvs=(2*n*m-n**2-n)+10
c
        i_iiu=ivvs+lvvs
        l_iiu=(n**2-3*n+2)+10
c
        i_iiv=i_iiu+l_iiu
        l_iiv=(2*n*m-n**2-n)+10

        ida=i_iiv+l_iiv
        lda=n+2
c
        isuba=ida+lda
        lsuba=n+2
c
        iw=isuba+lsuba
        lw=n*m+18*n+100
c
        ltot=lw+iw
        call prinf('and ltot=*',ltot,1)
c
        call svd_nonpivot00(a,n,m,u,v,rlams,
     1      w(ida),w(isuba),w(iw),v,
     2      w(iuus),w(ivvs),w(i_iiu),w(i_iiv) )
c
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_nonpivot00(a,n,m,u,v,rlams,
     1      da,suba,w,v_as,
     2      uus,vvs,iiu,iiv)
        implicit real *8 (a-h,o-z)
        dimension a(n,m),v(m,m),u(n,n),v_as(n,m),da(1),suba(1),
     1      rlams(1),w(1),uus(2,1),iiu(2,1),vvs(2,1),iiv(2,1)
c
c       This subroutine constructs a partial SVD of the
c       user-supplied n \times m matrix a. 
c
        call svd_nonpivot_copy(a,v_as,n*m)
c
c       bidiagonalize the matrix a
c
        call svd_nonpivot_bidiag(v_as,n,m,uus,iiu,niu,
     1      vvs,iiv,niv,w)
c
c       SVD the obtained bidiagonal matrix
c
        do 2200 i=1,n-1
c
        da(i)=v_as(i,i)
        suba(i)=v_as(i+1,i)
 2200 continue
c
        da(n)=v_as(n,n)
c
        call qleigen_bidiag_vects(ier,n,da,suba,
     1      rlams,u,v,w)
c
c       multiply the matrices u,v bidiagonalizing the matrix a
c       by the matrices uk,vk diagonalizing the bidiagonal 
c       matrix
c

        call svd_nonpivot_combine(n,m,u,v,rlams,w,
     1      uus,iiu,niu,vvs,iiv,niv)
c
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_nonpivot_combine(n,m,u,v,rlams,w,
     1      uus,ius,niu,vvs,ivs,niv)
        implicit real *8 (a-h,o-z)
        dimension v(m,m),rlams(1),u(n,n),w(n,n),uus(2,1),
     1      ius(2,1),vvs(2,1),ivs(2,1)
c
c       multiply u^* by uk
c
        call svd_nonpivot_transp(u,n,n,w)
        call svd_nonpivot_bidiag_apply22(n,u,uus,ius,niu)
        call svd_nonpivot_transp(u,n,n,w)
c
c       multiply vk^* by v^*
c
        call svd_nonpivot_copy(v,w,m*n)
        call svd_nonpivot_setzero(v,m*n)
        do 1600 i=1,n
        do 1400 j=1,n
c
        v(j,i)=w(j,i)
 1400 continue
 1600 continue
c
        call svd_nonpivot_transp(v,m,n,w)
        call svd_nonpivot_bidiag_apply33(n,m,v,vvs,ivs,niv)
        call svd_nonpivot_transp(v,n,m,w)
c
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_nonpivot_bidiag_apply33(n,m,v,vvs,ivs,niv)
        implicit real *8 (a-h,o-z)
        dimension ub(2,2),v(m,m),vvs(2,1),ivs(2,1)
c
        do 2400 ii=niv,1,-1
c
        ub(1,1)=vvs(1,ii)
        ub(2,1)=vvs(2,ii)
        ub(2,2)=-ub(1,1)
        ub(1,2)=ub(2,1)
c
        j1=ivs(1,ii)
        j2=ivs(2,ii)
c
        i0=1
        call svd_nonpivot_cols_rotate(ub,v,n,j1,j2,i0)
 2400 continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_nonpivot_bidiag_apply22(n,u,uus,ius,niu)
        implicit real *8 (a-h,o-z)
        save
        dimension ub(2,2),u(n,n),uus(2,1),ius(2,1)
c
        do 2400 ii=niu,1,-1
c
        ub(1,1)=uus(1,ii)
        ub(2,1)=uus(2,ii)
        ub(2,2)=-ub(1,1)
        ub(1,2)=ub(2,1)
c
        j1=ius(1,ii)
        j2=ius(2,ii)
c
        i0=1
        call svd_nonpivot_cols_rotate(ub,u,n,j1,j2,i0)
 2400 continue
c
        return
        end

c 
c 
c 
c 
c 
        subroutine svd_nonpivot_bidiag(a,n,m,uus,ius,niu,
     1      vvs,ivs,niv,uut)
        implicit real *8 (a-h,o-z)
        dimension a(n,m),bb(3),ub(2,2),uus(2,1),
     1      ius(2,1),vvs(2,1),ivs(2,1),uut(2,1)
c
c       This subroutine cbidiagonalizes the the user-supplied 
c       matrix a(n,m). It also returns the orthogonal matrices 
c       v(m,m), u(n,n), such that 
c 
c                 b=u(n,n) * a(n,m) * v(m,m),
c
c       where a is the original matrix, and b is the didiagonal
c       one. Please note that the b is LOWER bidiagonal, i.e. 
c       it consists of a diagonal and a SUBdiagonal.
c
c                 Input paramaters:
c
c  a - the n \times m - matrix to be partially diagonalized.
c       Please note that on output, it is replaced with the 
c       partially bodiagonalized matrix

c
c        set v, ul to identity
c
        iiu=0
        iiv=0
c
        ich=0
        do 2600 ii=1,n
c
        ich=ich+1
        if(ich .eq. 10) then
            call prinf('ii=*',ii,1)
            ich=0
        endif
c
c       reduce the ii-th row to its ii-th element
c
        do 2200 i=ii+1,m
c
        bb(1)=a(ii,ii)
        bb(2)=a(ii,i)
        call svd_nonpivot_rotfnd(bb,ub)
c
        i0=ii
        call svd_nonpivot_cols_rotate(ub,a,n,ii,i,i0)
c
        iiv=iiv+1
        vvs(1,iiv)=ub(1,1)
        vvs(2,iiv)=ub(2,1)
c
        ivs(1,iiv)=ii
        ivs(2,iiv)=i
c
 2200 continue
c
c       reduce the ii-th column to its ii+1-st element
c
        icol=ii
        jj0=ii+2
        njj=n-1-ii
c
        if(njj .le. 0) goto 2600
c
        call svd_nonpivot_rows_block(a,n,m,icol,ii,jj0,njj,
     1      uus,ius,iiu,uut)
c
 2600 continue
c
        niu=iiu
        niv=iiv
c
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_nonpivot_rows_block(a,n,m,icol,ii,jj0,njj,
     1      uus,ius,iiu,uut)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),bb(3),ub(2,2),uus(2,1),ius(2,1),uut(2,1)
c 
c       construct the vector of values to be eliminated
c
        val0=a(ii+1,icol)
c
c       construct the array of eliminating matrices
c
        do 1600 i=1,njj
c
        val1=a(jj0+i-1,icol)


        bb(1)=val0
        bb(2)=val1
        call svd_nonpivot_rotfnd(bb,ub)
c
        d1=ub(1,1)*val0+ub(1,2)*val1
        d2=ub(2,1)*val0+ub(2,2)*val1
c
        val0=d1
c
        iiu=iiu+1
        uus(1,iiu)=ub(1,1)
        uus(2,iiu)=ub(2,1)
c
        ius(1,iiu)=ii+1
        ius(2,iiu)=jj0+i-1
c
        uut(1,i)=ub(1,1)
        uut(2,i)=ub(2,1)
c
 1600 continue
c
c       now, apply the obtained matrices to the rows of the matrix a
c
        do 3000 j=ii,m,4
c
        if(j .lt. m-5) goto 2000

        do 1800 i=1,njj
c
        jj=jj0+i-1
c
        d1=uut(1,i)*a(ii+1,j)+uut(2,i)*a(jj,j)
        a(jj,j)=uut(2,i)*a(ii+1,j)-uut(1,i)*a(jj,j)
        a(ii+1,j)=d1
c
        if(j+1 .gt. m) goto 1650
        d1=uut(1,i)*a(ii+1,j+1)+uut(2,i)*a(jj,j+1)
        a(jj,j+1)=uut(2,i)*a(ii+1,j+1)-uut(1,i)*a(jj,j+1)
        a(ii+1,j+1)=d1
 1650 continue
c
        if(j+2 .gt. m) goto 1700
c
        d1=uut(1,i)*a(ii+1,j+2)+uut(2,i)*a(jj,j+2)
        a(jj,j+2)=uut(2,i)*a(ii+1,j+2)-uut(1,i)*a(jj,j+2)
        a(ii+1,j+2)=d1
 1700 continue
c
        if(j+3 .gt. m) goto 1800
c
        d1=uut(1,i)*a(ii+1,j+3)+uut(2,i)*a(jj,j+3)
        a(jj,j+3)=uut(2,i)*a(ii+1,j+3)-uut(1,i)*a(jj,j+3)
        a(ii+1,j+3)=d1
 1800 continue
c
        goto 3000
c
 2000 continue
c
c       we are far away from the right border of a. Drop 
c       all precautions
c

        do 2600 i=1,njj
c
        jj=jj0+i-1
c
        d1=uut(1,i)*a(ii+1,j)+uut(2,i)*a(jj,j)
        a(jj,j)=uut(2,i)*a(ii+1,j)-uut(1,i)*a(jj,j)
        a(ii+1,j)=d1
c
        d1=uut(1,i)*a(ii+1,j+1)+uut(2,i)*a(jj,j+1)
        a(jj,j+1)=uut(2,i)*a(ii+1,j+1)-uut(1,i)*a(jj,j+1)
        a(ii+1,j+1)=d1
c
        d1=uut(1,i)*a(ii+1,j+2)+uut(2,i)*a(jj,j+2)
        a(jj,j+2)=uut(2,i)*a(ii+1,j+2)-uut(1,i)*a(jj,j+2)
        a(ii+1,j+2)=d1
c
        d1=uut(1,i)*a(ii+1,j+3)+uut(2,i)*a(jj,j+3)
        a(jj,j+3)=uut(2,i)*a(ii+1,j+3)-uut(1,i)*a(jj,j+3)
        a(ii+1,j+3)=d1
c
 2600 continue
c
 3000 continue
c
        return
        end
c 
c 
c 
c
c 
        subroutine svd_nonpivot_cols_rotate(u,a,n,ii,jj,i0)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),u(2,2)
c 
        do 1200 i=i0,n
c 
        d1=u(1,1)*a(i,ii)+u(1,2)*a(i,jj)
        d2=u(2,1)*a(i,ii)+u(2,2)*a(i,jj)
c 
        a(i,ii)=d1
        a(i,jj)=d2
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_nonpivot_copy(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        real *8 x(1),y(1)
c 
        do 2400 i=1,n
        y(i)=x(i)
 2400 continue
        return
        end  
c 
c 
c 
c 
c 
        subroutine svd_nonpivot_setzero(x,n)
        implicit real *8 (a-h,o-z)
        save
        dimension x(n)
c 
        do 2000 i=1,n
        x(i)=0
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_nonpivot_transp(a,n,m,w)
        implicit real *8 (a-h,o-z)
        dimension a(m,n),w(n,m)
c
        call svd_nonpivot_copy(a,w,n*m)
        do 1400 i=1,n
        do 1200 j=1,m
c
        a(j,i)=w(i,j)
 1200 continue
 1400 continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_nonpivot_rotfnd(a,u)
        implicit real *8 (a-h,o-z)
        save
        dimension a(2),u(2,2)
c 
        u21=-a(2)
        u22=a(1)
c 
        d=sqrt(u22*u22+u21*u21)
c 
        if(d .eq. 0) then
c 
            u(2,2)=1
            u(1,2)=0
            u(1,1)=1
            u(2,1)=0
            return
        endif
c 
        u(2,2)=u22/d
        u(2,1)=u21/d
c 
        u(1,1)=-u(2,2)
        u(1,2)=u(2,1)
        return
        end
