!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D linear elastic continuum element
!
!    The file also needs the following subrouines:
!          abq_UEL_2D_integrationpoints           - defines integration points for 2D continuum elements
!          abq_UEL_2D_shapefunctions              - defines shape functions for 2D continuum elements
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
    !
      INCLUDE 'ABA_PARAM.INC'
    !
    !
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1   SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2   DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3   JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4   PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

    !
    !       Variables that must be computed in this routine
    !       RHS(i)                     Right hand side vector.  In EN234_FEA the dimensions are always RHS(MLVARX,1)
    !       AMATRX(i,j)                Stiffness matrix d RHS(i)/ d DU(j)
    !       SVARS(1:NSVARS)            Element state variables.  Must be updated in this routine
    !       ENERGY(1:8)
    !                                  Energy(1) Kinetic Energy
    !                                  Energy(2) Elastic Strain Energy
    !                                  Energy(3) Creep Dissipation
    !                                  Energy(4) Plastic Dissipation
    !                                  Energy(5) Viscous Dissipation
    !                                  Energy(6) Artificial strain energy
    !                                  Energy(7) Electrostatic energy
    !                                  Energy(8) Incremental work done by loads applied to the element
    !       PNEWDT                     Allows user to control ABAQUS time increments.
    !                                  If PNEWDT<1 then time step is abandoned and computation is restarted with
    !                                  a time increment equal to PNEWDT*DTIME
    !                                  If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
    !
    !       Variables provided for information
    !       NDOFEL                     Total # DOF for the element
    !       NRHS                       Dimension variable
    !       NSVARS                     Total # element state variables
    !       PROPS(1:NPROPS)            User-specified properties of the element
    !       NPROPS                     No. properties
    !       JPROPS(1:NJPROPS)          Integer valued user specified properties for the element
    !       NJPROPS                    No. integer valued properties
    !       COORDS(i,N)                ith coordinate of Nth node on element
    !       MCRD                       Maximum of (# coords,minimum of (3,#DOF)) on any node
    !       U                          Vector of DOF at the end of the increment
    !       DU                         Vector of DOF increments
    !       V                          Vector of velocities (defined only for implicit dynamics)
    !       A                          Vector of accelerations (defined only for implicit dynamics)
    !       JTYPE                      Integer identifying element type (the number n in the Un specification in the input file)
    !       TIME(1:2)                  TIME(1)   Current value of step time
    !                                  TIME(2)   Total time
    !       DTIME                      Time increment
    !       KSTEP                      Current step number (always 1 in EN234_FEA)
    !       KINC                       Increment number
    !       JELEM                      User assigned element number in ABAQUS (internally assigned in EN234_FEA)
    !       PARAMS(1:3)                Time increment parameters alpha, beta, gamma for implicit dynamics
    !       NDLOAD                     Number of user-defined distributed loads defined for this element
    !       JDLTYP(1:NDLOAD)           Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
    !       ADLMAG(1:NDLOAD)           Distributed load magnitudes
    !       DDLMAG(1:NDLOAD)           Increment in distributed load magnitudes
    !       PREDEF(1:2,1:NPREDF,1:NNODE)   Predefined fields.
    !       PREDEF(1,...)              Value of predefined field
    !       PREDEF(2,...)              Increment in predefined field
    !       PREDEF(1:2,1,k)            Value of temperature/temperature increment at kth node
    !       PREDEF(1:2,2:NPREDF,k)     Value of user defined field/field increment at kth node (not used in EN234FEA)
    !       NPREDF                     Number of predefined fields (1 for en234FEA)
    !       LFLAGS                     Control variable
    !       LFLAGS(1)                  Defines procedure type
    !       LFLAGS(2)                  0 => small displacement analysis  1 => Large displacement (NLGEOM option)
    !       LFLAGS(3)                   1 => Subroutine must return both RHS and AMATRX (always true in EN234FEA)
    !                                   2 => Subroutine must return stiffness AMATRX = -dF/du
    !                                   3 => Subroutine must return daming matrix AMATRX = -dF/dudot
    !                                   4 => Subroutine must return mass matrix AMATRX = -dF/duddot
    !                                   5 => Define the RHS only
    !                                   6 => Define the mass matrix for the initial acceleration calculation
    !                                   100 => Define perturbation quantities for output
    !       LFLAGS(4)                   0 => General step   1 => linear perturbation step
    !       LFLAGS(5)                   0 => current approximation to solution based on Newton correction; 1 => based on extrapolation
    !       MLVARX                      Dimension variable (equal to NDOFEL in EN234FEA)
    !       PERIOD                      Time period of the current step
    !
    !
    ! Local Variables
      integer      :: i,j,n_points,kint
    !
      double precision  ::  xi(2,9)                          ! Area integration points
      double precision  ::  w(9)                             ! Area integration weights
      double precision  ::  N(9)                             ! 2D shape functions
      double precision  ::  dNdxi(9,2)                       ! 2D shape function derivatives
      double precision  ::  dNdx(9,2)                        ! Spatial derivatives
      double precision  ::  Nbar(9)
      double precision  ::  dNbardxi(9,2)
      double precision  ::  dNbardx(9,2)
      double precision  ::  dxdxi(2,2)                       ! Derivative of spatial coords wrt normalized coords

      double precision  ::  sol(9), dsol(9),ep11,ep22,ep12                   ! Sol vector contains [mu, c, dmudx1, dmudx2, dcdx1, dcdx2]
      double precision  ::  q(9)                              ! q vector defined in class
      double precision  ::  D(9,9)                            ! D matrix defined in class
      double precision  ::  B(9,24)             ! p = B*U
      double precision  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
      double precision  ::  diffusion_coeft,kappa,theta,E,xnu,WC,Om       ! Material properties
      double precision  ::  c,d11,d12,dfdc,fc                                 ! concentration
      double precision  ::  st(4)                         ! Stress vector contains [s11, s22, s33, s12]

    !
    !     Example ABAQUS UEL implementing 2D phase field model

      if (NNODE == 3) n_points = 4
      if (NNODE == 4) n_points = 4
      if (NNODE == 6) n_points = 4
      if (NNODE == 8) n_points = 4
      if (NNODE == 9) n_points = 9

      call abq_UEL_2D_integrationpoints(n_points, NNODE, xi, w)

      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0

      E = PROPS(1)
      xnu = PROPS(2)
      Om = PROPS(3)
      WC = PROPS(4)
      kappa = PROPS(5)
      diffusion_coeft = PROPS(6)
      theta = PROPS(7)

      d11 = ( E/( 1.d0+xnu ) )*( 1+ ( xnu/(1.d0-2*xnu) ) )
      d12 = E*xnu/((1.d0+xnu)*(1.d0-2*xnu))
!       write(1,*) U(1:24)
    !     --  Loop over integration points
      do kint = 1, n_points
        call abq_UEL_2D_shapefunctions(xi(1:2,kint),NNODE,N,dNdxi)
        dxdxi = matmul(coords(1:2,1:NNODE),dNdxi(1:NNODE,1:2))
        !invert 2d for dxi/dx
        determinant = dxdxi(1,1)*dxdxi(2,2)-dxdxi(1,2)*dxdxi(2,1)
        dxidx(1,1:2) =  [ dxdxi(2,2),-dxdxi(1,2)]/determinant
        dxidx(2,1:2) =  [-dxdxi(2,1),dxdxi(1,1) ]/determinant
        !dNdx calculation
        dNdx(1:NNODE,1:2) = matmul(dNdxi(1:NNODE,1:2),dxidx)

        call abq_UEL_2D_shapefunctions(xi(1:2,kint),4,Nbar,dNbardxi)
        dxdxi = matmul(coords(1:2,1:4),dNbardxi(1:4,1:2))
        !invert 2d for dxi/dx
        determinant = dxdxi(1,1)*dxdxi(2,2)-dxdxi(1,2)*dxdxi(2,1)
        dxidx(1,1:2) =  [ dxdxi(2,2),-dxdxi(1,2)]/determinant
        dxidx(2,1:2) =  [-dxdxi(2,1),dxdxi(1,1) ]/determinant
        !dNdx calculation
        dNbardx(1:4,1:2) = matmul(dNbardxi(1:4,1:2),dxidx)

        !B matrix calculation
        B = 0.d0
        B(1,1:13:4) = dNdx(1:4,1)
        B(1,17:(2*NNODE+8)-1:2) = dNdx(5:NNODE,1)
        B(2,2:14:4) = dNdx(1:4,2)
        B(2,18:(2*NNODE+8):2) = dNdx(5:NNODE,2)
        B(3,1:13:4) = dNdx(1:4,2)
        B(3,17:(2*NNODE+8)-1:2) = dNdx(5:NNODE,2)
        B(3,2:14:4) = dNdx(1:4,1)
        B(3,18:(2*NNODE+8):2) = dNdx(5:NNODE,1)
        B(4,3:15:4) = Nbar(1:4)
        B(5,4:16:4) = Nbar(1:4)
        B(6,3:15:4) = dNbardx(1:4,1)
        B(7,3:15:4) = dNbardx(1:4,2)
        B(8,4:16:4) = dNbardx(1:4,1)
        B(9,4:16:4) = dNbardx(1:4,2)

        !p vector calculation
        sol = matmul(B(1:9,1:2*NNODE+8),U(1:2*NNODE+8))      ! The p vector at the end of the step
        dsol = matmul(B(1:9,1:2*NNODE+8),DU(1:2*NNODE+8,1))    ! Increment in the p vector
        c = sol(5)

!        write(1,*) 'c',c

        !Stress calculation
        st = 0.d0
        st(1)=E/(1.d0+xnu)*(sol(1)+xnu/(1.d0-2*xnu)*(sol(1)+sol(2)))
     1   - E/(1-2*xnu)*Om/3.d0*c
        st(2)=E/(1.d0+xnu)*(sol(2)+xnu/(1.d0-2*xnu)*(sol(1)+sol(2)))
     1   - E/(1-2*xnu)*Om/3.d0*c
        st(3)=E/(1.d0+xnu)*xnu/(1.d0-2*xnu)*(sol(1)+sol(2))
     1   - E/(1-2*xnu)*Om/3.d0*c
        st(4) = E/(1.d0+xnu)*0.5*sol(3)

        !q vector calculation
!        fc = log(c)-log(1-c) + WC*(1-2*c)
!         fc =WC*(14.346d0*c**5-35.045d0*c**4+29.376d0*c**3-9.693d0*c**2
!     1   + 1.03d0*c-0.007603d0)
!         fc =WC*(11.648d0*c**6-29.706d0*c**5+25.95d0*c**4-8.468d0*c**3
!     1   + 0.4452d0*c**2+0.13462*c-0.000115d0)
!        fc = 2*WC*c*(c-1.d0)*(2.d0*c-1.d0)
!        fc = 8.d0*c**7-21.d0*c**6+19.5d0*c**5-7.5d0*c**4+c**3
        q(1) = st(1)
        q(2) = st(2)
        q(3) = st(4)
        q(4)=
     1  sol(4)-fc-Om/3.d0*(st(1)+st(2)+st(3))
        q(5) = dsol(5)/DTIME
        q(6) = -kappa*sol(8)
        q(7) = -kappa*sol(9)
        q(8) = diffusion_coeft*(sol(6)+(theta-1.d0)*dsol(6))
        q(9) = diffusion_coeft*(sol(7)+(theta-1.d0)*dsol(7))

!        write(1,*) 'q',q(4)

        !Complete D matrix
!        dfdc = 1/(c-c**2) - 2*WC
!      dfdc=71.73d0*c**4-140.18d0*c**3+88.128d0*c**2-19.386d0*c+1.03d0
!        dfdc = 69.888d0*c**5-148.53d0*c**4+103.8d0*c**3-25.404d0*c**2
!     1   + 0.8904d0*c+0.13462
!        dfdc = (12.d0*c**2.d0-12.d0*c+2.d0)
!        dfdc = 56.d0*c**6-126.d0*c**5+97.5d0*c**4-30.d0*c**3+3.d0*c**2
        D = 0.d0
        D(1,1) = d11
        D(1,2) = d12
        D(2,2) = d11
        D(2,1) = d12
        D(3,3) = 0.5*E/(1.d0+xnu)
        D(4,1) = -Om/3.d0*(d11+2.d0*d12)
        D(4,2) = -Om/3.d0*(d11+2.d0*d12)
        D(4,4) = 1.d0
        D(1,5) = -Om/3.d0*(E/(1.d0-2.d0*xnu))
        D(2,5) = -Om/3.d0*(E/(1.d0-2.d0*xnu))
        D(4,5) = -dfdc + Om**2.d0/3.d0*(E/(1.d0-2.d0*xnu))
        D(5,5) = 1.d0/DTIME
        D(6,8) = -kappa
        D(7,9) = -kappa
        D(8,6) = theta*diffusion_coeft
        D(9,7) = theta*diffusion_coeft

        !Calculation of RHS and Stiffness matrix
        RHS(1:2*NNODE+8,1) = RHS(1:2*NNODE+8,1)
     1   - matmul(transpose(B(1:9,1:2*NNODE+8)),q)*w(kint)*determinant

        AMATRX(1:2*NNODE+8,1:2*NNODE+8) =AMATRX(1:2*NNODE+8,1:2*NNODE+8)
     1   + matmul(transpose(B(1:9,1:2*NNODE+8)),
     2             matmul(D,B(1:9,1:2*NNODE+8)))*w(kint)*determinant

        if (NSVARS>=n_points*4) then   ! Store stress at each integration point (if space was allocated to do so)
            SVARS(4*kint-3:4*kint) = st(1:4)
        end if
      end do


      return

      END SUBROUTINE UEL
