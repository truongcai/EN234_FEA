!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D linear elastic continuum element
!
!    The file also contains the following subrouines:
!          abq_UEL_3D_integrationpoints           - defines integration ponits for 3D continuum elements
!          abq_UEL_3D_shapefunctions              - defines shape functions for 3D continuum elements
!          abq_UEL_invert3D                       - computes the inverse and determinant of a 3x3 matrix
!          abq_facenodes_3D                       - returns list of nodes on the face of a 3D element
!
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL_HW7(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
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
      integer      :: i,j,n_points,kint, nfacenodes, ipoin, ie
      integer      :: face_node_list(8)                       ! List of nodes on an element face
    !
      double precision  ::  xi(3,64)                          ! Volumetric Integration points
      double precision  ::  w(64)                             ! Integration weights
      double precision  ::  N(20)                             ! 3D Shape functions
      double precision  ::  dNdxi(20,3)                       ! 3D Shape function derivatives
      double precision  ::  dxdxi(3,3)                        ! Derivative of position wrt normalized coords
      double precision  ::  dNdx(20,3)                        ! Derivative of shape functions wrt reference coords
      double precision  ::  dNdy(20,3)                        ! Derivative of shape functions wrt deformed coords
    !
    !   Variables below are for computing integrals over element faces
      double precision  ::  face_coords(3,8)                  ! Coords of nodes on an element face
      double precision  ::  xi2(2,9)                          ! Area integration points
      double precision  ::  N2(9)                             ! 2D shape functions
      double precision  ::  dNdxi2(9,2)                       ! 2D shape function derivatives
      double precision  ::  norm(3)                           ! Normal to an element face
      double precision  ::  dxdxi2(3,2)                       ! Derivative of spatial coord wrt normalized areal coord
    !
      double precision  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
      double precision  ::  F(3,3)                            ! Deformation gradient
      double precision  ::  Finv(3,3)                         ! Inverse of deformation gradient
      double precision  ::  B(3,3)                            ! C-G deformation tensor
      double precision  ::  C(3,3)                            ! C-G deformation tensor
      double precision  ::  JJ                                ! det(F)
      double precision  ::  H(6,9)
      double precision  ::  D(6,6)                            ! Material tangent
      double precision  ::  Bstar(9,60)                       ! F = Bstar*(dof_total)
      double precision  ::  Smat(3*NNODE,3*NNODE)
      double precision  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
      double precision  ::  qmat(3,3)
      double precision  ::  qvec(9)
      double precision  ::  PK(3,3)
      double precision  ::  kmat(20,20)
      double precision  ::  Cauchy(3,3)

    !
    !     Example ABAQUS UEL implementing 3D linear elastic elements
    !     El props are:

    !     PROPS(1)         Young's modulus
    !     PROPS(2)         Poisson's ratio

      if (NNODE == 4) n_points = 1               ! Linear tet
      if (NNODE == 10) n_points = 4              ! Quadratic tet
      if (NNODE == 8) n_points = 8               ! Linear Hex
      if (NNODE == 20) n_points = 27             ! Quadratic hex

      call abq_UEL_3D_integrationpoints(n_points, NNODE, xi, w)

      if (MLVARX<3*NNODE) then
        write(6,*) ' Error in abaqus UEL '
        write(6,*) ' Variable MLVARX must exceed 3*NNODE'
        write(6,*) ' MLVARX = ',MLVARX,' NNODE = ',NNODE
        stop
      endif

      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0


      ENERGY(1:8) = 0.d0

    !     --  Loop over integration points
      do kint = 1, n_points
        call abq_UEL_3D_shapefunctions(xi(1:3,kint),NNODE,N,dNdxi)
        dxdxi = matmul(coords(1:3,1:NNODE),dNdxi(1:NNODE,1:3))
        call abq_UEL_invert3d(dxdxi,dxidx,determinant)
        dNdx(1:NNODE,1:3) = matmul(dNdxi(1:NNODE,1:3),dxidx)

        ! Assembling Bstar matrix

        Bstar = 0.d0
        Bstar(1,1:3*NNODE-2:3) = dNdx(1:NNODE,1)
        Bstar(2,2:3*NNODE-1:3) = dNdx(1:NNODE,2)
        Bstar(3,3:3*NNODE:3)   = dNdx(1:NNODE,3)
        Bstar(4,1:3*NNODE-2:3) = dNdx(1:NNODE,2)
        Bstar(5,2:3*NNODE-1:3) = dNdx(1:NNODE,1)
        Bstar(6,1:3*NNODE-2:3) = dNdx(1:NNODE,3)
        Bstar(7,3:3*NNODE:3)   = dNdx(1:NNODE,1)
        Bstar(8,2:3*NNODE-1:3) = dNdx(1:NNODE,3)
        Bstar(9,3:3*NNODE:3)   = dNdx(1:NNODE,2)

        ! Caculate the deformation gradient
        do i = 1,3
            ie = 3*(NNODE-1)+i
            F(i,1:3) = matmul(U(i:ie:3),dNdx(1:NNODE,1:3))
            F(i,i) = F(i,i) + 1.d0
        end do
        C = matmul(transpose(F),F)
        ! Calculate F inverse
        call abq_UEL_invert3d(F,Finv,JJ)

        ! Calculate 2nd PK stress and D
        call fung(PROPS(1:NPROPS),NPROPS,F,JJ,PK,D)

        !Assembling q vector
        qmat = matmul(PK,transpose(F))
        qvec = 0.d0
        qvec(1) = qmat(1,1)
        qvec(2) = qmat(2,2)
        qvec(3) = qmat(3,3)
        qvec(4) = qmat(2,1)
        qvec(5) = qmat(1,2)
        qvec(6) = qmat(3,1)
        qvec(7) = qmat(1,3)
        qvec(8) = qmat(3,2)
        qvec(9) = qmat(2,3)


        !Assembling RHS
        RHS(1:3*NNODE,1) = RHS(1:3*NNODE,1)
     1   - matmul(transpose(Bstar(1:9,1:3*NNODE)),qvec)*
     2                                          w(kint)*determinant

        !Assembling H matrix
        H = 0.d0
        H(1,1:9) = [F(1,1),0.d0,0.d0,0.d0,F(2,1),0.d0,F(3,1),0.d0,0.d0]
        H(2,1:9) = [0.d0,F(2,2),0.d0,F(1,2),0.d0,0.d0,0.d0,0.d0,F(3,2)]
        H(3,1:9) = [0.d0,0.d0,F(3,3),0.d0,0.d0,F(1,3),0.d0,F(2,3),0.d0]
        H(4,1:9) = [F(1,2),F(2,1),0.d0,
     1                          F(1,1),F(2,2),0.d0,F(3,2),0.d0,F(3,1)]
        H(5,1:9) = [F(1,3),0.d0,F(3,1),
     1                          0.d0,F(2,3),F(1,1),F(3,3),F(2,1),0.d0]
        H(6,1:9) = [0.d0,F(2,3),F(3,2),
     1                           F(1,3),0.d0,F(1,2),0.d0,F(2,2),F(3,3)]

       AMATRX(1:3*NNODE,1:3*NNODE) = AMATRX(1:3*NNODE,1:3*NNODE)
     1 + matmul( matmul(transpose(Bstar(1:9,1:3*NNODE)),transpose(H)),
     2   matmul(D,matmul(H,Bstar(1:9,1:3*NNODE))) )*w(kint)*determinant

!       Geometric stiffness

        kmat = matmul(matmul(dNdx(1:NNODE,1:3),PK),
     1                             transpose(dNdx(1:NNODE,1:3)))
        Smat = 0.d0
        do i = 0,NNODE-1
          Smat(3*i+1,1:3*NNODE-2:3) = kmat(i+1,1:NNODE)
          Smat(3*i+2,2:3*NNODE-1:3) = kmat(i+1,1:NNODE)
          Smat(3*i+3,3:3*NNODE:3) = kmat(i+1,1:NNODE)
        end do

        AMATRX(1:3*NNODE,1:3*NNODE) = AMATRX(1:3*NNODE,1:3*NNODE) +
     1    Smat*w(kint)*determinant

        !Compute Cauchy stress as (6x1) vector
        Cauchy = matmul(F,matmul(PK,transpose(F)))/JJ
        stress = 0.d0
        stress(1) = Cauchy(1,1)
        stress(2) = Cauchy(2,2)
        stress(3) = Cauchy(3,3)
        stress(4) = Cauchy(1,2)
        stress(5) = Cauchy(1,3)
        stress(6) = Cauchy(2,3)


        if (NSVARS>=n_points*6) then   ! Store Cauchy stress at each integration point (if space was allocated to do so)
            SVARS(6*kint-5:6*kint) = stress(1:6)
        endif
      end do


      PNEWDT = 1.d0          ! This leaves the timestep unchanged (ABAQUS will use its own algorithm to determine DTIME)
    !
    !   Apply distributed loads
    !
    !   Distributed loads are specified in the input file using the Un option in the input file.
    !   n specifies the face number, following the ABAQUS convention.
    !
    !   This is coded to apply nominal tractions to the element face (the residual force does not change as the element deforms)
    !
    !
!      do j = 1,NDLOAD
!
!        call abq_facenodes_3D(NNODE,iabs(JDLTYP(j,1)),
!     1                                     face_node_list,nfacenodes)
!
!        do i = 1,nfacenodes
!            face_coords(1:3,i) = coords(1:3,face_node_list(i))
!        end do
!
!        if (nfacenodes == 3) n_points = 3
!        if (nfacenodes == 6) n_points = 4
!        if (nfacenodes == 4) n_points = 4
!        if (nfacenodes == 8) n_points = 9
!
!        call abq_UEL_2D_integrationpoints(n_points, nfacenodes, xi2, w)
!
!        do kint = 1,n_points
!            call abq_UEL_2D_shapefunctions(xi2(1:2,kint),
!     1                        nfacenodes,N2,dNdxi2)
!            dxdxi2 = matmul(face_coords(1:3,1:nfacenodes),
!     1                           dNdxi2(1:nfacenodes,1:2))
!            norm(1)=(dxdxi2(2,1)*dxdxi2(3,2))-(dxdxi2(2,2)*dxdxi2(3,1))
!            norm(2)=(dxdxi2(1,1)*dxdxi2(3,2))-(dxdxi2(1,2)*dxdxi2(3,1))
!            norm(3)=(dxdxi2(1,1)*dxdxi2(2,2))-(dxdxi2(1,2)*dxdxi2(2,1))
!
!            do i = 1,nfacenodes
!                ipoin = 3*face_node_list(i)-2
!                RHS(ipoin:ipoin+2,1) = RHS(ipoin:ipoin+2,1)
!     1                 - N2(1:nfacenodes)*adlmag(j,1)*norm(1:3)*w(kint)      ! Note determinant is already in normal
!            end do
!        end do
!      end do

      return

      END SUBROUTINE UEL_HW7


      subroutine fung(element_properties,n_properties,F,J,PK,D)

       implicit none

       integer, intent(in)           :: n_properties
       double precision, intent(in)  :: element_properties(n_properties)
       double precision, intent(in)  :: F(3,3)
       double precision, intent(in)  :: J
       double precision, intent(out) :: PK(3,3)
       double precision, intent(out) :: D(6,6)

       double precision :: C(3,3)
       double precision :: Cbar(3,3)
       double precision :: Cinv(3,3)
       double precision :: CQvec(6)
       double precision :: Cvec(6)
       double precision :: Civec(6)
       double precision :: mu
       double precision :: K
       double precision :: ss
       double precision :: G(6,6)
       double precision :: Q
       double precision :: Pdum1(6)
       double precision :: Pdum2(3,3)
       double precision :: P(3,3)
       double precision :: Pvec(6)
       double precision :: Om(6,6)
       double precision :: T1(6,6)
       double precision :: T2(6,6)
       double precision :: T3(6,6)
       double precision :: T4(6,6)

       integer :: i

       !  This subroutine calculates the 2nd PK stress and D
       mu = element_properties(1)
       K  = element_properties(2)
       G = 0.d0
       G(1,1) = element_properties(3)
       G(2,2) = element_properties(4)
       G(3,3) = element_properties(5)
       G(4,4) = element_properties(6)
       G(5,5) = element_properties(6)
       G(6,6) = element_properties(6)

       C = matmul(transpose(F),F)
       Cvec = 0.d0
       Cvec(1) = C(1,1)
       Cvec(2) = C(2,2)
       Cvec(3) = C(3,3)
       Cvec(4) = 2*C(1,2)
       Cvec(5) = 2*C(1,3)
       Cvec(6) = 2*C(2,3)
       call abq_UEL_invert3d(C,Cinv,ss)      ! ss is just a dummy variable here

       ss = J**(-2.d0/3.d0)
       Cbar = C*ss
       CQvec = 0.d0
       CQvec(1) = Cbar(1,1)-1
       CQvec(2) = Cbar(2,2)-1
       CQvec(3) = Cbar(3,3)-1
       CQvec(4) = 2*Cbar(1,2)
       CQvec(5) = 2*Cbar(1,3)
       CQvec(6) = 2*Cbar(2,3)

       Q = 0.25*dot_product(CQvec,matmul(G,CQvec))

       Pdum1 = matmul(G,CQvec)
       Pdum2 = 0.d0
       Pdum2(1,1) = Pdum1(1)
       Pdum2(2,2) = Pdum1(2)
       Pdum2(3,3) = Pdum1(3)
       Pdum2(1,2) = Pdum1(4)
       Pdum2(2,1) = Pdum1(4)
       Pdum2(1,3) = Pdum1(5)
       Pdum2(3,1) = Pdum1(5)
       Pdum2(2,3) = Pdum1(6)
       Pdum2(3,2) = Pdum1(6)

       ! P here is calculated as 3x3 matrix
       P=0.5*ss*(Pdum2-1.d0/3.d0
     1          *dot_product(Cvec,matmul(G,CQvec))*Cinv)


!      Stress calculation PK is 3x3
       PK = mu*exp(Q)*P + K*J*(J-1)*Cinv

       Civec = 0.d0
       Civec(1) = Cinv(1,1)
       Civec(2) = Cinv(2,2)
       Civec(3) = Cinv(3,3)
       Civec(4) = Cinv(1,2)
       Civec(5) = Cinv(1,3)
       Civec(6) = Cinv(2,3)

       Pvec = 0.d0
       Pvec(1) = P(1,1)
       Pvec(2) = P(2,2)
       Pvec(3) = P(3,3)
       Pvec(4) = P(1,2)
       Pvec(5) = P(1,3)
       Pvec(6) = P(2,3)

       !Omega calculation
       Om = 0.d0
       !Fill off-diagonal elements upper triangle
       Om(1,2) = Cinv(1,2)*Cinv(1,2)
       Om(1,3) = Cinv(1,3)*Cinv(1,3)
       Om(1,4) = Cinv(1,1)*Cinv(1,2)
       Om(1,5) = Cinv(1,1)*Cinv(1,3)
       Om(1,6) = Cinv(1,2)*Cinv(1,3)
       Om(2,3) = Cinv(2,3)*Cinv(2,3)
       Om(2,4) = Cinv(2,1)*Cinv(2,2)
       Om(2,5) = Cinv(2,1)*Cinv(2,3)
       Om(2,6) = Cinv(2,2)*Cinv(2,3)
       Om(3,4) = Cinv(3,1)*Cinv(3,2)
       Om(3,5) = Cinv(3,1)*Cinv(3,3)
       Om(3,6) = Cinv(3,2)*Cinv(3,3)
       Om(4,5) = 0.5*(Cinv(1,1)*Cinv(2,3)+Cinv(1,3)*Cinv(1,2))
       Om(4,6) = 0.5*(Cinv(1,2)*Cinv(2,3)+Cinv(1,3)*Cinv(2,2))
       Om(5,6) = 0.5*(Cinv(1,2)*Cinv(3,3)+Cinv(1,3)*Cinv(2,3))

       !Fill all off-diagonal elements
       Om = Om + transpose(Om)

       !Fill diagonal elements
       Om(1,1) = Cinv(1,1)*Cinv(1,1)
       Om(2,2) = Cinv(2,2)*Cinv(2,2)
       Om(3,3) = Cinv(3,3)*Cinv(3,3)
       Om(4,4) =0.5*(Cinv(1,1)*Cinv(2,2)+Cinv(1,2)*Cinv(1,2))
       Om(5,5) =0.5*(Cinv(1,1)*Cinv(3,3)+Cinv(1,3)*Cinv(1,3))
       Om(6,6) =0.5*(Cinv(2,2)*Cinv(3,3)+Cinv(2,3)*Cinv(2,3))

       T1 = G-1.d0/3.d0*( spread(matmul(G,Cvec),dim=2,ncopies=6)*
     1              spread(Civec,dim=1,ncopies=6) +
     2              spread(Civec,dim=2,ncopies=6)*
     3              spread(matmul(G,Cvec),dim=1,ncopies=6) )

       T2 = -1.d0/(ss*3.d0)*(dot_product(Cvec,matmul(G,CQvec))*Om)
     1      +1.d0/9.d0*(dot_product(Cvec,matmul(G,Cvec)))*
     2              spread(Civec,dim=2,ncopies=6)*
     3              spread(Civec,dim=1,ncopies=6)

       T3 = 2*( spread(Pvec,dim=2,ncopies=6)*
     1          spread((Pvec-1/3.d0*Civec),dim=1,ncopies=6) ) -
     2      ss/3.d0*( spread(Civec,dim=2,ncopies=6)*
     3              spread(matmul(G,CQvec),dim=1,ncopies=6) )

       T4 = K*J*( (2*J-1)*spread(Civec,dim=2,ncopies=6)*
     1              spread(Civec,dim=1,ncopies=6) +
     2              2*(J-1)*Om )

       D = mu*exp(Q)*J**(-4.d0/3.d0)*(T1+T2) + mu*exp(Q)*T3 + T4

       return

      end subroutine fung

!      subroutine abq_UEL_invert3d(A,A_inverse,determinant)
!
!      double precision, intent(in) :: A(3,3)
!      double precision, intent(out) :: A_inverse(3,3)
!      double precision, intent(out) :: determinant
!
!      double precision COFACTOR(3,3)
!
!!   Compute inverse and determinant of 3x3 matrix
!
!      determinant =   A(1,1)*A(2,2)*A(3,3)
!     1   - A(1,1)*A(2,3)*A(3,2)
!     2   - A(1,2)*A(2,1)*A(3,3)
!     3   + A(1,2)*A(2,3)*A(3,1)
!     4   + A(1,3)*A(2,1)*A(3,2)
!     5   - A(1,3)*A(2,2)*A(3,1)
!
!      IF (determinant==0.d0) THEN
!        write(6,*) ' Error in subroutine abq_UEL_inver3d'
!        write(6,*) ' A 3x3 matrix has a zero determinant'
!        stop
!      endif
!      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
!      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
!      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
!      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
!      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
!      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
!      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
!      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
!      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))
!
!      A_inverse = transpose(COFACTOR) / determinant
!
!
!      end subroutine abq_UEL_invert3d
