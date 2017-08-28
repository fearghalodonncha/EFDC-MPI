SUBROUTINE MHKPWRDIS  
!**********************************************************************C  
! **  SUBROUTINE MHKPWRDIS CALCULATES POWER DISSIPATION FROM MARINE HYDROKINETIC  
! **  DEVICES AS A FUNCTION OF WATER VELOCITY  
!
!     04/2010  Bill Arnold and Scott James
!              Sandia National Laboratories
!  
!**********************************************************************C  
  USE GLOBAL
  IMPLICIT NONE
  REAL::FLOWSPEED=0.0,THRSTCOEF=0.0,ZTOP=0.0,ZBOTTOM=0.0,MAXFLOWSPD=-999.
  REAL::FS_WF,FS_NF,FS_EF,FS_SF,VELUP,UVELUP,VVELUP,FMHK,FSUP,USPD,VSPD,MAXSPD
  REAL::SUMLAYM,SUMLAYS,SUMNEGLAYM,SUMNEGLAYS
  REAL::LAYFRACM(KC),LAYFRACS(KC),NEGLAYFRACM(KC),NEGLAYFRACS(KC)
  REAL::UATVFACE,VATUFACE,UATVFACEN,VATUFACEE,SPDN,SPDS,SPDE,SPDW,AWEIGHTX1,AWEIGHTX2,AWEIGHTY1,AWEIGHTY2
  INTEGER::MHKCOUNT,M,L,K,LW,LE,LS,LN,LNW,LSE,LNE,FACE
!  INTEGER::LAYCOUNTM(TCOUNT),LAYCOUNTS(TCOUNT),LAYCOUNTLASTM(TCOUNT),LAYCOUNTLASTS(TCOUNT)
  LOGICAL::STATUS
  LOGICAL,SAVE::FIRSTTIME=.FALSE.
!**********************************************************************C 
  INQUIRE(FILE='POWERBUG.DAT',EXIST=STATUS)
  IF(STATUS.AND.DEBUG.AND..NOT.FIRSTTIME)THEN
    OPEN(357,FILE='POWERBUG.DAT')
    FIRSTTIME=.TRUE.
  ENDIF
  MHKCOUNT=0 !initialize the running count of the number of cells with MHKdevices (will equal TCOUNT)
  FXMHKE(:)=0.0;FYMHKE(:)=0.0;FXSUPE(:)=0.0;FYSUPE(:)=0.0 !initialize x,y-forces from MHK device/support for the external solution
  FXMHK(:,:)=0.0;FYMHK(:,:)=0.0;FXSUP(:,:)=0.0;FYSUP(:,:)=0.0 !initialize x,y-forces from the MHK device/support for the internal solution
  PMHK(:,:)=0.0;PSUP(:,:)=0.0 !initialize power extracted from the MHK device/support from flow
  DO L=2,LA  !loop through model area
    IF(.NOT.LMASKDRY(L).OR.MVEGL(L)<90)CYCLE !if the cell is dry or does not contain an MHK, skip this cell
    MHKCOUNT=MHKCOUNT+1
    LW=LWEST(L)      !west cell, I-1,J
    LE=LEAST(L)      !east cell, I+1,J
    LS=LSC(L)   !south cell, I,J-1
    LN=LNC(L)   !north cell, I,J+1
    LNW=LNWC(L) !northwest cell, I-1,J+1
    LNE=LNEC(L) !northeast cell, I+1,J+1
    LSE=LSEC(L) !southeast cell, I+1,J-1
    M=MVEGL(L)-90 !This was put in to have MHKs be vegetative inputs > 90; this is the mhktype
    LAYFRACM(:)=0.0;NEGLAYFRACM(:)=0.0 !initialize the layer fraction variable MHK
    LAYFRACS(:)=0.0;NEGLAYFRACS(:)=0.0 !initialize the layer fraction variable support
    MAXFLOWSPD=-999. !initialize maximum flow speed for error checking
    DO K=1,KC !MHK device layer filler - which layers does the device occupy and at what fraction
      ZTOP=HP(L)*Z(K)+BELV(L) !layer top elevation
      IF(ZTOP<ZMINMHK(M,L))CYCLE !layer is below device
      ZBOTTOM=HP(L)*Z(K-1)+BELV(L) !layer bottom elevation
      IF(ZBOTTOM>ZMAXMHK(M,L))CYCLE !layer is above device
      IF(ZTOP>=ZMAXMHK(M,L).AND.ZBOTTOM<ZMINMHK(M,L))THEN !device is wholly contained in this layer (special case)
        LAYFRACM(K)=(ZMAXMHK(M,L)-ZMINMHK(M,L))/(HP(L)*DZC(K)) !calculate fraction of layer that is occupied
        EXIT
      ENDIF
      IF(ZMAXMHK(M,L)>=ZTOP.AND.ZMINMHK(M,L)<ZBOTTOM)THEN !this layer is fully occupied by the device
        LAYFRACM(K)=1.0
        CYCLE
      ENDIF
      IF(ZBOTTOM<ZMINMHK(M,L).AND.ZMAXMHK(M,L)>=ZTOP)THEN !this layer is partially occupied by the device (bottom)
        LAYFRACM(K)=(ZTOP-ZMINMHK(M,L))/(HP(L)*DZC(K)) !calculate the fraction of layer that is occupied
        CYCLE
      ENDIF
      IF(ZTOP>=ZMAXMHK(M,L).AND.ZMINMHK(M,L)<ZBOTTOM)THEN !this layer is partially occupied by the device (top)
        LAYFRACM(K)=(ZMAXMHK(M,L)-ZBOTTOM)/(HP(L)*DZC(K)) !calculate the fraction of layer that is occupied
        CYCLE
      ENDIF
    ENDDO
    NEGLAYFRACM(:)=LAYFRACM(:)-1.0
    SUMLAYM=SUM(LAYFRACM(1:KC));SUMNEGLAYM=SUM(NEGLAYFRACM(1:KC)) !Sum of MHK layer fractions
    DO K=1,KC !MHK support layer filler - which layers does the support occupy and at what fraction
      ZTOP=HP(L)*Z(K)+BELV(L) !layer top elevation
      IF(ZTOP<ZMINSUP(M,L))CYCLE !layer is below support
      ZBOTTOM=HP(L)*Z(K-1)+BELV(L) !layer bottom elevation
      IF(ZBOTTOM>ZMAXSUP(M,L))CYCLE !layer is above support
      IF(ZTOP>=ZMAXSUP(M,L).AND.ZBOTTOM<ZMINSUP(M,L))THEN !support is wholly contained in this layer (special case)
        LAYFRACS(K)=(ZMAXSUP(M,L)-ZMINSUP(M,L))/(HP(L)*DZC(K)) !calculate fraction of layer that is occupied
        EXIT
      ENDIF
      IF(ZMAXSUP(M,L)>=ZTOP.AND.ZMINSUP(M,L)<ZBOTTOM)THEN !this layer is fully occupied by the support
        LAYFRACS(K)=1.0
        CYCLE
      ENDIF
      IF(ZBOTTOM<ZMINSUP(M,L).AND.ZMAXSUP(M,L)>=ZTOP)THEN !this layer is partially occupied by the support (bottom)
        LAYFRACS(K)=(ZTOP-ZMINSUP(M,L))/(HP(L)*DZC(K)) !calculate the fraction of layer that is occupied
        CYCLE
      ENDIF
      IF(ZTOP>=ZMAXSUP(M,L).AND.ZMINSUP(M,L)<ZBOTTOM)THEN !this layer is partially occupied by the support (top)
        LAYFRACS(K)=(ZMAXSUP(M,L)-ZBOTTOM)/(HP(L)*DZC(K)) !calculate the fraction of layer that is occupied
        CYCLE
      ENDIF
    ENDDO
    NEGLAYFRACS(:)=LAYFRACS(:)-1.0  
    SUMLAYS=SUM(LAYFRACS(1:KC));SUMNEGLAYS=SUM(NEGLAYFRACS(1:KC)) !Sum of support layer fractions
    DO K=1,KC
      ZTOP=HP(L)*Z(K)+BELV(L) !layer top elevation
      ZBOTTOM=HP(L)*Z(K-1)+BELV(L) !layer bottom elevation
      IF(ZMAXMHK(M,L)>HP(L).OR.ZMAXSUP(M,L)>HP(L))THEN !Protruding from water?
        PRINT*,'MHK DEVICE PROTRUDING FROM WATER',HP(L),ZMAXMHK(M,L),L,IL(L),JL(L)
        STOP
      ELSEIF(ZMINMHK(M,L)<BELV(L))THEN !Below bed elevation?
        PRINT*,'MKH DEVICE IN SEDIMENT',BELV(L),ZMINMHK(M,L),L,IL(L),JL(L)
        STOP
      ENDIF
      IF((ZBOTTOM>ZMAXMHK(M,L).AND.ZBOTTOM>ZMAXSUP(M,L)).OR.(ZTOP<ZMINMHK(M,L).AND. &
ZTOP<ZMINSUP(M,L)).OR.(U(L,K)==0.0.AND.V(L,K)==0.0))CYCLE !no MHK or support or velocity in this layer
      USPD=ABS(0.5*(U(L,K)+U(LE,K)))      !I,J cell center u-speed
      VSPD=ABS(0.5*(V(L,K)+V(LN,K)))      !I,J cell center v-speed
      FLOWSPEED=SQRT(USPD*USPD+VSPD*VSPD) !I,J cell center speed
      MAXFLOWSPD=MAX(MAXFLOWSPD,FLOWSPEED) !check that there is some flow somewhere in this column of water
!       write(*,*) 'Flowspeed =',flowspeed,upstream,XPAR(IL(L)),YPAR(JL(L)),K,N
      IF(FLOWSPEED<1.0E-03)CYCLE !skip calcs in this layer if the flow is very slow (< 1 mm/s)
      IF(UPSTREAM==1)THEN !use the upstream flowspeed to assess power extraction
        UATVFACE= 0.25*(U(L,K)+U(LE,K)+U(LS,K)+U(LSE,K))  !u-velocity at south face (the v-face)
        VATUFACE= 0.25*(V(L,K)+V(LW,K)+V(LN,K)+V(LNW,K))  !v-velocity at west face  (the u-face)
        UATVFACEN=0.25*(U(L,K)+U(LE,K)+U(LN,K)+U(LNE,K))  !u-velocity at north face (the u-north-face)
        VATUFACEE=0.25*(V(L,K)+V(LE,K)+V(LN,K)+V(LNE,K))  !v-velocity at east face  (the v-east-face)
        FS_WF=U(L ,K);FS_EF=U(LE,K) !velocities on the west/east   faces (u velocities into the cell)
        FS_NF=V(LN,K);FS_SF=V(LS,K) !velocities on the north/south faces (v velocities into the cell)
        SPDN=SQRT(UATVFACEN*UATVFACEN+V(LN,K)*V(LN,K)) !speed at north face
        SPDS=SQRT(UATVFACE*UATVFACE+V(L,K)*V(L,K)) !speed at south face
        SPDE=SQRT(U(LE,K)*U(LE,K)+VATUFACEE*VATUFACEE) !speed at east face
        SPDW=SQRT(U(L,K)*U(L,K)+VATUFACE*VATUFACE) !speed at west face
        IF(FS_NF>-0.01)SPDN=0.0 !flow is OUT of north face
        IF(FS_SF<0.01)SPDS=0.0 !flow is OUT of south face
        IF(FS_WF<0.01)SPDW=0.0 !flow is OUT of west face
        IF(FS_EF>-0.01)SPDE=0.0 !flow is OUT of east face
        MAXSPD=MAX(SPDN,SPDS,SPDE,SPDW) !identify maximum speed
        VVELUP = MAX(SPDN,SPDS)
        IF(V(L,K)<0)VVELUP=-VVELUP
        UVELUP = MAX(SPDW,SPDE)
        IF(U(L,K)<0)UVELUP=-UVELUP
        IF(MAXSPD==SPDN)THEN !what face is it on?
          VELUP=SQRT((0.25*(U(LN,K)+U(LN+1,K)+U(LNC(LN),K)+U(LNC(LN)+1,K)))**2+V(LN,K)**2)
          FACE=1 !North
        ELSEIF(MAXSPD==SPDS)THEN
          VELUP=SQRT((0.25*(U(LS,K)+U(LS+1,K)+U(LSC(LS),K)+U(LSC(LS)+1,K)))**2+V(LS,K)**2)
          FACE=2 !South
        ELSEIF(MAXSPD==SPDE)THEN
          VELUP=SQRT(U(LE+1,K)**2+(0.25*(V(LE,K)+V(LE+1,K)+V(LN+1,K)+V(MIN(LC,LN+2),K)))**2)
          FACE=3 !East
        ELSE
          VELUP=SQRT(U(LW  ,K)**2+(0.25*(V(LW,K)+V(LW-1,K)+V(LN-1,K)+V(LN-2,K)))**2)
          FACE=4 !West
        ENDIF
      ELSEIF(UPSTREAM==0)THEN !use the local cell's flowspeed to assess power extraction
        UVELUP=0.5*(U(L,K)+U(LE,K))
        VVELUP=0.5*(V(L,K)+V(LN,K))
        VELUP=FLOWSPEED
      ENDIF
      IF(LAYFRACM(K)>0.0)THEN !MHK device exists in this layer 
        IF(VELUP<VMINCUT(M))THEN !no power generation
          THRSTCOEF=0.0 !no need for these calcs
        ELSEIF(VELUP<VMAXCUT(M))THEN !optimal power generation
          THRSTCOEF=CTMHK(M)
        ELSE !superoptimal flow speed limits power generation to VMAXCUT
          THRSTCOEF=CTMHK(M)
          VELUP=VMAXCUT(M)
        ENDIF
        FMHK=0.5*LAYFRACM(K)*THRSTCOEF*VELUP*HP(L)*DZC(K)*WIDTHMHK(M) !total "force" from MHK due to PMHK generation [m^3/s]
        PMHK(L,K)=FMHK*FLOWSPEED*FLOWSPEED !ThrustCoef*|u|u^2*area [m^5/s^3] (room for improvement as the area is ASSUMED square)
        AWEIGHTX1=DYU(L)*HU(L)/(DYU(L)*HU(L)+DYU(LE)*HU(LE)) !area-weighting for west face
        AWEIGHTX2=1.0-AWEIGHTX1 !area-weight for east face
        FXMHK(L ,K)=FXMHK(L ,K)+AWEIGHTX1*SUB(L )*FMHK*USPD/FLOWSPEED !SUB(L)*FMHK(L,K)*(Uvel/q) [m/s]
        FXMHK(LE,K)=FXMHK(LE,K)+AWEIGHTX2*SUB(LE)*FMHK*USPD/FLOWSPEED !distribute forces on each U-face of the cell
        AWEIGHTY1=DXV(L)*HV(L)/(DXV(L)*HV(L)+DXV(LN)*HV(LN)) !area-weight for south face
        AWEIGHTY2=1.0-AWEIGHTY1 !area-weight for north face
        FYMHK(L ,K)=FYMHK(L ,K)+AWEIGHTY1*SVB(L )*FMHK*VSPD/FLOWSPEED !x and y components of "forces" [m/s]
        FYMHK(LN,K)=FYMHK(LN,K)+AWEIGHTY2*SVB(LN)*FMHK*VSPD/FLOWSPEED !distribute forces on each V-face of the cell
        IF(BSC>0.0)THEN !if variable density, take it into account
          PMHK(L,K)=PMHK(L,K)*(B(L,K)+1.0)*1000.0
        ELSE
          PMHK(L,K)=PMHK(L,K)*1024. !density of seawater is ~1024kg/m^3
        ENDIF
        IF(DEBUG.AND.UPSTREAM/=0)WRITE(357,'(3(I2,1X),5(E12.5,1X))')IL(L), &
           JL(L),K,PMHK(L,K),VELUP, USPD,  VSPD,HP(L)
        IF(DEBUG.AND.     UPSTREAM/=0)WRITE(357,'(3(I2,1X),5(E12.5,1X))')IL(L), &
           JL(L),K,PMHK(L,K),VELUP,UVELUP,VVELUP,HP(L)
!         if(l==360)then
!            write(642,'(i1,6(1x,1pe11.4))'),k,uvelup,vvelup,velup,flowspeed,u(l,k),u(le,k)
!         endif
      ENDIF
      IF(LAYFRACS(K)>0.0)THEN !MHK support exists in this layer
        FSUP=0.5*LAYFRACM(K)*CDSUP(M)*FLOWSPEED*HP(L)*DZC(K)*WIDTHSUP(M) !calculate the force on a cell [m^3/s]
        PSUP(L,K)=FSUP*FLOWSPEED*FLOWSPEED !0.5*rho*C_d*|u|u^2 [m^5/s^3]
        AWEIGHTX1=DYU(L)*HU(L)/(DYU(L)*HU(L)+DYU(LE)*HU(LE)) !area-weighting for west face
        AWEIGHTX2=1.0-AWEIGHTX1 !area-weight for east face
        FXSUP(L ,K)=FXSUP(L ,K)+AWEIGHTX1*SUB(L )*FSUP*USPD/FLOWSPEED  ![m/s] 
        FXSUP(LE,K)=FXSUP(LE,K)+AWEIGHTX2*SUB(LE)*FSUP*USPD/FLOWSPEED  !distribute forces on both U-faces of the cell 
        AWEIGHTY1=DXV(L)*HV(L)/(DXV(L)*HV(L)+DXV(LN)*HV(LN)) !area-weight for south face
        AWEIGHTY2=1.0-AWEIGHTY1 !area-weight for north face
        FYSUP(L ,K)=FYSUP(L ,K)+AWEIGHTY1*SVB(L )*FSUP*VSPD/FLOWSPEED  ![m/s]
        FYSUP(LN,K)=FYSUP(LN,K)+AWEIGHTY2*SVB(LN)*FSUP*VSPD/FLOWSPEED  !distribute forces on both V-faces of the cell
        IF(BSC>0.0)THEN !if variable density, take it into account
          PSUP(L,K)=PSUP(L,K)*(B(L,K)+1.0)*1000.0
        ELSE
          PSUP(L,K)=PSUP(l,K)*1024.
       ENDIF
      ENDIF
    ENDDO
    IF(MAXFLOWSPD<1.0E-3)CYCLE !if there is no flow in this entire column of water, skip the calcs
    FXMHKE(L)=FXMHKE(L)+SUM(FXMHK(L,1:KC));FYMHKE(L)=FYMHKE(L)+SUM(FYMHK(L,1:KC)) !total MHK resistance in the column [m/s] used in FUHDXE & FVHDYE
    FXSUPE(L)=FXSUPE(L)+SUM(FXSUP(L,1:KC));FYSUPE(L)=FYSUPE(L)+SUM(FYSUP(L,1:KC)) !total support resistance in the column [m/s] used in FUHDXE & FVHDYE
    IF(FXMHKE(L)+FYMHKE(L)==0.0)CYCLE !At cold startup there is no power and this avoids NaNs
    DO K=1,KC
      IF(SUMLAYS==0.0)THEN
        FX(L,K)=FX(L,K)+PB_COEF*(FXMHKE(L)*(LAYFRACM(K)/SUMLAYM-NEGLAYFRACM(K)/&
SUMNEGLAYM))*SUM(U(L,1:KC)*DZC(1:KC)) !pull x-force out of MHK layer for internal mode (no support structure)
        FY(L,K)=FY(L,K)+PB_COEF*(FYMHKE(L)*(LAYFRACM(K)/SUMLAYM-NEGLAYFRACM(K)/&
SUMNEGLAYM))*SUM(V(L,1:KC)*DZC(1:KC)) !pull x-force out of MHK layer for internal mode (no support structure)
      ELSE
        FX(L,K)=FX(L,K)+(PB_COEF*FXMHKE(L)*(LAYFRACM(K)/SUMLAYM-NEGLAYFRACM(K)/&
SUMNEGLAYM)+FXSUPE(L)*(LAYFRACS(K)/SUMLAYS-NEGLAYFRACS(K)/SUMNEGLAYS))*SUM(U(L,1:KC)*DZC(1:KC)) !pull x-force out of MHK/support layer for internal mode
        FY(L,K)=FY(L,K)+(PB_COEF*FYMHKE(L)*(LAYFRACM(K)/SUMLAYM-NEGLAYFRACM(K)/&
SUMNEGLAYM)+FYSUPE(L)*(LAYFRACS(K)/SUMLAYS-NEGLAYFRACS(K)/SUMNEGLAYS))*SUM(V(L,1:KC)*DZC(1:KC)) !pull x-force out of MHK/support layer for internal mode
      ENDIF
    ENDDO
    FXMHKE(LE)=FXMHKE(LE)+SUM(FXMHK(LE,1:KC));FYMHKE(LN)=FYMHKE(LN)+SUM(FYMHK(LN,1:KC)) !Make sure to include adjacent cells as we only calculate on MHK cells
    FXSUPE(LE)=FXSUPE(LE)+SUM(FXSUP(LE,1:KC));FYSUPE(LN)=FYSUPE(LN)+SUM(FYSUP(LN,1:KC)) !Both x and y and north and east
    DO K=1,KC
      IF(SUMLAYS==0.0)THEN
        FX(LE,K)=FX(LE,K)+PB_COEF*(FXMHKE(LE)*(LAYFRACM(K)/SUMLAYM- & 
NEGLAYFRACM(K)/SUMNEGLAYM))*SUM(U(LE,1:KC)*DZC(1:KC)) !pull x-force out of MHK layer for internal mode
        FY(LN,K)=FY(LN,K)+PB_COEF*(FYMHKE(LN)*(LAYFRACM(K)/SUMLAYM- &
NEGLAYFRACM(K)/SUMNEGLAYM))*SUM(V(LN,1:KC)*DZC(1:KC)) !pull x-force out of MHK layer for internal mode
      ELSE
        FX(LE,K)=FX(LE,K)+(PB_COEF*FXMHKE(LE)*(LAYFRACM(K)/SUMLAYM-&
NEGLAYFRACM(K)/SUMNEGLAYM)+FXSUPE(LE)*(LAYFRACS(K)/SUMLAYS-&
NEGLAYFRACS(K)/SUMNEGLAYS))*SUM(U(LE,1:KC)*DZC(1:KC)) !pull x-force out of MHK/support layer for internal mode
        FY(LN,K)=FY(LN,K)+(PB_COEF*FYMHKE(LN)*(LAYFRACM(K)/SUMLAYM- &
NEGLAYFRACM(K)/SUMNEGLAYM)+FYSUPE(LN)*(LAYFRACS(K)/SUMLAYS- &
NEGLAYFRACS(K)/SUMNEGLAYS))*SUM(V(LN,1:KC)*DZC(1:KC)) !pull x-force out of MHK/support layer for internal mode
      ENDIF
    ENDDO
    EMHK(MHKCOUNT,L)=EMHK(MHKCOUNT,L)+DT*SUM(PMHK(L,1:KC))*2.7778E-10 !factor converts to MW-hr
    ESUP(MHKCOUNT,L)=ESUP(MHKCOUNT,L)+DT*SUM(PSUP(L,1:KC))*2.7778E-10 !factor converts to MW-hr
  ENDDO
  RETURN
END SUBROUTINE MHKPWRDIS



