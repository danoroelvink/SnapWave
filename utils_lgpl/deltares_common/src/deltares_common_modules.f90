module precision_basics
!----- LGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2022.                                
!                                                                               
!  This library is free software; you can redistribute it and/or                
!  modify it under the terms of the GNU Lesser General Public                   
!  License as published by the Free Software Foundation version 2.1.                 
!                                                                               
!  This library is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
!  Lesser General Public License for more details.                              
!                                                                               
!  You should have received a copy of the GNU Lesser General Public             
!  License along with this library; if not, see <http://www.gnu.org/licenses/>. 
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D" and "Deltares"    
!  are registered trademarks of Stichting Deltares, and remain the property of  
!  Stichting Deltares. All rights reserved.                                     
!                                                                               
!-------------------------------------------------------------------------------
!  $Id: precision_basics.F90 68181 2021-01-20 13:41:21Z leander $
!  $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/branches/research/Deltares/20210510_UNST-4987_4990_SO_3d_waves_morphology/src/utils_lgpl/deltares_common/packages/deltares_common/src/precision_basics.F90 $
!!--description-----------------------------------------------------------------
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
implicit none
!
! parameters, used in conversions: sp=single precision, hp=high (double) precision
!
integer, parameter :: sp=kind(1.0e00)
integer, parameter :: hp=kind(1.0d00)
!
! double precision integers:
!
integer, parameter :: long = SELECTED_INT_KIND(16)
!
! interfaces
!
  interface comparereal
     module procedure comparerealdouble
     module procedure comparerealsingle
  end interface
! 

contains

function comparerealdouble(val1, val2, eps)
!!--description-----------------------------------------------------------------
!
! Compares two double precision numbers
! Allow users to define the value of eps. If not, eps equals to the default machine eps
!
! Return value: -1 if val1 < val2
!                0 if val1 = val2
!               +1 if val1 > val2
!
!!--pseudo code and references--------------------------------------------------
!
! The functionality in this subroutine is copied from subroutine Ifdbl,
! written by Jaap Zeekant.
!
! eps must be machine precision dependent.
! eps may not be given by the user! See what happens when
! val1 = -666.0, val2 = -999.0, eps = 0.5
!
!!--declarations----------------------------------------------------------------
    implicit none
!
! Return value
!
integer :: comparerealdouble
!
! Global variables
!
real(hp), intent(in)           :: val1
real(hp), intent(in)           :: val2
real(hp), optional, intent(in) :: eps
!
! Local variables
!
real(hp) :: eps0
real(hp) :: value
!
!! executable statements -------------------------------------------------------
!
if (present(eps)) then
    eps0 = eps
else 
    eps0 = 2.0_hp * epsilon(val1)
endif
!
if (abs(val1)<1.0_hp .or. abs(val2)<1.0_hp) then
   value = val1 - val2
else
   value = val1/val2 - 1.0_hp
endif
!
if (abs(value)<eps0) then
   comparerealdouble = 0
elseif (val1<val2) then
   comparerealdouble = -1
else
   comparerealdouble = 1
endif
end function comparerealdouble



function comparerealsingle(val1, val2,eps)
!!--description-----------------------------------------------------------------
!
! REMARK: THE NAME OF THIS FUNCTION IS WRONG!
!         The name should be comparefp
!
! Compares two real numbers of type fp
! Allow users to define the value of eps. If not, eps equals to the default machine eps
!
! Return value: -1 if val1 < val2
!                0 if val1 = val2
!               +1 if val1 > val2
!
!!--pseudo code and references--------------------------------------------------
!
! The functionality in this subroutine is copied from subroutine Ifflt,
! written by Jaap Zeekant.
!
! eps must be machine precision dependent.
! eps may not be given by the user! See what happens when
! val1 = -666.0, val2 = -999.0, eps = 0.5
!
!!--declarations----------------------------------------------------------------
implicit none
!
! Return value
!
integer :: comparerealsingle
!
! Global variables
!
real(sp), intent(in)           :: val1
real(sp), intent(in)           :: val2
real(sp), optional, intent(in) :: eps
!
! Local variables
!
real(sp) :: eps0
real(sp) :: value
!
!! executable statements -------------------------------------------------------
!
!  
if (present(eps)) then
    eps0 = eps
else
    eps0 = 2.0_sp * epsilon(val1)
endif
!
if (abs(val1)<1.0_sp .or. abs(val2)<1.0_sp) then
   value = val1 - val2
else
   value = val1/val2 - 1.0_sp
endif
!
if (abs(value)<eps0) then
   comparerealsingle = 0
elseif (val1<val2) then
   comparerealsingle = -1
else
   comparerealsingle = 1
endif
end function comparerealsingle


end module precision_basics
   
   
   module precision
   !----- LGPL --------------------------------------------------------------------
   !
   !  Copyright (C)  Stichting Deltares, 2011-2022.
   !
   !  This library is free software; you can redistribute it and/or
   !  modify it under the terms of the GNU Lesser General Public
   !  License as published by the Free Software Foundation version 2.1.
   !
   !  This library is distributed in the hope that it will be useful,
   !  but WITHOUT ANY WARRANTY; without even the implied warranty of
   !  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   !  Lesser General Public License for more details.
   !
   !  You should have received a copy of the GNU Lesser General Public
   !  License along with this library; if not, see <http://www.gnu.org/licenses/>.
   !
   !  contact: delft3d.support@deltares.nl
   !  Stichting Deltares
   !  P.O. Box 177
   !  2600 MH Delft, The Netherlands
   !
   !  All indications and logos of, and references to, "Delft3D" and "Deltares"
   !  are registered trademarks of Stichting Deltares, and remain the property of
   !  Stichting Deltares. All rights reserved.
   !
   !-------------------------------------------------------------------------------
   !  $Id: precision.F90 68181 2021-01-20 13:41:21Z leander $
   !  $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/branches/research/Deltares/20210510_UNST-4987_4990_SO_3d_waves_morphology/src/utils_lgpl/deltares_common/packages/deltares_common/src/precision.F90 $
   !!--description-----------------------------------------------------------------
   !
   ! This module contains the parameters used to switch easily from
   ! single precision mode to double precision mode.
   !
   !
   ! See also precision.h file for C-code (DD)
   ! See also tri-dyn.igd file for connection with esm
   !
   ! sp: single precision
   ! hp: high (or double) precision
   ! fp: flexible precision, single or double
   !     fp is the used precision
   !
   ! SWITCHING FROM SINGLE PRECISION   FP   TO DOUBLE PRECISION:
   ! 1) File libsrc\flow_modsrc\precision.f90
   !    - Comment out the following line:
   !      INTEGER, PARAMETER :: FP=SP
   !    - Activate the following line:
   !      INTEGER, PARAMETER :: FP=HP
   ! 2) File include\flow\tri-dyn.igd
   !    - Comment out the following line:
   !      equivalence ( r(0),  rbuf(0))
   !    - Activate the following line:
   !      equivalence ( r(0),  dbuf(0))
   ! 3) File include\hydra\precision.h
   !    - Comment out the following line:
   !      #undef FLOW_DOUBLE_PRECISION
   !    - Activate the following line:
   !      #define FLOW_DOUBLE_PRECISION
   !
   ! SWITCHING FROM SINGLE PRECISION BODSED/DPS TO DOUBLE PRECISION:
   ! 1) File libsrc\flow_modsrc\precision.f90
   !    - Comment out the following line:
   !      integer, parameter :: prec=sp
   !    - Activate the following line:
   !      integer, parameter :: prec=hp
   ! 2) File include\flow\tri-dyn.igd
   !    - Comment out the following line:
   !      equivalence ( d(0),  rbuf(0))
   !    - Activate the following line:
   !      equivalence ( d(0),  dbuf(0))
   ! 3) File libsrc\flow_dd\hyexth\precision.h
   !    - Comment out the following line:
   !      #undef PREC_DOUBLE_PRECISION
   !    - Activate the following line:
   !      #define PREC_DOUBLE_PRECISION
   !
   !!--pseudo code and references--------------------------------------------------
   ! NONE
   !!--declarations----------------------------------------------------------------
   use precision_basics
   use iso_c_binding
   implicit none
   !
   ! fp is the generally used precision in Delft3D-FLOW
   !
   integer, parameter :: fp=hp
   !integer, parameter :: fp=sp
   !
   ! prec is used to switch bodsed/dps from sp to hp
   !
   integer, parameter :: prec=hp
   !integer, parameter :: prec=sp
   !
   ! old hp's in sobek that should stay sp
   !
   integer, parameter :: fhp=sp
   !
   ! length of integers which are esm/fsm pointers
   ! = 4 for 32bit versions
   ! = 8 for 64bit versions
   !
   integer, parameter :: pntrsize=c_size_t

   end module precision

   module mathconsts
   !----- LGPL --------------------------------------------------------------------
   !
   !  Copyright (C)  Stichting Deltares, 2011-2022.
   !
   !  This library is free software; you can redistribute it and/or
   !  modify it under the terms of the GNU Lesser General Public
   !  License as published by the Free Software Foundation version 2.1.
   !
   !  This library is distributed in the hope that it will be useful,
   !  but WITHOUT ANY WARRANTY; without even the implied warranty of
   !  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   !  Lesser General Public License for more details.
   !
   !  You should have received a copy of the GNU Lesser General Public
   !  License along with this library; if not, see <http://www.gnu.org/licenses/>.
   !
   !  contact: delft3d.support@deltares.nl
   !  Stichting Deltares
   !  P.O. Box 177
   !  2600 MH Delft, The Netherlands
   !
   !  All indications and logos of, and references to, "Delft3D" and "Deltares"
   !  are registered trademarks of Stichting Deltares, and remain the property of
   !  Stichting Deltares. All rights reserved.
   !
   !-------------------------------------------------------------------------------
   !  $Id: mathconsts.f90 68181 2021-01-20 13:41:21Z leander $
   !  $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/branches/research/Deltares/20210510_UNST-4987_4990_SO_3d_waves_morphology/src/utils_lgpl/deltares_common/packages/deltares_common/src/mathconsts.f90 $
   !!--description-----------------------------------------------------------------
   !
   ! This module defines some general mathematical constants like pi and
   ! conversion factors from degrees to radians and from (earth) days to
   ! seconds.
   !
   ! This module does NOT include physical constants like earth radius and
   ! gravity, or application dependent constants like missing values.
   !
   !!--pseudo code and references--------------------------------------------------
   ! NONE
   !!--declarations----------------------------------------------------------------
   use precision
   implicit none

   !
   ! single precision constants
   !
   real(sp), parameter :: ee_sp      = exp(1.0_sp)                      !< ee = 2.718281...
   real(sp), parameter :: pi_sp      = 4.0_sp*atan(1.0_sp)              !< pi = 3.141592...
   real(sp), parameter :: twopi_sp   = 8.0_sp*atan(1.0_sp)              !< 2pi
   real(sp), parameter :: sqrt2_sp   = sqrt(2.0_sp)                     !< sqrt(2)
   real(sp), parameter :: degrad_sp  = 4.0_sp*atan(1.0_sp)/180.0_sp     !< conversion factor from degrees to radians (pi/180)
   real(sp), parameter :: raddeg_sp  = 180.0_sp/(4.0_sp*atan(1.0_sp))   !< conversion factor from radians to degrees (180/pi)
   real(sp), parameter :: daysec_sp  = 24.0_sp*60.0_sp*60.0_sp          !< conversion factor from earth day to seconds
   real(sp), parameter :: yearsec_sp = 365.0_sp*24.0_sp*60.0_sp*60.0_sp !< conversion factor from earth year to seconds (non-leap)
   real(sp), parameter :: eps_sp     = epsilon(1.0_sp)                  !< epsilon for sp

   !
   ! flexible precision constants
   !
   real(fp), parameter :: ee         = exp(1.0_fp)                      !< ee = 2.718281.../2.71828182845904...
   real(fp), parameter :: pi         = 4.0_fp*atan(1.0_fp)              !< pi = 3.141592.../3.14159265358979...
   real(fp), parameter :: twopi      = 8.0_fp*atan(1.0_fp)              !< 2pi
   real(fp), parameter :: sqrt2      = sqrt(2.0_fp)                     !< sqrt(2)
   real(fp), parameter :: degrad     = 4.0_fp*atan(1.0_fp)/180.0_fp     !< conversion factor from degrees to radians (pi/180)
   real(fp), parameter :: raddeg     = 180.0_fp/(4.0_fp*atan(1.0_fp))   !< conversion factor from radians to degrees (180/pi)
   real(fp), parameter :: daysec     = 24.0_fp*60.0_fp*60.0_fp          !< conversion factor from earth day to seconds
   real(fp), parameter :: yearsec    = 365.0_fp*24.0_fp*60.0_fp*60.0_fp !< conversion factor from earth year to seconds (non-leap)
   real(fp), parameter :: eps_fp     = epsilon(1.0_fp)                  !< epsilon for fp

   !
   ! high precision constants
   !
   real(hp), parameter :: ee_hp        = exp(1.0_hp)                      !< ee = 2.71828182845904...
   real(hp), parameter :: pi_hp        = 4.0_hp*atan(1.0_hp)              !< pi = 3.14159265358979...
   real(hp), parameter :: twopi_hp     = 8.0_hp*atan(1.0_hp)              !< 2pi
   real(hp), parameter :: sqrt2_hp     = sqrt(2.0_hp)                     !< sqrt(2)
   real(hp), parameter :: degrad_hp    = 4.0_hp*atan(1.0_hp)/180.0_hp     !< conversion factor from degrees to radians (pi/180)
   real(hp), parameter :: raddeg_hp    = 180.0_hp/(4.0_hp*atan(1.0_hp))   !< conversion factor from radians to degrees (180/pi)
   real(hp), parameter :: daysec_hp    = 24.0_hp*60.0_hp*60.0_hp          !< conversion factor from earth day to seconds
   real(hp), parameter :: yearsec_hp   = 365.0_hp*24.0_hp*60.0_hp*60.0_hp !< conversion factor from earth year to seconds (non-leap)
   real(hp), parameter :: eps_hp       = epsilon(1.0_hp)                  !< epsilon for hp

   contains

   !> Obsolete initialization method.
   subroutine init_mathconsts()
   !
   end subroutine init_mathconsts

   end module mathconsts

   module physicalconsts
   !----- LGPL --------------------------------------------------------------------
   !
   !  Copyright (C)  Stichting Deltares, 2011-2022.
   !
   !  This library is free software; you can redistribute it and/or
   !  modify it under the terms of the GNU Lesser General Public
   !  License as published by the Free Software Foundation version 2.1.
   !
   !  This library is distributed in the hope that it will be useful,
   !  but WITHOUT ANY WARRANTY; without even the implied warranty of
   !  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   !  Lesser General Public License for more details.
   !
   !  You should have received a copy of the GNU Lesser General Public
   !  License along with this library; if not, see <http://www.gnu.org/licenses/>.
   !
   !  contact: delft3d.support@deltares.nl
   !  Stichting Deltares
   !  P.O. Box 177
   !  2600 MH Delft, The Netherlands
   !
   !  All indications and logos of, and references to, "Delft3D" and "Deltares"
   !  are registered trademarks of Stichting Deltares, and remain the property of
   !  Stichting Deltares. All rights reserved.
   !
   !-------------------------------------------------------------------------------
   !  $Id: physicalconsts.f90 68181 2021-01-20 13:41:21Z leander $
   !  $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/branches/research/Deltares/20210510_UNST-4987_4990_SO_3d_waves_morphology/src/utils_lgpl/deltares_common/packages/deltares_common/src/physicalconsts.f90 $
   !!--description-----------------------------------------------------------------
   !
   ! This module defines some general physical constants
   !
   !!--pseudo code and references--------------------------------------------------
   ! NONE
   !!--declarations----------------------------------------------------------------
   use precision
   implicit none
   private
   !
   ! high precision constants
   !

   real(kind=hp), parameter, public :: earth_radius = 6378137_hp        !< earth radius (m)
   real(kind=hp), parameter, public :: dtol_pole    = 0.0001_hp         !< pole tolerance in degrees
   real(kind=hp), parameter, public :: CtoKelvin    = 273.15_hp         !< conversion offset between Celsius and Kelvin
   real(kind=hp), parameter, public :: stf          = 5.6705085e-8_hp   !< Stefan's constant =5.6705085e-8 [W/m^2/K^4]
   !! (see 19308-part-iv-physical-processes.pdf from ECMWF;
   !!  it differs slightly from the value after the redefinition of SI in 2019)

   end module physicalconsts