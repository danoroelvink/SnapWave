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
module m_ec_triangle           ! original name : m_triangle
   use precision, only : hp
   implicit none
   real(kind=hp)   , allocatable :: XCENT(:), YCENT(:)
   integer, allocatable          :: INDX(:,:)
   integer, allocatable          :: EDGEINDX(:,:)
   integer, allocatable          :: TRIEDGE(:,:)
   integer                       :: NUMTRI
   integer                       :: NUMTRIINPOLYGON
   integer                       :: NUMEDGE
   integer, parameter            :: ITYPE = 2 ! 1 = ORIGINAL FORTRAN ROUTINE, 2 = NEW C ROUTINE

   integer                       :: jagetwf = 0    ! if 1, also assemble weightfactors and indices in:
   integer, allocatable          :: indxx(:,:)     ! to be dimensioned by yourselves 3,*
   real(kind=hp)   , allocatable :: wfxx (:,:)

   type t_nodi
      integer                    :: NUMTRIS       ! total number of TRIANGLES ATtached to this node
      integer, allocatable       :: TRINRS(:)     ! numbers of ATTACHED TRIANGLES
   end type t_nodi

   type (t_nodi), dimension(:), allocatable :: NODE

   integer, dimension(:,:), allocatable :: LNtri  ! triangles connected to edges, dim(2,numedges)

   integer                              :: IDENT   ! identifier
   integer, dimension(:),   allocatable :: imask   ! mask array for triangles

end module m_ec_triangle