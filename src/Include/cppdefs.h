!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                Choose a pre-defined model application               !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!---------------------------------------------------------------------!
!                                                                     !
!    This part of the program is called using the command:            !
!                    #include "cppdefs.h"                             !
!                                                                     !
!    Choose the C-preprocessing options by using the command          !
!                    #define   to activate option or                  !
!                    #undef    to deactivate option.                  !
!                                                                     !
!---------------------------------------------------------------------!
!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!      Detailed description of all available CPP options.             !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
! Note: Many definition are out of date, see each key for details.    !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!
!      ________________________________________________________
!     |                                                        |
!     |   All the test cases                                   |
!     |________________________________________________________|

#      undef KeyTESTChannel         /* Application to clapage test-case*/
!      ________________________________________________________
!     |                                                        |
!     |   Parallelization option                               |
!     |________________________________________________________|

#      undef KeyParallel          /* Parallelization */
#      undef KeyDisplayParallel   /* Display parallel commends */
!      ________________________________________________________
!     |                                                        |
!     |   Interpolation 2D                                     |
!     |________________________________________________________|

#      undef KeyInterpoNew      /* Vertex interpolation: distance weighting, new formulation Jul2014*/
#      undef KeyInterpoArea     /* Vertex interpolation: area weighting*/
#      undef KeyInterpoDist     /* Vertex interpolation: distance C-V weighting*/
#      undef KeyInterpoDist2    /* Vertex interpolation: distance C-V weighting more points*/
#      undef KeyInterpoLSM      /* Vertex interpolation: using the LSM*/
#      undef KeyInterpoMix      /* Vertex interpolation: using a mix of LSM & dist*/
#      undef KeyUseInterGhost   /* Using ghost points during vertex interpolation*/
!      ________________________________________________________
!     |                                                        |
!     |   Matrix for Pressure Poisson Equation                 |
!     |________________________________________________________|
#      undef KeyPPECenter
#      undef KeyPPERHS
!      ________________________________________________________
!     |                                                        |
!     |   Display options                                      |
!     |________________________________________________________|

#      undef KeyDbg             /* Print information during execution */
#      undef KeyDbgPBC          /* Print information during execution */
#      undef KeyDisplay         /* Print lines-information during execution */
!      ________________________________________________________
!     |                                                        |
!     |  Turbulence models                                     |
!     |________________________________________________________|
#      undef KeyLES   /* large eddy simulation */
#      undef KeyLESVan /* Van Driest near wall damping function */
!      ________________________________________________________
!     |                                                        |
!     |  Boundary conditions                                   |
!     |________________________________________________________|
#      undef KeyTESTpBC /* periodic bc support*/


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                           PROBLEM CHOICE                            !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |                  Choose a test case                    |
!     |________________________________________________________|
#            define KeyTESTChannel
!     --------------------------------------------------
!     Run in parallel:
#            define KeyParallel
!     --------------------------------------------------
!     Type of interpolation:
#             define KeyInterpoNew
!     --------------------------------------------------
!     Type of Momentum Discretization:
!#             define KeyImplicit
!#             define KeyAdvUpwind
!#             define KeyOldPressure
!#            define  KeyFullStress
#             define KeyKimChoi
!#             define KeyMahesh
!     --------------------------------------------------
!     Type of PPE Matrix:
#             define KeyPPECenter
#             define KeyPPERHS
!     --------------------------------------------------
!     Print information during execution
!#            define KeyDbg
!#            define KeyDbgPBC
!#            define KeyDisplay
!#            define KeyDisplayParallel
!     --------------------------------------------------
!     Turbulence Models:
!#            define KeyLES
!#            define KeyLESVan
!     --------------------------------------------------
!     Boundary condition
#             define KeyTESTpBC
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                                  END                                !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

