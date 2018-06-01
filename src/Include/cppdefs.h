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
!
!      ________________________________________________________
!     |                                                        |
!     |   All the test cases                                   |
!     |________________________________________________________|

#      undef KeyTestOnlyPoisson /* 3D Poisson equation*/
#      undef KeyTaylorVortex    /* NS: 3D Taylor vortex problem*/
#      undef KeyStaticCylinder  /* NS: Static cylinder problem*/
#      undef KeyStaticChannel   /* Channel flow with periodic BC*/
#      undef KeyStandingWave    /* FS: standing wave problem*/
#      undef KeyEstuaryGironde  /* Application: Estuary Gironde*/

!      ________________________________________________________
!     |                                                        |
!     |   Compiler                                             |
!     |________________________________________________________|

#      undef KeyXLF             /* XLF compiler */
#      undef KeyG77             /* GNU compiler */
!      ________________________________________________________
!     |                                                        |
!     |   Parallelization option                               |
!     |________________________________________________________|

#      undef KeyParallel          /* Parallelization */
#      undef KeyDisplayParallel   /* Display parallel commends */
#      undef KeyCUDA              /* Parallelization by CUDA */
!      ________________________________________________________
!     |                                                        |
!     |   Mesh and Domain                                      |
!     |________________________________________________________|

#      undef KeyMeshType1       /* Mesh type 1: square divided in two triangles */
#      undef KeyMeshType2       /* Mesh type 2: square divided in four triangles */
#      undef KeyMeshType3       /* Mesh type 3: BlueKenue unstructured grid */
#      undef KeyMesh1Example3   /* General Mesh designed for the Example 3: NN=16 */
#      undef KeyMesh2Example3   /* General Mesh designed for the Example 3: NN=32 */
#      undef KeyMesh3Example3   /* General Mesh designed for the Example 3: NN=64 */
#      undef KeyMesh4Example3   /* General Mesh designed for the Example 3: NN=128*/
#      undef KeyDomainGeneral   /* Mesh in a general domain */
#      undef KeyMeshFreeSurface /* Mesh used for the free-surface example */

!      ________________________________________________________
!     |                                                        |
!     |   Free-surface                                         |
!     |________________________________________________________|

#      undef KeyFixedFreeSurface  /* Rigid surface */
#      undef KeyFroudeEta         /* Froude number used (otherwise gravity) */
#      undef KeyMaximumEta        /* Set a maximum value of eta */
#      undef KeyFluxLimiterWL     /* FluxLimiter for the Water-Level */
#      undef KeyStanWaveExample2D  /* Standing wave Example 2D */
#      undef KeyStanWaveExample3D  /* Standing wave Example 3D */
!      ________________________________________________________
!     |                                                        |
!     |   Advection-Diffusion Scheme                           |
!     |________________________________________________________|

#      undef KeyAdvImplicit       /* Solving advection Implicit */
#      undef KeyAdvSemiImplicit   /* Solving advection Semi-implicit */
#      undef KeyAdvExplicit       /* Solving advection Explicit */
#      undef KeyAdvFullExplicit   /* Solving advection Full Explicit */

#      undef KeyDiffExplicit      /* Solving diffusion explicit */
#      undef KeyDiffImplicit      /* Solving diffusion implicit */
#      undef KeyDiffFullExplicit  /* Solving advection full explicit */

!      ----------------------------
!      Diffusion coefficients
#      undef KeyDiffusiveRe       /* Using only Reynold's number*/
#      undef KeyDiffusiveLES      /* Using turbulence LES*/
#      undef KeyDiffusiveDNS      /* Using turbulence DNS*/

!      ----------------------------
!      Turbulence LES
#      undef  KeyLESWallSigma    /* LES: wall distance = abs(sigma) */
#      undef  KeyLESWallGeneral  /* LES: wall distance general calculation */
!      .......................
#      undef KeyLESOptEqual      /* LES: same scalar for all coeff. */
#      undef KeyLESOptDifferent  /* LES: different scalar for each coeff. */
!      .......................
#      undef KeyLESMason         /* LES: Mason Damping function */
#      undef KeyLESVan           /* LES: Van Driest Damping function */

!      ________________________________________________________
!     |                                                        |
!     |   Boundary conditions (Jan 2018)                       |
!     |________________________________________________________|

!      ----------------------------
!      Vertical: Bottom
#      undef KeyBCbottomSlip
#      undef KeyBCbottomNoSlip
#      undef KeyBCbottomExact_1st
#      undef KeyBCbottomExact_2nd
!      ----------------------------
!      Vertical: Top
#      undef KeyBCtopSlip
#      undef KeyBCtopNoSlip
#      undef KeyBCtopExact_1st
#      undef KeyBCtopExact_2nd
!      ----------------------------
!      Horizontal: Wall
#      undef KeyBCwallXSlip
#      undef KeyBCwallYSlip
#      undef KeyBCwallNoSlip
#      undef KeyBCwallExact
#      undef KeyBCwallImpermeability
!      ----------------------------
!      Horizontal: Inflow & Outflow
#      undef KeyBCinflowExact
#      undef KeyBCoutflowNeumann
!      ----------------------------
!      Periodic
#      undef KeyBCperiodicX          /* wall.up <=> wall.down */
#      undef KeyBCperiodicY          /* inflow  <=> outflow   */
#      undef KeyBCperiodicZ          /* top     <=> bottom    */
!      ----------------------------
!      Inflow
#      undef KeyInflowConstant       /* Inflow Constant profile */
#      undef KeyInflowPoiseuille     /* Inflow Poiseuille profile */
!      ________________________________________________________
!     |                                                        |
!     |   Decoupling Pressure-Velocity Methods                 |
!     |________________________________________________________|

!      ----------------------------
!      Dynamic pressure calculation
#      undef DynamicPressure
!      ----------------------------
!      Pressure decoupling method
#      undef Key_PredictorCorrector  /* Predictor-Corrector with div. */
#      undef Key_ProjectionFirst     /* Projection First order */
#      undef Key_ProjectionSecond    /* Projection Second order */
!      ----------------------------
!      Pressure boundary condition
#      undef Key_DirichletBCp        /* Dirichlet BC for the pressure */
#      undef Key_NeumannBCp          /* Neumann   BC for the pressure */
#      undef Key_MixBCp              /* Mix: Neumann + Dirichlet top BC for p */
#      undef Key_PeriodicBCp         /* Periodic BC for p (XY direction) */
#      undef KeyUseNeumannBCvetex    /* Apply Neumann BC at vertex points */

!      ________________________________________________________
!     |                                                        |
!     |   Linear system solvers                                |
!     |________________________________________________________|

#      undef KeyIterationWithVertex /* Vertex update at the pressure iterations*/
#      undef KeySOR             /* Succesive-Over-Relaxation */
#      undef KeySiOR            /* Simultaneous-Over-Relaxation */
#      undef KeyJSOR            /* Partitioning Jacobi SOR */
#      undef KeyPSOR            /* Partition Succesive-Over-Relaxation */
#      undef KeyMultiSOR        /* Multicoloring Succesive-Over-Relaxation */
#      undef KeyGMRES           /* GMRES solver */
#      undef KeyPseudoTime      /* Pseudo time solver */
#         undef KeyAuxJacobi       /* Auxiliar functions*/
#         undef KeyAuxJSOR
#         undef KeyAuxMCSOR
#         undef KeyAuxPSOR
#      undef KeySchwarzSOR     /* Schwarz solver */

#     undef KeySORcuda           /* CUDA SOR solver */
#     undef KeyPDSORcuda         /* CUDA PD-SOR solver */
#     undef KeyJacobi_cuda       /* CUDA_Jacobi solver */
#     undef KeyJSOR_cuda         /* CUDA-JSOR solver */
#     undef KeyMSOR_cuda         /* CUDA-Multicolor solver */

!      ________________________________________________________
!     |                                                        |
!     |   Computing options and parameters                     |
!     |________________________________________________________|

#      undef KeyMask            /* Mask velocity on a predifined area of the mesh */
#      undef KeyFluxLimiter     /* Using the flux limiter in the edge approx */
#      undef KeyRigid           /* No free surface calculation */
#      undef KeyOrderingIndex   /* Re-order cell-center index using the regions*/
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
!     |   Display options                                      |
!     |________________________________________________________|

#      undef KeyDbg             /* Print information during execution */
#      undef KeyDisplay         /* Print lines-information during execution */
#      undef KeySaveVaryEps     /* Save results varying different tolerance=eps */
#      undef KeySaveMatrixA     /* Save the matrix A to know how sparse it is*/
#      undef KeySave1DReference /* Save 1D results for free-surface */
#      undef KeySave1DRef3Points/* Save 1D results for FS 3 points */
!      ________________________________________________________
!     |                                                        |
!     |   Save options                                         |
!     |________________________________________________________|

#      undef KeySave1DReference     /* Save one-dimensional data */
#      undef KeySaveStatistics      /* Save turbulence statistics */
#      undef KeyReStartStatistics   /* Restart saving turbulence statistics */

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                           PROBLEM CHOICE                            !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |                  Choose a compiler                     |
!     |________________________________________________________|

#     define KeyG77

!      ________________________________________________________
!     |                                                        |
!     |   ***********    Choose a test Case   **************   |
!     |________________________________________________________|

!#     define KeyEstuaryGironde
#      define KeyStaticChannel    /* Serial,MPI,CUDA*/
!#     define KeyStaticCylinder   /* Serial,MPI,CUDA*/
!#     define KeyStandingWave     /* Serial,MPI,CUDA*/
!#     define KeyTaylorVortex     /* Serial,MPI,CUDA*/
!#     define KeyTestOnlyPoisson  /* Serial,MPI,CUDA*/

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                 General Options for all problems                    !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |                     Run in parallel                    |
!     |________________________________________________________|

#      define KeyParallel
!#     define KeyCUDA
!      -------------
!#     define KeyOrderingIndex
!      ________________________________________________________
!     |                                                        |
!     |                Linear Solver (pressure)                |
!     |________________________________________________________|

#      define KeyIterationWithVertex
!      -------------
#      if defined(KeyParallel)
#         define KeyMultiSOR
!#        define KeyJSOR
!#        define KeyPSOR
!#        define KeySchwarzSOR
!#        define KeyJacobi
!      -------------
#      elif defined(KeyCUDA)
#         define KeyMSOR_cuda
!#        define KeyJacobi_cuda
!#        define KeyJSOR_cuda
!#        define KeyPDMSOR_cuda
!      -------------
#      else
!#        define KeyMultiSOR
#         define KeySOR
!#        define KeySiOR
!#        define KeyGMRES
#      endif
!      -------------
!#     define KeyPseudoTime
!#         define KeyAuxMCSOR
!#         define KeyAuxPSOR
!#         define KeyAuxJSOR
!#         define KeyAuxJacobi
!      ________________________________________________________
!     |                                                        |
!     |                  Type of interpolation                 |
!     |________________________________________________________|

#     define KeyInterpoNew
!#    define KeyInterpoNewPressureLSM
!     ---------------------
!#    define KeyInterpoDist
!#    define KeyInterpoLSM
!#    define KeyInterpoMix
!#    define KeyUseInterGhost /*Use ghost points during interpolation*/

!      ________________________________________________________
!     |                                                        |
!     |           Other options (saving & displaying)          |
!     |________________________________________________________|

!     --------------------------------------------------
!     Save extra information
#            define KeySave1DReference
!#            define KeySaveVaryEps
!#            define KeySaveMatrixA
!     --------------------------------------------------
!     Domain (only use to save data):
#            define KeyDomainGeneral
!#            define KeyMeshType2
!#            define KeyMesh3Example3
!#            define KeyMeshFreeSurface
!     --------------------------------------------------
!     Print information during execution
!#            define KeyDbg
!#            define KeyDisplay
!#            define KeyDisplayParallel

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                   Options for each problem case.                    !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |                    ESTUARY GIRONDE                     |
!     |________________________________________________________|

#     if defined(KeyEstuaryGironde)
!        ------------------------------------------------------
!        [1] Free-surface
#        define KeyFixedFreeSurface
!#       define KeyFroudeEta
!        ----------------------------
!        Free-surface: Advection
!#       define KeyFreeFullExplicit
#        define KeyFreeExplicit
!#       define KeyFreeImplicit
!            .......................
#            ifdef KeyFreeImplicit
#               define KeyFS_PoissonJSOR
!#              define KeyFS_PoissonJacobi
#            endif
!            .......................
!#       define KeyFluxLimiterWL
!        ------------------------------------------------------
!        [2] Velocity: Advection
!#       define KeyFluxLimiter
!#       define KeyAdvFullExplicit
!#       define KeyAdvExplicit
#        define KeyAdvSemiImplicit
!        ----------------------------
!        Velocity: Diffusion coeff. & turbulence
#        define KeyDiffusion
#        define KeyDiffusiveRe
!        ----------------------------
!        BC Velocity
!#       define KeyBCbottomNoSlip
#        define KeyBCtopSlip
#        define KeyBCwallNoSlip
#        define KeyBCinflowExact
#        define KeyBCoutflowNeumann
!            .......................
#            ifdef KeyBCinflowExact
#               define KeyInflowConstant
!#              define KeyInflowPoiseuille
#            endif
!            .......................
!        ------------------------------------------------------
!        [3] Pressure: dynamic pressure calculation
#        define DynamicPressure
!        ----------------------------
!        BC Pressure
#        define Key_MixBCp
!        ----------------------------
!        Decoupling: pressure & vel
#        define Key_ProjectionFirst
!#       define Key_ProjectionSecond
!#       define Key_PredictorCorrector
#     endif

!      ________________________________________________________
!     |                                                        |
!     |                    STATIC CHANNEL                      |
!     |________________________________________________________|

#     if defined(KeyStaticChannel)
!        ------------------------------------------------------
!        [1] Free-surface
#        define KeyFixedFreeSurface
!#       define KeyFroudeEta
!        ------------------------------------------------------
!        [2] Velocity: Advection
#        define KeyAdvSemiImplicit
!        ----------------------------
!        Velocity: Diffusion coefficient & turbulence
#        define KeyDiffusion
#        define KeyDiffusiveRe
!        ----------------------------
!        BC Velocity
#        define KeyBCbottomNoSlip
#        define KeyBCtopSlip /* KeyBCtopSlip or KeyBCtopNoSlip */
#        define KeyBCperiodicX /* wall.up <=> wall.down */
#        define KeyBCperiodicY /* inflow  <=> outflow   */
!            .......................
#            ifdef KeyBCinflowExact
#               define KeyInflowConstant
!#              define KeyInflowPoiseuille
#            endif
!            .......................
!        ------------------------------------------------------
!        [3] Pressure: dynamic pressure calculation
#        define DynamicPressure
!        ----------------------------
!        Decoupling: pressure & vel
#        define Key_ProjectionFirst
!        ----------------------------
!        BC Pressure
!#       define Key_NeumannBCp
#        define Key_PeriodicBCp
!        ----------------------------
!        Body force
#        define KeyBodyForce
!        ------------------------------------------------------
!        [4] Other: saving
#        define KeySaveStatistics
!#       define KeyReStartStatistics
#     endif

!      ________________________________________________________
!     |                                                        |
!     |                    STATIC CYLINDER                     |
!     |________________________________________________________|

#     if defined(KeyStaticCylinder)
!        ------------------------------------------------------
!        [1] Free-surface
#        define KeyFixedFreeSurface
!#       define KeyMaximumEta
#        define KeyFroudeEta
!        ----------------------------
!        Free-surface: Advection
!#       define KeyFreeFullExplicit
#        define KeyFreeExplicit
!#       define KeyFreeImplicit
!            .......................
#            ifdef KeyFreeImplicit
#               define KeyFS_PoissonJSOR
!#              define KeyFS_PoissonJacobi
#            endif
!            .......................
!#       define KeyFluxLimiterWL
!        ------------------------------------------------------
!        [2] Velocity: Advection
!#       define KeyFluxLimiter
!#       define KeyAdvFullExplicit
!#       define KeyAdvExplicit
#        define KeyAdvSemiImplicit
!        ----------------------------
!        Velocity: Diffusion coefficient & turbulence
#        define KeyDiffusion
#        define KeyDiffusiveRe
!#       define KeyDiffusiveLES
!            .......................
#            ifdef KeyDiffusiveLES
#               define  KeyLESWallSigma
!#              define  KeyLESWallGeneral
!               ------
#               define KeyLESOptEqual
!#              define KeyLESOptDifferent
!               ------
#               define KeyLESMason
!#              define KeyLESVan
#            endif
!            .......................
!        ----------------------------
!        Velocity: BC
!#       define KeyBCbottomSlip
#        define KeyBCbottomNoSlip
#        define KeyBCtopSlip
#        define KeyBCwallXSlip
#        define KeyBCoutflowNeumann
#        define KeyBCinflowExact
!            .......................
#            ifdef KeyBCinflowExact
#               define KeyInflowConstant
!#              define KeyInflowPoiseuille
#            endif
!            .......................
!        ------------------------------------------------------
!        [3] Pressure: dynamic pressure calculation
#        define DynamicPressure
!        ----------------------------
!        Pressure: BC
#        ifdef KeyFixedFreeSurface
#           define Key_NeumannBCp
#        else
#           define Key_MixBCp
#        endif
!        ----------------------------
!        Decoupling: pressure & vel
!#       define Key_PredictorCorrector /* Not working so far*/
#        define Key_ProjectionFirst
!#       define Key_ProjectionSecond
!        ------------------------------------------------------
!        [4] Other: saving
#        define KeySave1DRef3Points /* 3 reference FS points */
#     endif

!      ________________________________________________________
!     |                                                        |
!     |                STANDING WAVE TEST CASE                 |
!     |________________________________________________________|

#     if defined(KeyStandingWave)
#        define KeyStanWaveExample2D
!#       define KeyStanWaveExample3D
!        ------------------------------------------------------
!        [1] Free-surface
!#       define KeyFixedFreeSurface
!#       define KeyMaximumEta
!        ----------------------------
!        Free-surface: Advection
!#       define KeyFreeFullExplicit
#        define KeyFreeExplicit
!#       define KeyFreeImplicit
!            .......................
#            ifdef KeyFreeImplicit
#               define KeyFS_PoissonJSOR
!#              define KeyFS_PoissonJacobi
#            endif
!            .......................
!#       define KeyFluxLimiterWL
!        ------------------------------------------------------
!        [2] Velocity: Advection
!#       define KeyFluxLimiter
!#       define KeyAdvFullExplicit
!#       define KeyAdvExplicit
#        define KeyAdvSemiImplicit
!        ----------------------------
!        Velocity: Diffusion coeff. & turbulence
!#       define KeyDiffusion
!#       define KeyDiffusiveRe
!        ----------------------------
!        BC Velocity
#        define KeyBCbottomSlip
#        define KeyBCtopSlip
#        define KeyBCwallImpermeability
!        ------------------------------------------------------
!        [3] Pressure: dynamic pressure calculation
#        define DynamicPressure
!        ----------------------------
!        BC Pressure
#        define Key_MixBCp
!        ----------------------------
!        Decoupling: pressure & vel
!#       define Key_PredictorCorrector
#        define Key_ProjectionFirst
!#       define Key_ProjectionSecond
#     endif

!      ________________________________________________________
!     |                                                        |
!     |                     TAYLOR VORTEX                      |
!     |________________________________________________________|

#     if defined(KeyTaylorVortex)
!        ------------------------------------------------------
!        [1] Free-surface
#        define KeyFixedFreeSurface
!        ------------------------------------------------------
!        [2] Velocity: Advection
!#       define KeyFluxLimiter
!#       define KeyAdvFullExplicit
#        define KeyAdvSemiImplicit
!#       define KeyAdvExplicit  /*Bad approximations*/
!#       define KeyAdvImplicit  /*Not ready*/
!        ----------------------------
!        Velocity: Diffusion coefficient
#        define KeyDiffusion
#        define KeyDiffusiveRe
!        ----------------------------
!        BC Velocity
#        define KeyBCbottomExact_2nd
#        define KeyBCtopExact_2nd
#        define KeyBCwallExact
!        ------------------------------------------------------
!        [3] Pressure: dynamic pressure calculation
#        define DynamicPressure
!        ----------------------------
!        BC Pressure
#        define Key_DirichletBCp
!#       define Key_NeumannBCp
!#       define Key_MixBCp
!        ----------------------------
!        Decoupling: pressure & vel
#        define Key_PredictorCorrector
!#       define Key_ProjectionFirst
!#       define Key_ProjectionSecond
#     endif

!      ________________________________________________________
!     |                                                        |
!     |                   TEST ONLY POISSON                    |
!     |________________________________________________________|

#     if defined(KeyTestOnlyPoisson)
!        ------------------------------------------------------
!        [1] Free-surface
#        define KeyFixedFreeSurface
!        ------------------------------------------------------
!        [3] Pressure: dynamic pressure calculation
#        define DynamicPressure
!        ----------------------------
!        BC Pressure
#        define Key_DirichletBCp
!#       define Key_NeumannBCp
!#       define Key_MixBCp
!#       define Key_PeriodicBCp /* Not tested yet */
!#       define KeyUseNeumannBCvertex
#     endif

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                                  END                                !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

