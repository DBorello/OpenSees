# Microsoft Developer Studio Project File - Name="material" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=material - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "material.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "material.mak" CFG="material - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "material - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "material - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "material - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib"
# PROP Intermediate_Dir "..\..\obj\material"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "c:\Program Files\tcl\include" /I "..\..\..\src\renderer" /I "..\..\..\src\material\nD\soil" /I "..\..\..\src\material\nD\template3dep" /I "..\..\..\src\recorder\response" /I "..\..\..\src\material\backbone" /I "..\..\..\src\material\state" /I "..\..\..\src\material\state\strength" /I "..\..\..\src\material\state\deformation" /I "..\..\..\src\material\state\stiffness" /I "..\..\..\src\material\section\repres\section" /I "..\..\..\src\material\section\repres\cell" /I "..\..\..\src\material\section\repres\patch" /I "..\..\..\src\material\section\repres\reinfBar" /I "..\..\..\src\material\section\repres\reinfLayer" /I "..\..\..\src\material\section\fiber" /I "..\..\..\src\element\nonlinearBeamColumn\fiber" /I "..\..\..\src\element\nonlinearBeamColumn\matrixutil" /I "..\..\..\src\material\section" /I "..\..\..\src\handler" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\material\nD" /I "..\..\..\src\element" /I "..\..\..\src\actor\channel" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\matrix" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\tagged" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\domain\component" /I "..\..\..\src\material" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\material\nd\template3dep" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "material - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib"
# PROP Intermediate_Dir "..\..\obj\material"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "c:\Program Files\tcl\include" /I "..\..\..\src\renderer" /I "..\..\..\src\material\nD\soil" /I "..\..\..\src\material\nD\template3dep" /I "..\..\..\src\recorder\response" /I "..\..\..\src\material\backbone" /I "..\..\..\src\material\state" /I "..\..\..\src\material\state\strength" /I "..\..\..\src\material\state\deformation" /I "..\..\..\src\material\state\stiffness" /I "..\..\..\src\material\section\repres\section" /I "..\..\..\src\material\section\repres\cell" /I "..\..\..\src\material\section\repres\patch" /I "..\..\..\src\material\section\repres\reinfBar" /I "..\..\..\src\material\section\repres\reinfLayer" /I "..\..\..\src\material\section\fiber" /I "..\..\..\src\element\nonlinearBeamColumn\fiber" /I "..\..\..\src\element\nonlinearBeamColumn\matrixutil" /I "..\..\..\src\material\section" /I "..\..\..\src\handler" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\material\nD" /I "..\..\..\src\element" /I "..\..\..\src\actor\channel" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\matrix" /I "..\..\..\src" /I "..\..\..\src\actor\actor" /I "..\..\..\src\tagged" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\domain\component" /I "..\..\..\src\material" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\material\nd\template3dep" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /D "MHS_INSURE" /FR /FD /GZ /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "material - Win32 Release"
# Name "material - Win32 Debug"
# Begin Group "uniaxial"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\Concrete01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\Concrete01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ElasticMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ElasticMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ElasticPPMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ElasticPPMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\EPPGapMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\EPPGapMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\HardeningMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\HardeningMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\HystereticMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\HystereticMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ParallelMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ParallelMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\PathIndependentMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\PathIndependentMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\SeriesMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\SeriesMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\Steel01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\Steel01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\TclModelBuilderUniaxialMaterialCommand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\UniaxialMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\UniaxialMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ViscousMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\uniaxial\ViscousMaterial.h
# End Source File
# End Group
# Begin Group "nD"

# PROP Default_Filter ""
# Begin Group "elasticIsotropic"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropic3D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropic3D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicPlaneStrain2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicPlaneStrain2D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicPlaneStress2D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\ElasticIsotropicPlaneStress2D.h
# End Source File
# End Group
# Begin Group "j2Plasticity"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2AxiSymm.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2AxiSymm.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2PlaneStrain.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2PlaneStrain.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2PlaneStress.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2PlaneStress.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2Plasticity.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2Plasticity.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2ThreeDimensional.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\J2ThreeDimensional.h
# End Source File
# End Group
# Begin Group "interactionModels"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\BidirectionalMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\BidirectionalMaterial.h
# End Source File
# End Group
# Begin Group "soilModels"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\FluidSolidPorousMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\FluidSolidPorousMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\MultiYieldSurface.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\MultiYieldSurface.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\PressureDependMultiYield.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\PressureDependMultiYield.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\PressureIndependMultiYield.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\PressureIndependMultiYield.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\T2Vector.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\T2Vector.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\soil\TclUpdateMaterialStageCommand.cpp
# End Source File
# End Group
# Begin Group "template3dep"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\DP_PS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\DP_PS.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\DP_YS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\DP_YS.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\DP_YS01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\DP_YS01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_LEeq.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_LEeq.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_LEij.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_LEij.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_NLEeq.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_NLEeq.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_NLEij.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_NLEij.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_NLEijMD.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_NLEijMD.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_S.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_S.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_T.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EL_T.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EPState.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\EPState.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\MatPoint3D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\MatPoint3D.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\MD_PS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\MD_PS.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\MD_PS01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\MD_PS01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\MD_YS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\MD_YS.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\Template3Dep.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\Template3Dep.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\VM_PS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\VM_PS.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\VM_YS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\VM_YS.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\YS.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\Template3Dep\YS.h
# End Source File
# End Group
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\NDMaterial.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\NDMaterial.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\nD\TclModelBuilderNDMaterialCommand.cpp
# End Source File
# End Group
# Begin Group "section"

# PROP Default_Filter ""
# Begin Group "repres"

# PROP Default_Filter ""
# Begin Group "cell"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\cell\Cell.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\cell\Cell.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\cell\QuadCell.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\cell\QuadCell.h
# End Source File
# End Group
# Begin Group "patch"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\patch\CircPatch.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\patch\CircPatch.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\patch\Patch.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\patch\Patch.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\patch\QuadPatch.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\patch\QuadPatch.h
# End Source File
# End Group
# Begin Group "reinfBar"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\reinfBar\ReinfBar.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\reinfBar\ReinfBar.h
# End Source File
# End Group
# Begin Group "reinfLayer"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\reinfLayer\CircReinfLayer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\reinfLayer\CircReinfLayer.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\reinfLayer\ReinfLayer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\reinfLayer\ReinfLayer.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\reinfLayer\StraightReinfLayer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\reinfLayer\StraightReinfLayer.h
# End Source File
# End Group
# Begin Group "sect"

# PROP Default_Filter ".cpp; .h"
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\section\FiberSectionRepr.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\section\FiberSectionRepr.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\section\SectionRepres.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\repres\section\SectionRepres.h
# End Source File
# End Group
# End Group
# Begin Group "fiber"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\material\section\fiber\Fiber.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\fiber\Fiber.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\fiber\UniaxialFiber2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\fiber\UniaxialFiber2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\fiber\UniaxialFiber3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\fiber\UniaxialFiber3d.h
# End Source File
# End Group
# Begin Source File

SOURCE=..\..\..\SRC\material\section\ElasticSection2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\ElasticSection2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\ElasticSection3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\ElasticSection3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\FiberSection.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\FiberSection.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\GenericSection1d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\GenericSection1d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\GenericSectionNd.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\GenericSectionNd.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\SectionAggregator.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\SectionAggregator.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\SectionForceDeformation.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\SectionForceDeformation.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\section\TclModelBuilderSectionCommand.cpp
# End Source File
# End Group
# Begin Source File

SOURCE=..\..\..\SRC\material\Material.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\material\Material.h
# End Source File
# End Target
# End Project
