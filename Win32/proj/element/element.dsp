# Microsoft Developer Studio Project File - Name="element" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=element - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "element.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "element.mak" CFG="element - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "element - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "element - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "element - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib\release"
# PROP Intermediate_Dir "..\..\obj\element\release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MLd /W3 /GX /O2 /I "c:\Program Files\tcl\include" /I "..\..\..\src\element\dispBeamColumn" /I "..\..\..\src\element\brick" /I "..\..\..\src\element\shell" /I "..\..\..\src\element\8nbrick" /I "..\..\..\src\recorder\response" /I "..\..\..\src\element\nonlinearBeamColumn\quadRule" /I "..\..\..\src\material\nD" /I "..\..\..\src\element\fourNodeQuad" /I "..\..\..\src\element\damper" /I "..\..\..\src\coordTransformation" /I "..\..\..\src\element\beamWithHinges" /I "..\..\..\src\element\nonlinearBeamColumn\matrixutil" /I "..\..\..\src\element\zeroLength" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\element\feap" /I "..\..\..\src\handler" /I "..\..\..\src\element" /I "..\..\..\src\element\truss" /I "..\..\..\src\material\section" /I "..\..\..\src\element\beam3d" /I "..\..\..\src\element\beam2d" /I "..\..\..\src\material" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\matrix" /I "..\..\..\src\domain\load" /I "..\..\..\src\renderer" /I "..\..\..\src\actor\channel" /I "..\..\..\src\domain\node" /I "..\..\..\src\actor\actor" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src" /I "..\..\..\src\domain\domain" /I "..\..\..\src\material\nd\template3dep" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "element - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib\debug"
# PROP Intermediate_Dir "..\..\obj\element\debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /ZI /Od /I "..\..\..\src\element\material\nd\template3dep" /I "c:\Program Files\tcl\include" /I "..\..\..\src\element\dispBeamColumn" /I "..\..\..\src\element\brick" /I "..\..\..\src\element\shell" /I "..\..\..\src\element\8nbrick" /I "..\..\..\src\recorder\response" /I "..\..\..\src\element\nonlinearBeamColumn\quadRule" /I "..\..\..\src\material\nD" /I "..\..\..\src\element\fourNodeQuad" /I "..\..\..\src\element\damper" /I "..\..\..\src\coordTransformation" /I "..\..\..\src\element\beamWithHinges" /I "..\..\..\src\element\nonlinearBeamColumn\matrixutil" /I "..\..\..\src\element\zeroLength" /I "..\..\..\src\modelbuilder" /I "..\..\..\src\modelbuilder\tcl" /I "..\..\..\src\element\feap" /I "..\..\..\src\handler" /I "..\..\..\src\element" /I "..\..\..\src\element\truss" /I "..\..\..\src\material\section" /I "..\..\..\src\element\beam3d" /I "..\..\..\src\element\beam2d" /I "..\..\..\src\material" /I "..\..\..\src\material\uniaxial" /I "..\..\..\src\actor\objectBroker" /I "..\..\..\src\matrix" /I "..\..\..\src\domain\load" /I "..\..\..\src\renderer" /I "..\..\..\src\actor\channel" /I "..\..\..\src\domain\node" /I "..\..\..\src\actor\actor" /I "..\..\..\src\tagged" /I "..\..\..\src\domain\component" /I "..\..\..\src" /I "..\..\..\src\domain\domain" /I "..\..\..\src\material\nd\template3dep" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /FR /FD /GZ /c
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

# Name "element - Win32 Release"
# Name "element - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\SRC\element\Element.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\ElementalLoad.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\Information.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\TclElementCommands.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\SRC\element\Element.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\ElementalLoad.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\Information.h
# End Source File
# End Group
# Begin Group "truss"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\truss\CorotTruss.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\truss\CorotTruss.h
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\element\truss\CorotTrussSection.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\Src\element\truss\CorotTrussSection.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\truss\TclTrussCommand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\truss\Truss.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\truss\Truss.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\truss\TrussSection.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\truss\TrussSection.h
# End Source File
# End Group
# Begin Group "beam2d"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\beam2d\beam2d02.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beam2d\beam2d02.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beam2d\beam2d03.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beam2d\beam2d03.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beam2d\beam2d04.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beam2d\beam2d04.h
# End Source File
# End Group
# Begin Group "beam3d"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\beam3d\beam3d01.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beam3d\beam3d01.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beam3d\beam3d02.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beam3d\beam3d02.h
# End Source File
# End Group
# Begin Group "feap"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\feap\fElement.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\feap\fElement.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\feap\fElmt02.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\feap\fElmt02.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\feap\TclFeapElementCommand.cpp
# End Source File
# End Group
# Begin Group "zeroLength"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\zeroLength\TclZeroLength.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\zeroLength\ZeroLength.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\zeroLength\ZeroLength.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\zeroLength\ZeroLengthND.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\zeroLength\ZeroLengthND.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\zeroLength\ZeroLengthSection.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\zeroLength\ZeroLengthSection.h
# End Source File
# End Group
# Begin Group "beamWithHinges"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\beamWithHinges\BeamWithHinges2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beamWithHinges\BeamWithHinges2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beamWithHinges\BeamWithHinges3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beamWithHinges\BeamWithHinges3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beamWithHinges\TclBeamWithHingesBuilder.cpp
# End Source File
# End Group
# Begin Group "crdTransf"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\CrdTransf.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\CrdTransf.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\CrdTransf2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\CrdTransf2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\CrdTransf3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\CrdTransf3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\LinearCrdTransf2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\LinearCrdTransf2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\LinearCrdTransf3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\LinearCrdTransf3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\PDeltaCrdTransf2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\PDeltaCrdTransf2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\PDeltaCrdTransf3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\coordTransformation\PDeltaCrdTransf3d.h
# End Source File
# End Group
# Begin Group "fourNodeQuad"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\fourNodeQuad\ConstantPressureVolumeQuad.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\fourNodeQuad\ConstantPressureVolumeQuad.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\fourNodeQuad\EnhancedQuad.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\fourNodeQuad\EnhancedQuad.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\fourNodeQuad\FourNodeQuad.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\fourNodeQuad\FourNodeQuad.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\fourNodeQuad\TclFourNodeQuadCommand.cpp
# End Source File
# End Group
# Begin Group "elasticBeamColumn"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\beam2d\ElasticBeam2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beam2d\ElasticBeam2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beam3d\ElasticBeam3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\beam3d\ElasticBeam3d.h
# End Source File
# End Group
# Begin Group "eightNodeBrick"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\brick\BbarBrick.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\brick\BbarBrick.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\brick\Brick.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\brick\Brick.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\8nbrick\EightNodeBrick.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\8nbrick\EightNodeBrick.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\brick\shp3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\brick\TclBrickCommand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\8nbrick\TclEightNodeBrickCommand.cpp
# End Source File
# End Group
# Begin Group "shell"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\shell\R3vectors.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\shell\R3vectors.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\shell\ShellMITC4.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\shell\ShellMITC4.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\shell\TclShellCommand.cpp
# End Source File
# End Group
# Begin Group "dispBeamColumn"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\SRC\element\dispBeamColumn\DispBeamColumn2d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\dispBeamColumn\DispBeamColumn2d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\dispBeamColumn\DispBeamColumn3d.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\dispBeamColumn\DispBeamColumn3d.h
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\dispBeamColumn\TclDispBeamColumnCommand.cpp
# End Source File
# End Group
# Begin Source File

SOURCE=..\..\..\SRC\element\NewElement.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\SRC\element\NewElement.h
# End Source File
# End Target
# End Project
