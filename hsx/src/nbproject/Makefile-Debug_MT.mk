#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_CONF=Debug_MT
CND_DISTDIR=dist

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/2093546553/atom_GN.o \
	${OBJECTDIR}/_ext/2093546553/protein_GN.o \
	${OBJECTDIR}/_ext/2093546553/crystallizer_GN.o \
	${OBJECTDIR}/_ext/2093546553/atom_properties_GN.o \
	${OBJECTDIR}/_ext/2093546553/structure_additional_GN.o \
	${OBJECTDIR}/_ext/2093546553/molecule_GN.o \
	${OBJECTDIR}/_ext/2093546553/files_GN.o \
	${OBJECTDIR}/_ext/2093546553/delaunay_GN.o \
	${OBJECTDIR}/hotspotsX_mt.o \
	${OBJECTDIR}/_ext/2093546553/structure_GN.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-static -pedantic -Wall -m32 -D_LINUX_OS -pthread
CXXFLAGS=-static -pedantic -Wall -m32 -D_LINUX_OS -pthread

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-Debug_MT.mk dist/Debug_MT/GNU-Linux-x86/hotspotsx_mt

dist/Debug_MT/GNU-Linux-x86/hotspotsx_mt: ${OBJECTFILES}
	${MKDIR} -p dist/Debug_MT/GNU-Linux-x86
	${LINK.cc} -pthread -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hotspotsx_mt ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/2093546553/atom_GN.o: ../../C_libs/atom_GN.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2093546553
	${RM} $@.d
	$(COMPILE.cc) -g -I../../C_libs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2093546553/atom_GN.o ../../C_libs/atom_GN.cpp

${OBJECTDIR}/_ext/2093546553/protein_GN.o: ../../C_libs/protein_GN.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2093546553
	${RM} $@.d
	$(COMPILE.cc) -g -I../../C_libs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2093546553/protein_GN.o ../../C_libs/protein_GN.cpp

${OBJECTDIR}/_ext/2093546553/crystallizer_GN.o: ../../C_libs/crystallizer_GN.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2093546553
	${RM} $@.d
	$(COMPILE.cc) -g -I../../C_libs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2093546553/crystallizer_GN.o ../../C_libs/crystallizer_GN.cpp

${OBJECTDIR}/_ext/2093546553/atom_properties_GN.o: ../../C_libs/atom_properties_GN.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2093546553
	${RM} $@.d
	$(COMPILE.cc) -g -I../../C_libs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2093546553/atom_properties_GN.o ../../C_libs/atom_properties_GN.cpp

${OBJECTDIR}/_ext/2093546553/structure_additional_GN.o: ../../C_libs/structure_additional_GN.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2093546553
	${RM} $@.d
	$(COMPILE.cc) -g -I../../C_libs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2093546553/structure_additional_GN.o ../../C_libs/structure_additional_GN.cpp

${OBJECTDIR}/_ext/2093546553/molecule_GN.o: ../../C_libs/molecule_GN.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2093546553
	${RM} $@.d
	$(COMPILE.cc) -g -I../../C_libs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2093546553/molecule_GN.o ../../C_libs/molecule_GN.cpp

${OBJECTDIR}/_ext/2093546553/files_GN.o: ../../C_libs/files_GN.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2093546553
	${RM} $@.d
	$(COMPILE.cc) -g -I../../C_libs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2093546553/files_GN.o ../../C_libs/files_GN.cpp

${OBJECTDIR}/_ext/2093546553/delaunay_GN.o: ../../C_libs/delaunay_GN.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2093546553
	${RM} $@.d
	$(COMPILE.cc) -g -I../../C_libs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2093546553/delaunay_GN.o ../../C_libs/delaunay_GN.cpp

${OBJECTDIR}/hotspotsX_mt.o: hotspotsX_mt.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../../C_libs -MMD -MP -MF $@.d -o ${OBJECTDIR}/hotspotsX_mt.o hotspotsX_mt.cpp

${OBJECTDIR}/_ext/2093546553/structure_GN.o: ../../C_libs/structure_GN.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/2093546553
	${RM} $@.d
	$(COMPILE.cc) -g -I../../C_libs -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/2093546553/structure_GN.o ../../C_libs/structure_GN.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug_MT
	${RM} dist/Debug_MT/GNU-Linux-x86/hotspotsx_mt

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
