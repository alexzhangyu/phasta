# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/cmake/3.0.0/bin/cmake

# The command to remove a file.
RM = /usr/local/cmake/3.0.0/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /usr/local/cmake/3.0.0/share/cmake-3.0/Modules/FortranCInterface

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/CMakeFiles/FortranCInterface

# Include any dependencies generated for this target.
include CMakeFiles/myfort.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/myfort.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/myfort.dir/flags.make

CMakeFiles/myfort.dir/mysub.f.o: CMakeFiles/myfort.dir/flags.make
CMakeFiles/myfort.dir/mysub.f.o: /usr/local/cmake/3.0.0/share/cmake-3.0/Modules/FortranCInterface/mysub.f
	$(CMAKE_COMMAND) -E cmake_progress_report /fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/CMakeFiles/FortranCInterface/CMakeFiles $(CMAKE_PROGRESS_1)
	@echo "Building Fortran object CMakeFiles/myfort.dir/mysub.f.o"
	/usr/local/gcc/4.9.2/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /usr/local/cmake/3.0.0/share/cmake-3.0/Modules/FortranCInterface/mysub.f -o CMakeFiles/myfort.dir/mysub.f.o

CMakeFiles/myfort.dir/mysub.f.o.requires:
.PHONY : CMakeFiles/myfort.dir/mysub.f.o.requires

CMakeFiles/myfort.dir/mysub.f.o.provides: CMakeFiles/myfort.dir/mysub.f.o.requires
	$(MAKE) -f CMakeFiles/myfort.dir/build.make CMakeFiles/myfort.dir/mysub.f.o.provides.build
.PHONY : CMakeFiles/myfort.dir/mysub.f.o.provides

CMakeFiles/myfort.dir/mysub.f.o.provides.build: CMakeFiles/myfort.dir/mysub.f.o

CMakeFiles/myfort.dir/my_sub.f.o: CMakeFiles/myfort.dir/flags.make
CMakeFiles/myfort.dir/my_sub.f.o: /usr/local/cmake/3.0.0/share/cmake-3.0/Modules/FortranCInterface/my_sub.f
	$(CMAKE_COMMAND) -E cmake_progress_report /fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/CMakeFiles/FortranCInterface/CMakeFiles $(CMAKE_PROGRESS_2)
	@echo "Building Fortran object CMakeFiles/myfort.dir/my_sub.f.o"
	/usr/local/gcc/4.9.2/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /usr/local/cmake/3.0.0/share/cmake-3.0/Modules/FortranCInterface/my_sub.f -o CMakeFiles/myfort.dir/my_sub.f.o

CMakeFiles/myfort.dir/my_sub.f.o.requires:
.PHONY : CMakeFiles/myfort.dir/my_sub.f.o.requires

CMakeFiles/myfort.dir/my_sub.f.o.provides: CMakeFiles/myfort.dir/my_sub.f.o.requires
	$(MAKE) -f CMakeFiles/myfort.dir/build.make CMakeFiles/myfort.dir/my_sub.f.o.provides.build
.PHONY : CMakeFiles/myfort.dir/my_sub.f.o.provides

CMakeFiles/myfort.dir/my_sub.f.o.provides.build: CMakeFiles/myfort.dir/my_sub.f.o

CMakeFiles/myfort.dir/mymodule.f90.o: CMakeFiles/myfort.dir/flags.make
CMakeFiles/myfort.dir/mymodule.f90.o: /usr/local/cmake/3.0.0/share/cmake-3.0/Modules/FortranCInterface/mymodule.f90
	$(CMAKE_COMMAND) -E cmake_progress_report /fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/CMakeFiles/FortranCInterface/CMakeFiles $(CMAKE_PROGRESS_3)
	@echo "Building Fortran object CMakeFiles/myfort.dir/mymodule.f90.o"
	/usr/local/gcc/4.9.2/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /usr/local/cmake/3.0.0/share/cmake-3.0/Modules/FortranCInterface/mymodule.f90 -o CMakeFiles/myfort.dir/mymodule.f90.o

CMakeFiles/myfort.dir/mymodule.f90.o.requires:
.PHONY : CMakeFiles/myfort.dir/mymodule.f90.o.requires

CMakeFiles/myfort.dir/mymodule.f90.o.provides: CMakeFiles/myfort.dir/mymodule.f90.o.requires
	$(MAKE) -f CMakeFiles/myfort.dir/build.make CMakeFiles/myfort.dir/mymodule.f90.o.provides.build
.PHONY : CMakeFiles/myfort.dir/mymodule.f90.o.provides

CMakeFiles/myfort.dir/mymodule.f90.o.provides.build: CMakeFiles/myfort.dir/mymodule.f90.o

CMakeFiles/myfort.dir/my_module.f90.o: CMakeFiles/myfort.dir/flags.make
CMakeFiles/myfort.dir/my_module.f90.o: /usr/local/cmake/3.0.0/share/cmake-3.0/Modules/FortranCInterface/my_module.f90
	$(CMAKE_COMMAND) -E cmake_progress_report /fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/CMakeFiles/FortranCInterface/CMakeFiles $(CMAKE_PROGRESS_4)
	@echo "Building Fortran object CMakeFiles/myfort.dir/my_module.f90.o"
	/usr/local/gcc/4.9.2/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /usr/local/cmake/3.0.0/share/cmake-3.0/Modules/FortranCInterface/my_module.f90 -o CMakeFiles/myfort.dir/my_module.f90.o

CMakeFiles/myfort.dir/my_module.f90.o.requires:
.PHONY : CMakeFiles/myfort.dir/my_module.f90.o.requires

CMakeFiles/myfort.dir/my_module.f90.o.provides: CMakeFiles/myfort.dir/my_module.f90.o.requires
	$(MAKE) -f CMakeFiles/myfort.dir/build.make CMakeFiles/myfort.dir/my_module.f90.o.provides.build
.PHONY : CMakeFiles/myfort.dir/my_module.f90.o.provides

CMakeFiles/myfort.dir/my_module.f90.o.provides.build: CMakeFiles/myfort.dir/my_module.f90.o

# Object files for target myfort
myfort_OBJECTS = \
"CMakeFiles/myfort.dir/mysub.f.o" \
"CMakeFiles/myfort.dir/my_sub.f.o" \
"CMakeFiles/myfort.dir/mymodule.f90.o" \
"CMakeFiles/myfort.dir/my_module.f90.o"

# External object files for target myfort
myfort_EXTERNAL_OBJECTS =

libmyfort.a: CMakeFiles/myfort.dir/mysub.f.o
libmyfort.a: CMakeFiles/myfort.dir/my_sub.f.o
libmyfort.a: CMakeFiles/myfort.dir/mymodule.f90.o
libmyfort.a: CMakeFiles/myfort.dir/my_module.f90.o
libmyfort.a: CMakeFiles/myfort.dir/build.make
libmyfort.a: CMakeFiles/myfort.dir/link.txt
	@echo "Linking Fortran static library libmyfort.a"
	$(CMAKE_COMMAND) -P CMakeFiles/myfort.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/myfort.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/myfort.dir/build: libmyfort.a
.PHONY : CMakeFiles/myfort.dir/build

CMakeFiles/myfort.dir/requires: CMakeFiles/myfort.dir/mysub.f.o.requires
CMakeFiles/myfort.dir/requires: CMakeFiles/myfort.dir/my_sub.f.o.requires
CMakeFiles/myfort.dir/requires: CMakeFiles/myfort.dir/mymodule.f90.o.requires
CMakeFiles/myfort.dir/requires: CMakeFiles/myfort.dir/my_module.f90.o.requires
.PHONY : CMakeFiles/myfort.dir/requires

CMakeFiles/myfort.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/myfort.dir/cmake_clean.cmake
.PHONY : CMakeFiles/myfort.dir/clean

CMakeFiles/myfort.dir/depend:
	cd /fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/CMakeFiles/FortranCInterface && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /usr/local/cmake/3.0.0/share/cmake-3.0/Modules/FortranCInterface /usr/local/cmake/3.0.0/share/cmake-3.0/Modules/FortranCInterface /fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/CMakeFiles/FortranCInterface /fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/CMakeFiles/FortranCInterface /fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/CMakeFiles/FortranCInterface/CMakeFiles/myfort.dir/DependInfo.cmake
.PHONY : CMakeFiles/myfort.dir/depend

