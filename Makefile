# =============================================================================
# SnapWave Makefile (clean build layout)
#
# Outputs:
#   ./snapwave/build/bin/snapwave   (final executable)
#   ./snapwave/build/obj/...        (objects; removed if STRIP_BUILD=1)
#   ./snapwave/build/mod/...        (modules; removed if STRIP_BUILD=1)
#
# Usage:
#   make
#   make DEBUG=1
#   make STRIP_BUILD=1
#   make clean
# =============================================================================

# -----------------------------
# Toolchain
# -----------------------------
FC        := gfortran     # or ifort, gfortran, ...
CC        := gcc     # or gcc/clang, ...
NF_CONFIG := nf-config

# -----------------------------
# Output directories (everything under snapwave/build)
# -----------------------------
BUILD_ROOT := ./snapwave/lnx64
OBJ_DIR    := $(BUILD_ROOT)/obj
MOD_DIR    := $(BUILD_ROOT)/mod
BIN_DIR    := $(BUILD_ROOT)/bin
TARGET     := $(BIN_DIR)/snapwave

# -----------------------------
# NetCDF config
# -----------------------------
NC_FFLAGS := $(shell $(NF_CONFIG) --fflags 2>/dev/null || echo "")
NC_FLIBS  := $(shell $(NF_CONFIG) --flibs  2>/dev/null || echo "-lnetcdff -lnetcdf")

# -----------------------------
# Compiler family detection
# -----------------------------
FC_BASENAME := $(notdir $(FC))
IS_INTEL    := $(filter ifx ifort,$(FC_BASENAME))

# -----------------------------
# Common flags
# -----------------------------
# OpenMP + module output flag differs by compiler
ifeq ($(strip $(IS_INTEL)),)
  OMPFLAG := -fopenmp
  MODFLAG := -J$(MOD_DIR)
  PPFLAG  := -cpp
  LINEFLAG := -ffree-line-length-none
else
  OMPFLAG := -qopenmp
  MODFLAG := -module $(MOD_DIR)
  PPFLAG  := -fpp
  LINEFLAG :=
endif

# C flags (triangle needs TRILIBRARY)
CFLAGS ?= -O2 -DTRILIBRARY

# Fortran flags: base
FFLAGS_BASE := $(OMPFLAG) $(MODFLAG) -I$(MOD_DIR) $(NC_FFLAGS) $(PPFLAG) $(LINEFLAG)

# Debug vs release flags (avoid GNU-only flags on Intel)
ifeq ($(DEBUG),1)
  ifeq ($(strip $(IS_INTEL)),)
    FFLAGS ?= -g -O0 -fcheck=all -fbacktrace $(FFLAGS_BASE)
  else
    # Intel-friendly debug flags
    FFLAGS ?= -g -O0 -traceback -check all $(FFLAGS_BASE)
  endif
  CFLAGS += -g -O0
else
  FFLAGS ?= -O2 $(FFLAGS_BASE)
endif

# -----------------------------
# Sources / objects
# -----------------------------

# 1) Third party
TRIANGLE_SRC := third_party_open/triangle/triangle.c third_party_open/triangle/tricall2.c
TRIANGLE_OBJ := $(OBJ_DIR)/triangle.o $(OBJ_DIR)/tricall2.o

KDTREE_SRCS := $(shell find third_party_open/kdtree2 -type f \( -iname "*.f90" -o -iname "*.F90" \) ! -iname "*test*" ! -iname "*main*" -print | sort)
KDTREE_OBJS := $(patsubst third_party_open/kdtree2/%.f90,$(OBJ_DIR)/kdtree_%.o,$(KDTREE_SRCS))
KDTREE_OBJS := $(patsubst third_party_open/kdtree2/%.F90,$(OBJ_DIR)/kdtree_%.o,$(KDTREE_OBJS))

# 2) utils_lgpl (explicit order)
UTILS_FILES_ORDER := \
    utils_lgpl/deltares_common/src/deltares_common_modules.f90 \
    utils_lgpl/deltares_common/src/malloc.f90 \
    utils_lgpl/deltares_common/src/m_ec_triangle.f90 \
    utils_lgpl/kdtree_wrapper/src/kdtreeWrapper.f90

UTILS_OBJS := $(patsubst utils_lgpl/%.f90,$(OBJ_DIR)/utils_lgpl/%.o,$(UTILS_FILES_ORDER))

# 3) src (explicit order; module ordering matters)
SRC_FILES_ORDER := \
    src/snapwave_data.f90 \
    src/snapwave_date.f90 \
    src/snapwave_results.f90 \
    src/interp.F90 \
    src/snapwave_input.f90 \
    src/snapwave_windsource.f90 \
    src/snapwave_ncoutput.F90 \
    src/snapwave_domain.f90 \
    src/snapwave_boundaries.f90 \
    src/snapwave_obspoints.f90 \
    src/snapwave_solver.f90 \
    src/snapwave.f90

SRC_OBJS := $(patsubst src/%.f90,$(OBJ_DIR)/%.o,$(filter %.f90,$(SRC_FILES_ORDER))) \
            $(patsubst src/%.F90,$(OBJ_DIR)/%.o,$(filter %.F90,$(SRC_FILES_ORDER)))

ALL_OBJS := $(TRIANGLE_OBJ) $(KDTREE_OBJS) $(UTILS_OBJS) $(SRC_OBJS)

# -----------------------------
# Targets
# -----------------------------
.PHONY: all clean info directories
.DELETE_ON_ERROR:

all: directories $(TARGET)

directories:
	@mkdir -p $(BIN_DIR) $(OBJ_DIR) $(MOD_DIR)

info:
	@echo "FC:        $(FC)"
	@echo "CC:        $(CC)"
	@echo "DEBUG:     $(DEBUG)"
	@echo "STRIP_BUILD: $(STRIP_BUILD)"
	@echo "OBJ_DIR:   $(OBJ_DIR)"
	@echo "MOD_DIR:   $(MOD_DIR)"
	@echo "TARGET:    $(TARGET)"

# Link
$(TARGET): $(ALL_OBJS) | directories
	@echo "--> Linking executable: $@"
	$(FC) $(FFLAGS) -o $@ $(ALL_OBJS) $(NC_FLIBS)
	@if [ "$(STRIP_BUILD)" = "1" ]; then \
	  echo "--> STRIP_BUILD=1: removing intermediate obj/mod (keeping $@)"; \
	  rm -rf "$(OBJ_DIR)" "$(MOD_DIR)"; \
	fi

# -----------------------------
# Build rules
# -----------------------------

# Triangle (C)
$(OBJ_DIR)/triangle.o: third_party_open/triangle/triangle.c | directories
	@echo "Compiling C (Triangle): $<"
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/tricall2.o: third_party_open/triangle/tricall2.c | directories
	@echo "Compiling C (Tricall2): $<"
	$(CC) $(CFLAGS) -c $< -o $@

# Kdtree (Fortran)
$(OBJ_DIR)/kdtree_%.o: third_party_open/kdtree2/%.f90 | directories
	@mkdir -p $(dir $@)
	@echo "Compiling Fortran (Kdtree): $<"
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/kdtree_%.o: third_party_open/kdtree2/%.F90 | directories
	@mkdir -p $(dir $@)
	@echo "Compiling Fortran (Kdtree): $<"
	$(FC) $(FFLAGS) -c $< -o $@

# Utils
$(OBJ_DIR)/utils_lgpl/%.o: utils_lgpl/%.f90 | directories
	@mkdir -p $(dir $@)
	@echo "Compiling Utils: $<"
	$(FC) $(FFLAGS) -c $< -o $@

# SRC (Main)
$(OBJ_DIR)/%.o: src/%.f90 | directories
	@mkdir -p $(dir $@)
	@echo "Compiling SRC: $<"
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: src/%.F90 | directories
	@mkdir -p $(dir $@)
	@echo "Compiling SRC: $<"
	$(FC) $(FFLAGS) -c $< -o $@

# -----------------------------
# Explicit dependencies (module order enforcement)
# -----------------------------
$(UTILS_OBJS): $(KDTREE_OBJS) $(TRIANGLE_OBJ)
$(SRC_OBJS): $(UTILS_OBJS)

# Internal SRC deps
$(OBJ_DIR)/snapwave_input.o:      $(OBJ_DIR)/snapwave_data.o $(OBJ_DIR)/snapwave_date.o
$(OBJ_DIR)/snapwave_windsource.o: $(OBJ_DIR)/snapwave_data.o
$(OBJ_DIR)/snapwave_ncoutput.o:   $(OBJ_DIR)/snapwave_results.o $(OBJ_DIR)/snapwave_date.o $(OBJ_DIR)/snapwave_data.o
$(OBJ_DIR)/snapwave_domain.o:     $(OBJ_DIR)/snapwave_ncoutput.o $(OBJ_DIR)/snapwave_input.o $(OBJ_DIR)/interp.o \
                                  $(OBJ_DIR)/snapwave_results.o $(OBJ_DIR)/snapwave_data.o
$(OBJ_DIR)/snapwave_boundaries.o: $(OBJ_DIR)/snapwave_domain.o $(OBJ_DIR)/snapwave_data.o
$(OBJ_DIR)/snapwave_obspoints.o:  $(OBJ_DIR)/snapwave_data.o $(OBJ_DIR)/interp.o
$(OBJ_DIR)/snapwave_solver.o:     $(OBJ_DIR)/snapwave_domain.o $(OBJ_DIR)/snapwave_windsource.o $(OBJ_DIR)/snapwave_ncoutput.o \
                                  $(OBJ_DIR)/snapwave_data.o
$(OBJ_DIR)/snapwave.o:            $(OBJ_DIR)/snapwave_solver.o $(OBJ_DIR)/snapwave_obspoints.o $(OBJ_DIR)/snapwave_boundaries.o

# -----------------------------
# Clean
# -----------------------------
clean:
	rm -rf $(BUILD_ROOT)
