# =============================================================================
# COMPILER CONFIGURATION AND FLAGS
# =============================================================================
FC = ifort
CC = icx # Or gcc, if you prefer
NF_CONFIG = nf-config

# Output directories
BUILD_DIR := build
BIN_DIR := bin
TARGET := $(BIN_DIR)/snapwave

# NetCDF Configuration
NC_FFLAGS := $(shell $(NF_CONFIG) --fflags 2>/dev/null || echo "")
NC_FLIBS  := $(shell $(NF_CONFIG) --flibs  2>/dev/null || echo "-lnetcdff -lnetcdf")

# Compilation Flags
# Detect if ifort or gfortran is used for the module flag
ifeq ($(findstring ifort,$(FC)),ifort)
  OMPFLAG := -qopenmp
  MODFLAG := -module $(BUILD_DIR)
  # -DTRILIBRARY is necessary for triangle to compile as a library and not look for main()
  CFLAGS  := -O2 -DTRILIBRARY
else
  OMPFLAG := -fopenmp
  MODFLAG := -J$(BUILD_DIR)
  CFLAGS  := -O2 -DTRILIBRARY
endif

FFLAGS := -O2 $(OMPFLAG) $(MODFLAG) -I$(BUILD_DIR) $(NC_FFLAGS)

# =============================================================================
# SOURCE AND OBJECT DEFINITIONS
# =============================================================================

# 1. THIRD PARTY (Triangle in C and Kdtree in Fortran)
# -------------------------------------------------
TRIANGLE_SRC := third_party_open/triangle/triangle.c third_party_open/triangle/tricall2.c
TRIANGLE_OBJ := $(BUILD_DIR)/triangle.o $(BUILD_DIR)/tricall2.o

# Search for kdtree sources (excluding tests)
KDTREE_SRCS := $(shell find third_party_open/kdtree2 -type f \( -iname "*.f90" -o -iname "*.F90" \) ! -iname "*test*" ! -iname "*main*" -print | sort)
KDTREE_OBJS := $(patsubst third_party_open/kdtree2/%.f90, $(BUILD_DIR)/kdtree_%.o, $(KDTREE_SRCS))
KDTREE_OBJS := $(patsubst third_party_open/kdtree2/%.F90, $(BUILD_DIR)/kdtree_%.o, $(KDTREE_OBJS))

# 2. UTILS (LGPL Utility Libraries)
# -------------------------------------
# Define explicit order for utils_lgpl
UTILS_FILES_ORDER := \
    utils_lgpl/deltares_common/src/deltares_common_modules.f90 \
    utils_lgpl/deltares_common/src/malloc.f90 \
    utils_lgpl/deltares_common/src/m_ec_triangle.f90 \
    utils_lgpl/kdtree_wrapper/src/kdtreeWrapper.f90

# Map sources to objects maintaining directory structure in build/
UTILS_OBJS := $(patsubst utils_lgpl/%.f90, $(BUILD_DIR)/utils_lgpl/%.o, $(UTILS_FILES_ORDER))

# 3. SRC (Main Code in STRICT ORDER)
# -------------------------------------------
# Order matters for module generation
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

# Handle .f90 and .F90 extensions
SRC_OBJS := $(patsubst src/%.f90, $(BUILD_DIR)/%.o, $(filter %.f90,$(SRC_FILES_ORDER))) \
            $(patsubst src/%.F90, $(BUILD_DIR)/%.o, $(filter %.F90,$(SRC_FILES_ORDER)))

# Complete list of objects for final linking
ALL_OBJS := $(TRIANGLE_OBJ) $(KDTREE_OBJS) $(UTILS_OBJS) $(SRC_OBJS)

# =============================================================================
# COMPILATION RULES
# =============================================================================

.PHONY: all clean info directories

all: directories $(TARGET)

directories:
	@mkdir -p $(BIN_DIR) $(BUILD_DIR)

info:
	@echo "Compiler FC: $(FC)"
	@echo "Compiler CC: $(CC)"
	@echo "SRC Objects: $(SRC_OBJS)"

# Final Linking
$(TARGET): $(ALL_OBJS)
	@echo "--> Linking executable: $@"
	$(FC) $(FFLAGS) -o $@ $(ALL_OBJS) $(NC_FLIBS)

# --- Rules for Third Party ---

# Triangle (C)
$(BUILD_DIR)/triangle.o: third_party_open/triangle/triangle.c
	@echo "Compiling C (Triangle): $<"
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/tricall2.o: third_party_open/triangle/tricall2.c
	@echo "Compiling C (Tricall2): $<"
	$(CC) $(CFLAGS) -c $< -o $@

# Kdtree (Fortran)
$(BUILD_DIR)/kdtree_%.o: third_party_open/kdtree2/%.f90
	@mkdir -p $(dir $@)
	@echo "Compiling Fortran (Kdtree): $<"
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/kdtree_%.o: third_party_open/kdtree2/%.F90
	@mkdir -p $(dir $@)
	@echo "Compiling Fortran (Kdtree): $<"
	$(FC) $(FFLAGS) -c $< -o $@

# --- Rules for Utils ---

# Generic rule for utils (maintaining structure)
$(BUILD_DIR)/utils_lgpl/%.o: utils_lgpl/%.f90
	@mkdir -p $(dir $@)
	@echo "Compiling Utils: $<"
	$(FC) $(FFLAGS) -c $< -o $@

# --- Rules for SRC (Main) ---

# Generic rule for src .f90
$(BUILD_DIR)/%.o: src/%.f90
	@echo "Compiling SRC: $<"
	$(FC) $(FFLAGS) -c $< -o $@

# Generic rule for src .F90
$(BUILD_DIR)/%.o: src/%.F90
	@echo "Compiling SRC: $<"
	$(FC) $(FFLAGS) -c $< -o $@

# =============================================================================
# EXPLICIT DEPENDENCIES (To ensure module order)
# =============================================================================

# 1. Utils depends on Third Party (if they use kdtree, for example)
$(UTILS_OBJS): $(KDTREE_OBJS) $(TRIANGLE_OBJ)

# 2. SRC depends on Utils (and therefore on Third Party)
$(SRC_OBJS): $(UTILS_OBJS)

# 3. Internal SRC dependencies
$(BUILD_DIR)/snapwave_input.o:      $(BUILD_DIR)/snapwave_data.o $(BUILD_DIR)/snapwave_date.o
$(BUILD_DIR)/snapwave_windsource.o: $(BUILD_DIR)/snapwave_data.o
$(BUILD_DIR)/snapwave_ncoutput.o:   $(BUILD_DIR)/snapwave_results.o $(BUILD_DIR)/snapwave_date.o $(BUILD_DIR)/snapwave_data.o
$(BUILD_DIR)/snapwave_domain.o:     $(BUILD_DIR)/snapwave_ncoutput.o $(BUILD_DIR)/snapwave_input.o $(BUILD_DIR)/interp.o \
                                    $(BUILD_DIR)/snapwave_results.o $(BUILD_DIR)/snapwave_data.o
$(BUILD_DIR)/snapwave_boundaries.o: $(BUILD_DIR)/snapwave_domain.o $(BUILD_DIR)/snapwave_data.o
$(BUILD_DIR)/snapwave_obspoints.o:  $(BUILD_DIR)/snapwave_data.o $(BUILD_DIR)/interp.o
$(BUILD_DIR)/snapwave_solver.o:     $(BUILD_DIR)/snapwave_domain.o $(BUILD_DIR)/snapwave_windsource.o $(BUILD_DIR)/snapwave_ncoutput.o \
                                    $(BUILD_DIR)/snapwave_data.o
$(BUILD_DIR)/snapwave.o:            $(BUILD_DIR)/snapwave_solver.o $(BUILD_DIR)/snapwave_obspoints.o $(BUILD_DIR)/snapwave_boundaries.o

# =============================================================================
# CLEANING
# =============================================================================

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)
