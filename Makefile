# Using NVIDIA HPC Compiler nvc++
CXX = nvc++
CXXFLAGS = -O3 -std=c++20

CC = nvc++
CCFLAGS = -O3 -std=c11

NVHPC_OMP = -mp=multicore
NVHPC_STDPAR = -stdpar=gpu
NVHPC_STDEXEC = -stdpar=gpu --experimental-stdpar

# C++ source files
SRCS = $(wildcard ./Convergence/*.c*) \
       $(wildcard ./Normal/*.c*)

# Create names of executable binaries from .cpp or .c file
TARGETS = $(basename $(SRCS))

# Default targets
all: $(TARGETS)

# The rules for each source files
%: %.cpp
	@if echo $< | grep -q -- "-omp"; then \
		echo "Compiling $< with NVHPC_OMP flags"; \
		$(CXX) $(CXXFLAGS) $(NVHPC_OMP) $< -o $(notdir $@); \
	elif echo $< | grep -q -- "-stdpar"; then \
		echo "Compiling $< with NVHPC_STDPAR flags"; \
		$(CXX) $(CXXFLAGS) $(NVHPC_STDPAR) $< -o $(notdir $@); \
	elif echo $< | grep -q -- "-sender"; then \
		echo "Compiling $< with NVHPC_STDEXEC flags"; \
		$(CXX) $(CXXFLAGS) $(NVHPC_STDEXEC) $< -o $(notdir $@); \
	else \
		echo "Compiling $< with default C++ flags"; \
		$(CXX) $(CXXFLAGS) $< -o $(notdir $@); \
	fi
%: %.c
	@echo "Compiling $< with default C flags"; \
	$(CC) $(CCFLAGS) $< -o $(notdir $@);
	

# Clean up
.PHONY: clean
clean:
	rm -f $(notdir $(TARGETS))

# Help
.PHONY: help
help:
	@echo "Makefile for himenoBMT-par C++ Project"
	@echo "Usage:"
	@echo "  make        - build all binaries"
	@echo "  make clean  - remove generated files"
	@echo "  make help   - display this help message"
