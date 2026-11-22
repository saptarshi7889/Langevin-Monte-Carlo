CXX ?= g++
PYTHON ?= python3
EIGEN_INCLUDE ?= /usr/include/eigen3
CPPFLAGS += -Iinclude -I"$(EIGEN_INCLUDE)"
CXXFLAGS ?= -O2
CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic
BUILD := build
PROGRAMS := initial-incommensurate continue-incommensurate initial-commensurate continue-commensurate
HEADERS := include/kernel.hpp include/simulation.hpp

.PHONY: all test
all: $(addprefix $(BUILD)/,$(PROGRAMS))

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/initial-incommensurate: src/incommensurate/Langevin_simulation_initial_simulation.cpp $(HEADERS) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@
$(BUILD)/continue-incommensurate: src/incommensurate/Langevin_simulation_continuted_simulation.cpp $(HEADERS) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@
$(BUILD)/initial-commensurate: src/commensurate/Langevin_simulation_initial_commensurate.cpp $(HEADERS) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@
$(BUILD)/continue-commensurate: src/commensurate/Langevin_simulation_continued_commensurate.cpp $(HEADERS) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

$(BUILD)/test-incommensurate: tests/numerics.cpp src/incommensurate/Langevin_simulation_initial_simulation.cpp $(HEADERS) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@
$(BUILD)/test-commensurate: tests/numerics.cpp src/commensurate/Langevin_simulation_initial_commensurate.cpp $(HEADERS) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DCOMMENSURATE $< -o $@

test: all $(BUILD)/test-incommensurate $(BUILD)/test-commensurate
	$(BUILD)/test-incommensurate
	$(BUILD)/test-commensurate
	LANGEVIN_BUILD_DIR=$(BUILD) $(PYTHON) tests/test_cli.py
