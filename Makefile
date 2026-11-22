CXX ?= g++
PYTHON ?= python3
EIGEN_INCLUDE ?= /usr/include/eigen3
CPPFLAGS += -Iinclude -I"$(EIGEN_INCLUDE)"
CXXFLAGS ?= -O2
CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic
BUILD := build
PROGRAMS := initial-incommensurate continue-incommensurate initial-commensurate continue-commensurate
HEADERS := $(wildcard include/langevin/*.hpp)

.PHONY: all test examples
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

$(BUILD)/test-incommensurate: tests/numerics.cpp $(HEADERS) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@
$(BUILD)/test-commensurate: tests/numerics.cpp $(HEADERS) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DCOMMENSURATE $< -o $@

$(BUILD)/test-headers: tests/headers_main.cpp tests/headers_other.cpp $(HEADERS) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) tests/headers_main.cpp tests/headers_other.cpp -o $@

$(BUILD)/library-example: examples/library_usage.cpp $(HEADERS) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

examples: $(BUILD)/library-example

test: all $(BUILD)/test-incommensurate $(BUILD)/test-commensurate $(BUILD)/test-headers $(BUILD)/library-example
	$(BUILD)/test-incommensurate
	$(BUILD)/test-commensurate
	LANGEVIN_BUILD_DIR=$(BUILD) $(PYTHON) tests/test_cli.py
