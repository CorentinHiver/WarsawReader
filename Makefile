CXX := g++
STD := -std=c++17
CXXFLAGS := -Wall -Wextra
OPT := -O3 -march=native
DEPFLAGS := -MD -MP

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --glibs)

BUILD_DIR := build
TARGETS := CFDoptimizations writeTraces studyCFD caen2root #rootReaderExample

OBJS := $(addprefix $(BUILD_DIR)/,$(addsuffix .o,$(TARGETS)))
DEPS := $(OBJS:.o=.d)

all: $(TARGETS)
	@echo "Cleanup: Removing build directory..."
	rm -rf $(BUILD_DIR)

# Link
$(TARGETS): %: $(BUILD_DIR)/%.o
	$(CXX) -o $@ $< $(ROOTLIBS)

# Compile with explicit dependency file location
$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) -c $< -o $@ $(ROOTCFLAGS) $(CXXFLAGS) $(STD) $(OPT) \
	$(DEPFLAGS) -MF $(BUILD_DIR)/$*.d -MT $@

# Ensure directory exists
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Include dependencies if they exist
-include $(DEPS)

debug: OPT := -Og
debug: CXXFLAGS += -g -fno-omit-frame-pointer
debug: all

clean:
	rm -rf $(BUILD_DIR) $(TARGETS)