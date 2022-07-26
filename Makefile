CXX = g++
INCLUDE_DIR = -I./include -I./gflib/include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = rdeval
BUILD = build/bin
SOURCE = src
INCLUDE = include
LDFLAGS :=

GFLIB_DIR := $(CURDIR)/gflib
GFLIB_FILES := $(GFLIB_DIR)/$(SOURCE)/* $(GFLIB_DIR)/$(INCLUDE)/*

main: $(SOURCE)/main.cpp $(GFLIB_FILES) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SOURCE)/main.cpp -o $(BUILD)/$(TARGET)

$(GFLIB_FILES): gflib
	@# Do nothing


.PHONY: gflib
gflib:
	$(MAKE) -j -C $(GFLIB_DIR)
	
$(BUILD):
	-mkdir -p $@
	
clean:
	$(MAKE) -j -C $(GFLIB_DIR) clean
	$(RM) -r build
