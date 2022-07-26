CXX = g++
INCLUDE_DIR = -I./include -I./gfalibs/include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = rdeval
BUILD = build/bin
SOURCE = src
INCLUDE = include
BINDIR := $(BUILD)/.o

LDFLAGS := -pthread
LIBS = -lz

OBJS := main input reads
BINS := $(addprefix $(BINDIR)/, $(OBJS))

GFALIBS_DIR := $(CURDIR)/gfalibs
GFALIBS_FILES := $(GFALIBS_DIR)/$(SOURCE)/* $(GFALIBS_DIR)/$(INCLUDE)/*

head: $(BINS) $(GFALIBS_FILES) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD)/$(TARGET) $(BINDIR)/* $(GFALIBS_DIR)/*.o $(LIBS)

$(BINDIR)%: $(SOURCE)/%.cpp $(INCLUDE)/%.h | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $(SOURCE)/$(notdir $@).cpp -o $@

$(GFALIBS_FILES): gfalibs
	@# Do nothing


.PHONY: gfalibs
gfalibs:
	$(MAKE) -j -C $(GFALIBS_DIR)
	
$(BUILD):
	-mkdir -p $@
	
$(BINDIR):
	-mkdir -p $@
	
clean:
	$(MAKE) -j -C $(GFALIBS_DIR) clean
	$(RM) -r build
