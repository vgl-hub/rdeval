CXX = g++
INCLUDE_DIR = -I./include -I./gfalibs/include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = rdeval
TEST_TARGET = validate
BUILD = build/bin
SOURCE = src
INCLUDE = include
BINDIR := $(BUILD)/.o

LDFLAGS := -pthread
LIBS = -lz

OBJS := main input reads
BINS := $(addprefix $(BINDIR)/, $(OBJS))

#gfalibs
GFALIBS_DIR := $(CURDIR)/gfalibs

head: $(BINS) gfalibs | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD)/$(TARGET) $(BINDIR)/* $(GFALIBS_DIR)/*.o $(LIBS)
	
all: head validate

$(BINDIR)%: $(SOURCE)/%.cpp $(INCLUDE)/%.h | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $(SOURCE)/$(notdir $@).cpp -o $@

.PHONY: gfalibs
gfalibs:
	$(MAKE) -j -C $(GFALIBS_DIR)
	
validate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(TEST_TARGET) $(SOURCE)/$(TEST_TARGET).cpp $(LIBS)
	
$(BUILD):
	-mkdir -p $@
	
$(BINDIR):
	-mkdir -p $@
	
clean:
	$(MAKE) -j -C $(GFALIBS_DIR) clean
	$(RM) -r build
