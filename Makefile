CXX = g++
INCLUDE_DIR = -I./include -I./gfalibs/include -I./htslib
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = rdeval
TEST_TARGET = validate
GENERATE_TARGET = generate-tests
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

#htslib
HTSLIB_DIR := $(CURDIR)/htslib

head: $(BINS) gfalibs htslib | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD)/$(TARGET) $(BINDIR)/* $(GFALIBS_DIR)/*.o -L/Users/gformenti/Documents/GitHub/rdeval/htslib/ -lhts $(LIBS)
	
all: head validate regenerate

$(BINDIR)%: $(SOURCE)/%.cpp $(INCLUDE)/%.h | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $(SOURCE)/$(notdir $@).cpp -o $@

.PHONY: gfalibs
gfalibs:
	$(MAKE) -j -C $(GFALIBS_DIR)
	
.PHONY: htslib
htslib:
	$(MAKE) -j -C $(HTSLIB_DIR)
	
validate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(TEST_TARGET) $(SOURCE)/$(TEST_TARGET).cpp $(LIBS)

regenerate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(GENERATE_TARGET) $(SOURCE)/$(GENERATE_TARGET).cpp $(LIBS)
	
$(BUILD):
	-mkdir -p $@
	
$(BINDIR):
	-mkdir -p $@

	
clean:
	$(MAKE) -j -C $(GFALIBS_DIR) clean
	$(MAKE) -j -C $(HTSLIB_DIR) clean
	$(RM) -r build
