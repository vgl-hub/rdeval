CXX ?= g++
INCLUDE_DIR += -I./include -I./gfalibs/include -I./htslib
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS) $(CFLAGS)

TARGET = rdeval
TEST_TARGET = validate
GENERATE_TARGET = generate-tests
BUILD = build/bin
SOURCE = src
INCLUDE = include
BINDIR := $(BUILD)/.o

LDFLAGS += -pthread
LIBS = -lz -lcrypto -L./htslib -lhts

OBJS := main input reads
BINS := $(addprefix $(BINDIR)/, $(OBJS))

#gfalibs
GFALIBS_DIR := $(CURDIR)/gfalibs
#htslib
HTSLIB_DIR := $(CURDIR)/htslib

head: $(BINS) gfalibs htslib | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD)/$(TARGET) $(BINDIR)/* $(GFALIBS_DIR)/*.o $(LIBS)
	
all: head validate regenerate

$(BINDIR)%: $(SOURCE)/%.cpp $(INCLUDE)/%.h | $(BINDIR)
	$(CXX) $(CXXFLAGS) -c $(SOURCE)/$(notdir $@).cpp -o $@

.PHONY: gfalibs
gfalibs:
	$(MAKE) -j -C $(GFALIBS_DIR) -I../htslib
	
.PHONY: htslib
htslib:
ifeq ($(wildcard $(HTSLIB_DIR)/configure),)
	cd $(HTSLIB_DIR) && autoreconf -i && ./configure
else
	@echo "File exists!"
endif
	$(MAKE) -j -C $(HTSLIB_DIR)
	$(MAKE) -j -C $(HTSLIB_DIR) install
	
validate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(TEST_TARGET) $(SOURCE)/$(TEST_TARGET).cpp

regenerate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(GENERATE_TARGET) $(SOURCE)/$(GENERATE_TARGET).cpp
	
$(BUILD):
	-mkdir -p $@
	
$(BINDIR):
	-mkdir -p $@

	
clean:
	$(MAKE) -j -C $(GFALIBS_DIR) clean
	$(MAKE) -j -C $(HTSLIB_DIR) clean
	$(RM) -r build
