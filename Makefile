CXX ?= g++
INCLUDE_DIR += -I./include -I./gfalibs/include
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
LIBS = -lz -lcrypto -lhts

# Static linking settings
ifeq ($(UNAME_S),Linux)
    STATIC_LDFLAGS = -static -pthread
else
    STATIC_LDFLAGS = -pthread
endif

# Automatically populate static flags from pkg-config
PKG_CONFIG ?= pkg-config
export PKG_CONFIG_PATH
STATIC_CFLAGS = $(shell PKG_CONFIG_PATH="$(PKG_CONFIG_PATH)" $(PKG_CONFIG) --cflags --static htslib openssl zlib)
STATIC_LIBS_RAW = $(shell PKG_CONFIG_PATH="$(PKG_CONFIG_PATH)" $(PKG_CONFIG) --libs --static htslib openssl zlib)
STATIC_LIBS = $(shell echo $(STATIC_LIBS_RAW) | tr ' ' '\n' | sort -u | tr '\n' ' ')

OBJS := main input reads
BINS := $(addprefix $(BINDIR)/, $(OBJS))

#define
EVP := -DEVP

#gfalibs
GFALIBS_DIR := $(CURDIR)/gfalibs

head: $(BINS) gfalibs | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD)/$(TARGET) $(BINDIR)/* $(GFALIBS_DIR)/*.o $(LIBS)
	
# New static target
static: $(BINS) gfalibs | $(BUILD)
	$(CXX) $(CXXFLAGS) $(STATIC_LDFLAGS) -o $(BUILD)/$(TARGET) $(BINDIR)/* $(GFALIBS_DIR)/*.o $(STATIC_LIBS)
	
debug-pkg:
	@echo "PKG_CONFIG_PATH=$(PKG_CONFIG_PATH)"
	@$(PKG_CONFIG) --libs --static htslib
	
all: head validate regenerate

$(BINDIR)%: $(SOURCE)/%.cpp $(INCLUDE)/%.h | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(EVP) -c $(SOURCE)/$(notdir $@).cpp -o $@

.PHONY: gfalibs
gfalibs:
	$(MAKE) -j -C $(GFALIBS_DIR)

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
	$(RM) -r build
