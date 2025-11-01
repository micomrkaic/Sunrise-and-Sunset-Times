# Makefile for sun_times
# Default build: make (or make all)
# Install:       make install [PREFIX=/usr/local] [DESTDIR=...]
# Clean:         make clean

.POSIX:
.SUFFIXES:
.PHONY: all build release debug clean install uninstall

# ---- Project ----
SRC    := sun_times_v2.c
BIN    := sun_times

# ---- Toolchain (override on command line if needed) ----
override CC = gcc        # or clang
override CFLAGS += -std=c17

CC      ?= gcc
CSTD    ?= -std=c17
WARN    ?= -Wall -Wextra -Wpedantic -Wshadow -Wconversion
OPT     ?= -O2
DEFS    ?=
CFLAGS  ?= $(CSTD) $(WARN) $(OPT) $(DEFS)
LDFLAGS ?=
LDLIBS  ?= -lm

# ---- Install paths ----
PREFIX  ?= /usr/local
BINDIR  ?= $(PREFIX)/bin

# ---- Objects ----
OBJ := $(SRC:.c=.o)

# ---- Targets ----
all: build
build: release
release: $(BIN)

debug: CFLAGS := $(CSTD) $(WARN) -O0 -g $(DEFS)
debug: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $(OBJ) $(LDLIBS)

# Single-source pattern rule
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	@echo "Cleaning..."
	@rm -f $(OBJ) $(BIN)

install: $(BIN)
	@echo "Installing to $(DESTDIR)$(BINDIR)"
	@mkdir -p "$(DESTDIR)$(BINDIR)"
	@install -m 0755 $(BIN) "$(DESTDIR)$(BINDIR)/$(BIN)"

uninstall:
	@echo "Removing $(DESTDIR)$(BINDIR)/$(BIN)"
	@rm -f "$(DESTDIR)$(BINDIR)/$(BIN)"
