# Makefile for VcfAnalyte - Genetic Variant Analysis Tool:
#   make
#
# Usage          - Compile the program (optimized)
#   make debug    - Compile with debug symbols
#   make clean    - Remove compiled files
#   make test     - Run with sample data
#   make install  - Install to /usr/local/bin (requires sudo)

# Compiler settings
CC = gcc
CFLAGS = -std=c99 -Wall -Wextra -O2
LDFLAGS = -lm

# Target executable
TARGET = vcfanalyte

# Source files
SRCS = vcfanalyte.c
OBJS = $(SRCS:.c=.o)

# Test files
SAMPLE_VCF = sample.vcf
TEST_OUT = test_output.tsv
TEST_JSON = test_output.json

.PHONY: all clean test install help

# Default target
all: $(TARGET)

# Compile the main program
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

# Debug build
debug: CFLAGS += -g -DDEBUG
debug: $(TARGET)

# Clean build artifacts
clean:
	rm -f $(TARGET) $(OBJS) $(TEST_OUT) $(TEST_JSON)

# Run with sample data
test: $(TARGET)
	@echo "Running test with sample VCF..."
	./$(TARGET) $(SAMPLE_VCF) -s
	@echo ""
	@echo "Test with JSON output:"
	./$(TARGET) $(SAMPLE_VCF) --json -s > $(TEST_JSON)
	@echo "Output saved to $(TEST_JSON)"
	@echo ""
	@echo "Test with filtering:"
	./$(TARGET) $(SAMPLE_VCF) -t SNP -q 30 --stats
	@echo ""
	@echo "Test with custom chromosome:"
	./$(TARGET) $(SAMPLE_VCF) -c 1 -s
	@echo ""
	@echo "All tests completed!"

# Install to system
install: $(TARGET)
	cp $(TARGET) /usr/local/bin/
	chmod +x /usr/local/bin/$(TARGET)

# Show help
help:
	@echo "VcfAnalyte - Genetic Variant Analysis Tool"
	@echo ""
	@echo "Targets:"
	@echo "  make          - Build the program"
	@echo "  make debug    - Build with debug symbols"
	@echo "  make clean    - Remove build artifacts"
	@echo "  make test     - Run test with sample data"
	@echo "  make install  - Install to /usr/local/bin"
	@echo "  make help     - Show this help message"
