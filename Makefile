CC      ?= cc
PYTHON  ?= python3
# -D_DEFAULT_SOURCE: glibc hides M_PI / clock_gettime etc. under strict -std=c99.
CFLAGS  ?= -O2 -std=c99 -D_DEFAULT_SOURCE -Wall -Wextra -Wpedantic -Wshadow -Wconversion
# Sanitizer runs: some sandboxes need ASLR disabled, e.g.
#   make sanitize SAN_RUN="setarch $$(uname -m) -R" ASAN_OPTIONS=detect_leaks=0
SAN_RUN ?=
SAN_FLAGS := -O1 -g -std=c99 -D_DEFAULT_SOURCE -Wall -Wextra -fsanitize=address,undefined -fno-sanitize-recover=all
INC     := -Iinclude
BUILD   := build

LIB_SRC := src/energy.c
VIT_SRC := src/vitality.c
CLI_SRC := src/main.c
TEST_SRC := tests/test_energy.c
VTEST_SRC := tests/test_vitality.c

.PHONY: all lib cli test sanitize bench example check coverage docs clean

all: $(BUILD)/energy_cli $(BUILD)/test_energy $(BUILD)/test_vitality $(BUILD)/stream_demo $(BUILD)/vitality_demo

lib: $(BUILD)/libenergy.a

$(BUILD):
	@mkdir -p $(BUILD)

$(BUILD)/energy.o: $(LIB_SRC) include/energy.h | $(BUILD)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILD)/vitality.o: $(VIT_SRC) include/vitality.h include/energy.h | $(BUILD)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILD)/libenergy.a: $(BUILD)/energy.o $(BUILD)/vitality.o
	ar rcs $@ $^

$(BUILD)/energy_cli: $(CLI_SRC) $(BUILD)/libenergy.a include/energy.h include/vitality.h
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

$(BUILD)/test_energy: $(TEST_SRC) $(BUILD)/libenergy.a include/energy.h
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

$(BUILD)/test_vitality: $(VTEST_SRC) $(BUILD)/libenergy.a include/energy.h include/vitality.h
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

$(BUILD)/vitality_demo: examples/vitality_demo.c $(BUILD)/libenergy.a include/energy.h include/vitality.h
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

$(BUILD)/bench_energy: tests/bench_energy.c $(BUILD)/libenergy.a include/energy.h include/vitality.h
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

$(BUILD)/stream_demo: examples/stream_demo.c $(BUILD)/libenergy.a include/energy.h
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

bench: $(BUILD)/bench_energy
	@$(BUILD)/bench_energy

example: $(BUILD)/stream_demo $(BUILD)/vitality_demo
	@$(BUILD)/stream_demo
	@$(BUILD)/vitality_demo

test: $(BUILD)/test_energy $(BUILD)/test_vitality
	@$(BUILD)/test_energy
	@$(BUILD)/test_vitality

sanitize: | $(BUILD)
	$(CC) $(SAN_FLAGS) $(INC) $(LIB_SRC) $(VIT_SRC) $(CLI_SRC) -lm -o $(BUILD)/energy_cli_san
	$(CC) $(SAN_FLAGS) $(INC) $(LIB_SRC) $(TEST_SRC) -lm -o $(BUILD)/test_energy_san
	$(CC) $(SAN_FLAGS) $(INC) $(LIB_SRC) $(VIT_SRC) $(VTEST_SRC) -lm -o $(BUILD)/test_vitality_san
	@$(SAN_RUN) $(BUILD)/test_energy_san
	@$(SAN_RUN) $(BUILD)/test_vitality_san
	@echo "sanitizers: OK"

CHROME ?= $(shell command -v google-chrome 2>/dev/null || command -v chromium 2>/dev/null || echo "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")

docs:
	@test -x "$(CHROME)" || (echo "需要 Chrome/Chromium：make docs CHROME=/path/to/chrome"; exit 1)
	@mkdir -p $(BUILD)
	@"$(CHROME)" --headless=new --no-pdf-header-footer --no-sandbox --disable-gpu \
		--user-data-dir=$(CURDIR)/$(BUILD)/chrome-docs \
		--print-to-pdf=docs/technical.pdf "file://$(CURDIR)/docs/technical.html" 2>/dev/null
	@"$(CHROME)" --headless=new --no-pdf-header-footer --no-sandbox --disable-gpu \
		--user-data-dir=$(CURDIR)/$(BUILD)/chrome-docs \
		--print-to-pdf=docs/product.pdf "file://$(CURDIR)/docs/product.html" 2>/dev/null
	@"$(CHROME)" --headless=new --no-pdf-header-footer --no-sandbox --disable-gpu \
		--user-data-dir=$(CURDIR)/$(BUILD)/chrome-docs \
		--print-to-pdf=docs/vitality_technical.pdf "file://$(CURDIR)/docs/vitality_technical.html" 2>/dev/null
	@"$(CHROME)" --headless=new --no-pdf-header-footer --no-sandbox --disable-gpu \
		--user-data-dir=$(CURDIR)/$(BUILD)/chrome-docs \
		--print-to-pdf=docs/vitality_product.pdf "file://$(CURDIR)/docs/vitality_product.html" 2>/dev/null
	@echo "docs: docs/technical.pdf, docs/product.pdf, docs/vitality_technical.pdf, docs/vitality_product.pdf"

coverage: | $(BUILD)
	$(CC) -O0 -g --coverage -std=c99 -D_DEFAULT_SOURCE -Wall -Wextra $(INC) -c src/energy.c -o $(BUILD)/energy_cov.o
	$(CC) -O0 -g --coverage -std=c99 -D_DEFAULT_SOURCE -Wall -Wextra $(INC) -c src/vitality.c -o $(BUILD)/vitality_cov.o
	$(CC) -O0 -g --coverage -std=c99 -D_DEFAULT_SOURCE -Wall -Wextra $(INC) tests/test_energy.c $(BUILD)/energy_cov.o -lm -o $(BUILD)/test_cov
	$(CC) -O0 -g --coverage -std=c99 -D_DEFAULT_SOURCE -Wall -Wextra $(INC) tests/test_vitality.c $(BUILD)/energy_cov.o $(BUILD)/vitality_cov.o -lm -o $(BUILD)/test_vcov
	@$(BUILD)/test_cov >/dev/null
	@$(BUILD)/test_vcov >/dev/null
	@gcov -o $(BUILD) $(BUILD)/energy_cov.gcno | grep -A1 "energy.c"
	@gcov -o $(BUILD) $(BUILD)/vitality_cov.gcno | grep -A1 "vitality.c"

check: $(BUILD)/test_energy $(BUILD)/test_vitality $(BUILD)/energy_cli
	@$(BUILD)/test_energy
	@$(BUILD)/test_vitality
	@$(PYTHON) scripts/check_parity.py 3
	@echo "check: OK"

clean:
	rm -rf $(BUILD)
