CC      ?= cc
# -D_DEFAULT_SOURCE: glibc 在 -std=c99 下需要它才能看到 M_PI。
CFLAGS  ?= -O2 -std=c99 -D_DEFAULT_SOURCE -Wall -Wextra -Wpedantic -Wshadow -Wconversion
INC     := -Iinclude
BUILD   := build

LIB_SRC  := src/energy.c src/vitality.c src/evobeing_zhenghuimin.c
LIB_OBJ  := $(BUILD)/energy.o $(BUILD)/vitality.o $(BUILD)/evobeing_zhenghuimin.o
HDRS     := include/energy.h include/vitality.h include/evobeing_zhenghuimin.h src/energy_internal.h

.PHONY: all lib test example docs clean

all: $(BUILD)/energy_cli $(BUILD)/test_energy $(BUILD)/test_vitality $(BUILD)/test_evobeing $(BUILD)/stream_demo $(BUILD)/vitality_demo $(BUILD)/evobeing_demo

lib: $(BUILD)/libenergy.a

$(BUILD):
	@mkdir -p $(BUILD)

$(BUILD)/%.o: src/%.c $(HDRS) | $(BUILD)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILD)/libenergy.a: $(LIB_OBJ)
	ar rcs $@ $^

$(BUILD)/energy_cli: src/main.c $(BUILD)/libenergy.a $(HDRS)
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

$(BUILD)/test_energy: tests/test_energy.c $(BUILD)/libenergy.a $(HDRS)
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

$(BUILD)/test_vitality: tests/test_vitality.c $(BUILD)/libenergy.a $(HDRS)
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

$(BUILD)/test_evobeing: tests/test_evobeing.c $(BUILD)/libenergy.a $(HDRS)
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

$(BUILD)/stream_demo: examples/stream_demo.c $(BUILD)/libenergy.a $(HDRS)
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

$(BUILD)/vitality_demo: examples/vitality_demo.c $(BUILD)/libenergy.a $(HDRS)
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

$(BUILD)/evobeing_demo: examples/evobeing_demo.c $(BUILD)/libenergy.a $(HDRS)
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

test: $(BUILD)/test_energy $(BUILD)/test_vitality $(BUILD)/test_evobeing
	@$(BUILD)/test_energy
	@$(BUILD)/test_vitality
	@$(BUILD)/test_evobeing

example: $(BUILD)/evobeing_demo $(BUILD)/stream_demo $(BUILD)/vitality_demo
	@$(BUILD)/evobeing_demo
	@$(BUILD)/stream_demo
	@$(BUILD)/vitality_demo

CHROME ?= $(shell command -v google-chrome 2>/dev/null || command -v chromium 2>/dev/null || echo "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")

docs:
	@test -x "$(CHROME)" || (echo "需要 Chrome/Chromium：make docs CHROME=/path/to/chrome"; exit 1)
	@mkdir -p $(BUILD)
	@for d in technical product vitality_technical vitality_product; do \
		"$(CHROME)" --headless=new --no-pdf-header-footer --no-sandbox --disable-gpu \
			--user-data-dir=$(CURDIR)/$(BUILD)/chrome-docs \
			--print-to-pdf=docs/$$d.pdf "file://$(CURDIR)/docs/$$d.html" 2>/dev/null; \
	done
	@echo "docs: PDF 已更新"

clean:
	rm -rf $(BUILD)
