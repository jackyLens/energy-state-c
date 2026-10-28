CC      ?= cc
PYTHON  ?= python3
CFLAGS  ?= -O2 -std=c99 -Wall -Wextra -Wpedantic -Wshadow -Wconversion
INC     := -Iinclude
BUILD   := build

LIB_SRC := src/energy.c
CLI_SRC := src/main.c
TEST_SRC := tests/test_energy.c

.PHONY: all lib cli test sanitize bench example check coverage clean

all: $(BUILD)/energy_cli $(BUILD)/test_energy $(BUILD)/stream_demo

lib: $(BUILD)/libenergy.a

$(BUILD):
	@mkdir -p $(BUILD)

$(BUILD)/energy.o: $(LIB_SRC) include/energy.h | $(BUILD)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILD)/libenergy.a: $(BUILD)/energy.o
	ar rcs $@ $<

$(BUILD)/energy_cli: $(CLI_SRC) $(BUILD)/libenergy.a include/energy.h
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

$(BUILD)/test_energy: $(TEST_SRC) $(BUILD)/libenergy.a include/energy.h
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

$(BUILD)/bench_energy: tests/bench_energy.c $(BUILD)/libenergy.a include/energy.h
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

$(BUILD)/stream_demo: examples/stream_demo.c $(BUILD)/libenergy.a include/energy.h
	$(CC) $(CFLAGS) $(INC) $< $(BUILD)/libenergy.a -lm -o $@

bench: $(BUILD)/bench_energy
	@$(BUILD)/bench_energy

example: $(BUILD)/stream_demo
	@$(BUILD)/stream_demo

test: $(BUILD)/test_energy
	@$(BUILD)/test_energy

sanitize: | $(BUILD)
	$(CC) -O1 -g -std=c99 -Wall -Wextra $(INC) -fsanitize=address,undefined -fno-sanitize-recover=all \
		$(LIB_SRC) $(CLI_SRC) -lm -o $(BUILD)/energy_cli_san
	$(CC) -O1 -g -std=c99 -Wall -Wextra $(INC) -fsanitize=address,undefined -fno-sanitize-recover=all \
		$(LIB_SRC) $(TEST_SRC) -lm -o $(BUILD)/test_energy_san
	@$(BUILD)/test_energy_san
	@echo "sanitizers: OK"

coverage: | $(BUILD)
	$(CC) -O0 -g --coverage -std=c99 -Wall -Wextra $(INC) -c src/energy.c -o $(BUILD)/energy_cov.o
	$(CC) -O0 -g --coverage -std=c99 -Wall -Wextra $(INC) tests/test_energy.c $(BUILD)/energy_cov.o -lm -o $(BUILD)/test_cov
	@$(BUILD)/test_cov >/dev/null
	@gcov -o $(BUILD) $(BUILD)/energy_cov.gcno | grep -A1 "energy.c"

check: $(BUILD)/test_energy $(BUILD)/energy_cli
	@$(BUILD)/test_energy
	@$(PYTHON) scripts/check_parity.py 3
	@echo "check: OK"

clean:
	rm -rf $(BUILD)
