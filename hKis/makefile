all:
	$(MAKE) -C "/home/chenzh/experiments/cec/SAT-sweeping/hKis/build"
kissat:
	$(MAKE) -C "/home/chenzh/experiments/cec/SAT-sweeping/hKis/build" kissat
tissat:
	$(MAKE) -C "/home/chenzh/experiments/cec/SAT-sweeping/hKis/build" tissat
clean:
	rm -f "/home/chenzh/experiments/cec/SAT-sweeping/hKis"/makefile
	-$(MAKE) -C "/home/chenzh/experiments/cec/SAT-sweeping/hKis/build" clean
	rm -rf "/home/chenzh/experiments/cec/SAT-sweeping/hKis/build"
coverage:
	$(MAKE) -C "/home/chenzh/experiments/cec/SAT-sweeping/hKis/build" coverage
indent:
	$(MAKE) -C "/home/chenzh/experiments/cec/SAT-sweeping/hKis/build" indent
test:
	$(MAKE) -C "/home/chenzh/experiments/cec/SAT-sweeping/hKis/build" test
.PHONY: all clean coverage indent kissat test tissat
