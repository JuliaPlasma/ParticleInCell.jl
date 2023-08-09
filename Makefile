.PHONY: error
error:
	@echo "Please choose one of the following targets:"
	@echo "- test"
	@echo "  - test-basics"
	@echo "  - test-unit"
	@echo "  - test-integration"
	@echo "- format"
	@echo "- docs"
	@echo "- benchmark"
	@exit 2

.PHONY: test
test: test-basics test-unit test-integration

.PHONY: test-basics
test-basics:
	julia --project=. -e 'using ReTestItems; runtests(tags=:basics)'

.PHONY: test-unit
test-unit:
	julia --project=. -e 'using ReTestItems; runtests(tags=:unit)'

.PHONY: test-integration
test-integration:
	julia --project=. -e 'using ReTestItems; runtests(tags=:integration)'

.PHONY: format
format:
	julia --project=utils/format/ utils/format/format.jl

.PHONY: docs
docs:
	julia --project=docs docs/make.jl

.PHONY: benchmark
benchmark:
	julia --project=benchmark benchmark/runbenchmarks.jl

