.PHONY: test
test:
	julia --project=. -e 'using ReTestItems; runtests()'

.PHONY: format
format:
	julia --project=utils/format/ utils/format/format.jl

.PHONY: docs
docs:
	julia --project=docs docs/make.jl

.PHONY: benchmark
benchmark:
	julia --project=benchmark benchmark/runbenchmarks.jl

