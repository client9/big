build:
	go build .

test:
	go test .

fuzz:
	go test -fuzz=NatMulFFT

# runs the Net multiplication test
bench:
	go test -bench NatMul -benchmem

# runs the Nat multiplication tests without ASM optimizations
bench2:
	go test -tags math_big_pure_go -bench NatMul -benchmem

#
setup:
	go install golang.org/x/perf/cmd/benchstat@latest

bstat:
	go test -bench NatMul -benchmem -count=10 > old.txt
