build:
	go build .
	go test .


test:
	go test .

lint:
	gofmt -w -s *.go
	go vet *.go

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


prof:
	go test -cpuprofile cpu.prof -memprofile mem.prof -test.run=xxx -bench=BigOne

# Coverage report
# Go default is crap
# cov0 is red (not covered)
# cov8 is the green (covered)
cover:
	rm -f coverage*
	go test -coverpkg=. -coverprofile=coverage.out ./...
	go tool cover -html=coverage.out -o coverage.html-tmp
	cat coverage.html-tmp | sed 's/background: black/background: whitesmoke/g' | sed 's/80, 80, 80/0,0,0/g' | sed 's/Menlo/ui-monospace/g' | sed 's/bold/normal/g' | sed 's/rgb(192, 0, 0)/rgb(255,0,0);font-weight:bold;/g' > coverage.html
