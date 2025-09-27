package big

import (
	"fmt"
	"math/rand/v2"
	"testing"
)

func TestSSADump(t *testing.T) {
	ssaDump(10000)
}

func randnat(rnd *rand.Rand, n int) nat {
	v := make([]Word, n)
	for i := range v {
		v[i] = Word(rnd.Uint())
	}
	return nat(v).norm()
}

func benchmarkFFTk(b *testing.B, nwords int, k int) {
	x := rndNat(nwords)
	y := rndNat(nwords)
	var z nat
	b.ReportAllocs()
	for b.Loop() {
		ssaMulK(nil, k, z, x, y)
	}
}

// 1000 - 7
// 2000 - 7
// 2500 - 8
// 3000 - 8
// 4000 - 8
// 5000 - 9
// 6000 - 8
// 8000 - 8
// 9000 - 9
// 10000 - 10
// 16000 - 9
// 32000 - 10
// 64000 - 11
// 128000 - 12
// 256000 - 13
// 512000 - 13
// 1000000 - 14
// 2000000 - 14
// 4000000 - 15
// 8000000 - 15
func BenchmarkFFTK(b *testing.B) {
	n := 8500000
	for k := 13; k < 20; k++ {
		b.Run(fmt.Sprintf("%d", k), func(b *testing.B) {
			benchmarkFFTk(b, n, k)
		})
	}
}

func FuzzNatMulFFT(f *testing.F) {
	pcg := rand.NewPCG(1, 2)
	r := rand.New(pcg)

	// set initial "corpus"
	f.Add(uint64(123456789), uint64(784727674))

	f.Fuzz(func(t *testing.T, s1, s2 uint64) {
		pcg.Seed(s1, s2)
		a := randnat(r, r.IntN(100))
		b := randnat(r, r.IntN(100))
		z := nat{}
		z = ssaMul(nil, z, a, b)
		if z.cmp(testMul(a, b)) != 0 {
			t.Errorf("%s x %s failure", a.String(), b.String())
		}
	})
}
