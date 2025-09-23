package big

import (
	"math/rand/v2"
	"testing"
)

func randnat(rnd *rand.Rand, n int) nat {
	v := make([]Word, n)
	for i := range v {
		v[i] = Word(rnd.Uint())
	}
	return nat(v).norm()
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
		z = fftMul(z, a, b)
		if z.cmp(testMul(a, b)) != 0 {
			t.Errorf("%s x %s failure", a.String(), b.String())
		}
	})
}
