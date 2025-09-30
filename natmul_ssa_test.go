package big

import (
	"fmt"
	"math/rand/v2"
	"testing"
)

func myK(z, x, y nat) nat {
	m := len(x)
	n := len(y)
	z = z.make(m + n)

	stk := getStack()
	defer stk.free()

	// Let x = x1:x0 where x0 is the same length as y.
	// Compute z = x0*y and then add in x1*y in sections
	// if needed.
	karatsuba(stk, z[:2*n], x[:n], y)

	if n < m {
		clear(z[2*n:])
		defer stk.restore(stk.save())
		t := stk.nat(2 * n)
		for i := n; i < m; i += n {
			t = t.mul(stk, x[i:min(i+n, len(x))], y)
			addTo(z[i:], t)
		}
	}

	return z.norm()
}

func TestSSADump(t *testing.T) {
	ssaDump(1000000)
}

func BenchmarkBigOne(b *testing.B) {
	const nwords = 100000
	x := rndNat(nwords)
	y := rndNat(nwords)

	// size it
	z := nat{}.mul(nil, x, y)

	b.Run("SSA1", func(b *testing.B) {
		for b.Loop() {
			ssaMul(nil, z, x, y)
		}
	})
	b.Run("SSA2", func(b *testing.B) {
		for b.Loop() {
			ssaMul2(z, x, y)
		}
	})
}

func TestBigOne(t *testing.T) {
	t.Skip()
	fmt.Println("IN BIGONE")
	const nwords = 500000
	x := rndNat(nwords)
	y := rndNat(nwords)

	//z0 := nat{}
	//z0 = z0.mul(nil, x, y)

	//z0 = myK(z0, x, y)
	//z0 = testMul(x,y)
	fmt.Println("----")
	z1 := nat{}
	z1 = ssaMul(nil, z1, x, y)

	fmt.Println("----")
	z2 := nat{}
	z2 = ssaMul2(z2, x, y)
	/*
		if z0.cmp(z1) != 0 {
			t.Errorf("Z0 != Z1")
		}
		if z0.cmp(z2) != 0 {
			t.Errorf("Z0 != Z2")
		}
	*/
	if len(z1) != len(z2) {
		t.Errorf("Got different len")
	}
	if z1.cmp(z2) != 0 {
		t.Errorf("Z1 != Z2")
	}
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
	n := 1 << 21
	for k := 8; k < 16; k++ {
		b.Run(fmt.Sprintf("%d;k=%d", n, k), func(b *testing.B) {

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
		z = ssaMul2(z, a, b)
		if z.cmp(testMul(a, b)) != 0 {
			t.Errorf("%s x %s failure", a.String(), b.String())
		}
	})
}
