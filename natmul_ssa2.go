package big

import (
	"fmt"
)

type ssa struct {
	k  int
	n  int // words, ouput size
	N  int // n in bits  WHY
	l  int
	M  int
	Mp int

	nprime int // words, adjusted size

	A  nat
	B  nat
	Ap []nat
	Bp []nat

	T nat
}

func ssaMul2(z, x, y nat) nat {
	f := ssa{}
	n := len(x) + len(y)
	k := ssaMulParam(n)
	f.init(n, k)
	//f.dump()
	f.decompose(f.A, f.Ap, x)
	f.decompose(f.B, f.Bp, y)

	// when this function is called by nat.mul()
	// z.make(n) does nothing, as it has already been allocated
	z = z.make(n)

	return f.core(z).norm()
}

func (f *ssa) dump() {
	fmt.Println("SSA MUL 2")
	//fmt.Println("_W", _W, "bits")
	fmt.Println("n", f.n, "words")
	fmt.Println("k", f.k)
	fmt.Println("N", f.N, "bits", "N = M * 2ᵏ")
	fmt.Println("M", f.M, "bits")
	fmt.Println("l", f.l, "words")
	fmt.Println("nprime", f.nprime, "words")
	fmt.Println("Mp", f.Mp, "Nprime = Mp * 2ᵏ")
}

func (f *ssa) termMul() {

	if true && f.nprime >= 300 {
		next := ssa{}

		n := f.nprime // len(f.Ap[0])*2
		//:n := len(f.Ap[0])*2
		k := ssaMulParam(n)
		next.init(n, k)

		fmt.Println("RECURSE with nprime=", f.nprime, "k=", f.k)
		fmt.Printf("New nprime=%d, k=%d, K=%d\n", next.nprime, next.k, 1<<next.k)

		// recall f.Ap and f.Bp have same length
		for i := range len(f.Ap) {

			clear(next.A)
			clear(next.B)
			f.Ap[i] = fermatNormalize(f.Ap[i], f.nprime)
			f.Bp[i] = fermatNormalize(f.Bp[i], f.nprime)

			next.decompose(next.A, next.Ap, f.Ap[i])
			next.decompose(next.B, next.Bp, f.Bp[i])
			next.core(f.Ap[i])
			f.Ap[i][next.nprime] = 0
			//carry := next.core(f.Ap[i])
			//f.Ap[i] = fermatNormalize(carry, len(carry)-1)
		}
		return
	}

	// term multiplications (mod 2^Nprime+1)
	//  Ap[i] * Bp[i]
	//
	nprime := f.nprime

	//fmt.Println("V2 nprime=", nprime, "AP=", len(f.Ap))
	tp := f.T[:2*nprime]

	// recall f.Ap and f.Bp have same length
	for i := range len(f.Ap) {
		var cc Word
		a, b := f.Ap[i], f.Bp[i]

		//
		// a*b mod 2ᴺ+1
		//
		tp.mul(nil, a[:nprime], b[:nprime])
		if a[nprime] != 0 {
			cc = addVV(tp[nprime:], tp[nprime:], b[:nprime])
		}
		if b[nprime] != 0 {
			cc += addVV(tp[nprime:], tp[nprime:], a[:nprime]) + a[nprime]
		}
		if cc != 0 {
			addVW(tp, tp, cc)
		}

		// using tp, subtract first half from second half, and put result in a
		if subVV(a[:nprime], tp[:nprime], tp[nprime:]) != 0 && addVW(a[:nprime], a[:nprime], 1) != 0 {
			a[nprime] = 1
		} else {
			a[nprime] = 0
		}
	}
}

func (f *ssa) init(n, k int) {

	f.k = k
	K := 1 << uint(k) // K=2ᵏ

	// compute n such that:
	//
	//  1.  n > xlen+ylen (must be able to hold result)
	//  2.  n is a power of 2 (needed for FFT)
	//  3.  n is divisible by 2ᵏ (a Fermant Number, in Ring ℤ/Fnℤ )
	//
	n = 1 + ((n - 1) >> uint(f.k)) // ceil(n/2ᵏ)
	f.n = n << uint(f.k)           // ceil(n/2ᵏ) * 2ᵏ

	f.N = f.n * _W // total bits of z, output

	// Compute M and l, such that:
	//
	//    N = M * 2ᵏ bits
	//    n = l * 2ᵏ words
	//
	f.M = f.N >> uint(f.k) // divide N to 2ᵏ terms and per part is M bits
	f.l = f.n >> uint(f.k) // per term has l words

	// get prime for fft which is 2^Nprime+1.

	// maxLK is almost always K
	maxLK := max(K, _W) // ensure Nprime%_W = 0

	// Note: "Nprime" has nothing to do with prime numbers.
	//       Common notation of this value is N' ("n prime")
	Nprime := (1 + (2*f.M+f.k+2)/maxLK) * maxLK // total bits of prime
	f.nprime = Nprime / _W

	// divide Nprime to 2^k terms.
	// 2^(Mp*K)   mod 2^Nprime+1 = -1
	// 2^(2*Mp*K) mod 2^Nprime+1 = -1
	// Nprime = Mp * 2ᵏ
	f.Mp = Nprime >> uint(f.k)

	// temporarly storage
	f.T = make(nat, 2*(f.nprime+1))

	// allocation nat directly instead of "nat(nil).make" which
	// adds 4 extra words at the end.
	f.A = make(nat, K*(f.nprime+1))
	f.B = make(nat, K*(f.nprime+1))

	f.Ap = make([]nat, K)
	f.Bp = make([]nat, K)

	chunkLen := f.nprime + 1
	for i, start := 0, 0; i < len(f.Ap); i, start = i+1, start+chunkLen {
		f.Ap[i] = f.A[start : start+chunkLen]
		f.Bp[i] = f.B[start : start+chunkLen]
	}
}

func (f *ssa) decompose(a nat, ap []nat, x nat) {
	// Extend x to N bits then decompose it to K=2ᵏ terms
	// Each Ap[i] is a slice into A, with nprime+1 words
	// Loop putting l words of x into Ap[i]
	//

	for i, start := 0, 0; i < len(ap) && start < len(x); i, start = i+1, start+f.l {
		end := min(start+f.l, len(x))
		copy(ap[i], x[start:end])
	}
}

func (f *ssa) core(z nat) nat {
	//n := len(z)

	K := len(f.Ap)

	// Get convolution order for FFT
	fftOrder := fftOrderK(f.k)

	fftDirect(f.Ap, K, fftOrder, f.k, 2*f.Mp, f.nprime, 1, f.T)
	fftDirect(f.Bp, K, fftOrder, f.k, 2*f.Mp, f.nprime, 1, f.T)

	f.termMul()

	fftInverse(f.Ap, K, 2*f.Mp, f.nprime, f.T)

	for i := range K {
		fermatDiv2Exp(f.Bp[i], f.Ap[i], f.k, f.nprime)
	}

	// addition of terms in result p
	np1 := f.nprime + 1
	Km1 := K - 1
	pla := f.l*Km1 + np1

	// Reuse A for summation as it's no longer in use
	p := f.A[:pla]
	clear(p)

	for i, sh := Km1, f.l*Km1; i >= 0; i, sh = i-1, sh-f.l {
		t := p[sh:]
		if addVV(t[:np1], t[:np1], f.Bp[i][:np1]) != 0 {
			addVW(t[np1:pla-sh], t[np1:pla-sh], 1)
		}
	}

	// Copy final result to output z
	//
	copy(z, p[:len(z)])
	return z
}
