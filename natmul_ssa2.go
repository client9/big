package big

import (
	"fmt"
)

// addVVRagged adds two vectors together, assuming len(x) >= len(y)
// (add a short vector to a long vector)
func addVVRagged(z, x, y nat) Word {

	if len(z) != len(x) {
		panic("AddVVRAgged: output len != inputput len")
	}
	if len(x) < len(y) {
		panic("addVVRagged len(x) < len(y)")
	}
	n := len(y)
	carry := addVV(z[:n], x[:n], y)
	if len(x) > n {
		carry = addVW(z[n:], x[n:], carry)
	}
	return carry
}

// subVVRagged subtracts a shorter integer from a longer one.
func subVVRagged(z, x, y nat) Word {

	if len(z) != len(x) {
		panic("SubVVRAgged: output len != inputput len")
	}
	if len(x) < len(y) {
		panic("subVVRagged len(x) < len(y)")
	}
	n := len(y)
	carry := subVV(z[:n], x[:n], y)
	if len(x) > n {
		carry = subVW(z[n:], x[n:], carry)
	}
	return carry
}

func ssaMul3(z, x, y nat) nat {
	n := len(x) + len(y)
	k := ssaMulParam(n)
	n = scalen(n, k)

	ctx := &ssaContext{}
	ctx.recursiveThreshold = 64
	ssaContextInit(ctx, n, k)
	return ssaMulContext(ctx, z, x, y)
}

func ssaMulContext(ctx *ssaContext, z, x, y nat) nat {
	ssaDecompose(ctx, ctx.A, ctx.Ap, x)
	ssaDecompose(ctx, ctx.B, ctx.Bp, y)
	z = z.make(len(x) + len(y))
	ssaMain(ctx, z)
	return z.norm()
}

// scalen returns a scaled value of n, given k
//
// compute n such that:
//
//  1. n > xlen+ylen (must be able to hold result)
//  2. n is a power of 2 (needed for FFT)
//  3. n is divisible by 2ᵏ (a "Fermant Number", in Ring ℤ/Fnℤ )
//     rather n = v 2ᵏ  for some v.
//
// for FFT to work the input size needs to be expended
func scalen(n int, k int) int {
	n = 1 + ((n - 1) >> k) // ceil(n/2ᵏ)
	return n << k          // ceil(n/2ᵏ) * 2ᵏ
}

// lcmWordSize returns the least common multiple of the
// nat word size (_W) and 2^k
//
// _W << min(a,b)
// 5 + (_W >> 6)
// 32 = 2^5, 64 = 2^6
// W << (min(a, k) + a
func lcmWordSizeK(k int) int {
	b := k
	a := _W // word size in bits
	for a > 1 && b > 0 {
		a >>= 1
		b--
	}
	return a << k
}

// getNprimeMultipleOf
func nprimeMultipleOfNextK(nprime int) int {
	for {
		K := 1 << ssaMulParam(nprime)

		// nprime % K == 0
		if nprime&(K-1) == 0 {
			return nprime
		}
		nprime = (nprime + K - 1) & -K
	}
}

type ssaContext struct {
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

	recursiveThreshold int
}

func (f *ssaContext) dump() {
	fmt.Println("SSA MUL 3")
	fmt.Println("n", f.n, "words")
	fmt.Println("k", f.k)
	fmt.Println("N", f.N, "bits", "N = M * 2ᵏ")
	fmt.Println("M", f.M, "bits")
	fmt.Println("l", f.l, "words")
	fmt.Println("nprime", f.nprime, "words")
	fmt.Println("Mp", f.Mp, "Nprime = Mp * 2ᵏ")
}

func ssaContextInit(f *ssaContext, n, k int) {
	f.k = k
	K := 1 << uint(k) // K=2ᵏ
	f.n = n
	f.N = f.n * _W // total bits of z, output

	// Compute M and l, such that:
	//
	//    N = M * 2ᵏ bits
	//    n = l * 2ᵏ words
	//
	f.M = f.N >> uint(f.k) // divide N to 2ᵏ terms and per part is M bits
	f.l = f.n >> uint(f.k) // per term has l words

	// maxLK is almost always K
	maxLK := max(K, _W) // ensure Nprime%_W = 0
	max2 := lcmWordSizeK(k)
	if max2 != maxLK {
		fmt.Printf("MAXLK max=%d, lcm=%d\n", maxLK, max2)
		panic("MaxLK mismatch")
	}

	// Note: "Nprime" has nothing to do with prime numbers.
	//       Common notation of this value is N' ("n prime")
	Nprime := (1 + (2*f.M+f.k+2)/maxLK) * maxLK // total bits of prime

	Nprime2 := ((2*f.M + k + 2 + maxLK) / maxLK) * maxLK

	if Nprime != Nprime2 {
		fmt.Println("NPrime", Nprime, "Nprime2", Nprime2)
		panic("nprime calculation fail")
	}
	f.nprime = Nprime / _W

	if f.recursiveThreshold != 0 && f.nprime > f.recursiveThreshold {
		f.nprime = nprimeMultipleOfNextK(f.nprime)
		Nprime = f.nprime * _W
	}

	// divide Nprime to 2^k terms.
	// 2^(Mp*K)   mod 2^Nprime+1 = -1
	// 2^(2*Mp*K) mod 2^Nprime+1 = -1
	// Nprime = Mp * 2ᵏ
	f.Mp = Nprime >> uint(f.k)

	// temp storage
	f.T = make(nat, 2*(f.nprime+1))

	// allocation nat directly instead of "nat(nil).make" which
	// adds 4 extra words at the end.
	f.A = make(nat, K*(f.nprime+1))
	f.B = make(nat, K*(f.nprime+1))
	if (K * (f.nprime + 1)) != ((f.nprime + 1) << k) {
		panic("K calculation fail")
	}
	f.Ap = make([]nat, K)
	f.Bp = make([]nat, K)
}

// store in A/AP bits of X
func ssaDecompose(ctx *ssaContext, a nat, ap []nat, x nat) {
	// Extend x to N bits then decompose it to K=2ᵏ terms
	// Each Ap[i] is a slice into A, with nprime+1 words
	// Loop putting l words of x into Ap[i]
	//
	srcChunk := ctx.l
	Kl := len(ap) * srcChunk

	// TODO: original uses len(x) > Kl
	// however many times len(x) == Kl + 1
	//.  and x is alredy normalized this copy doesn't make sense
	if len(x) > Kl+1 {
		//fmt.Printf("Len(T)=%d, len(KL)=%d, len(x)=%d\n", len(ctx.T), Kl, len(x))
		z := make(nat, Kl+1)
		fermatNorm(z, x)
		/*
			if cc := subVVRagged(z[:Kl], x[:Kl], x[Kl:]); cc != 0 {
				addVW(z, z, cc)
			}
		*/
		//copy(x, z)
		//x = x[:Kl]
		x = z
	}

	clear(a)
	chunkLen := ctx.nprime + 1

	for i, start := 0, 0; i < len(ap); i, start = i+1, start+chunkLen {
		ap[i] = a[start : start+chunkLen]
	}

	tmp := ctx.T[:ctx.nprime+1]
	start := 0
	for i := range len(ap) {
		if start >= len(x) {
			break
		}
		end := min(start+srcChunk, len(x))
		// fermatMul2Exp input and output need to be same size
		// thus the copy of the slice into tmp
		// TODO:  seems like this could be fixed
		copy(tmp, x[start:end])
		clear(tmp[end-start:])
		fermatMul2Exp(ap[i], tmp, i*ctx.Mp, ctx.nprime)

		start += srcChunk
	}
	/*
		for i, start := 0, 0; i < len(ap) && start < len(x); i, start = i+1, start+srcChunk {
			end := min(start+srcChunk, len(x))
			// fermatMul2Exp input and output need to be same size
			// thus the copy of the slice into tmp
			//. TODO:  seems like this could be fixed
			clear(tmp)
			copy(tmp, x[start:end])
			fermatMul2Exp(ap[i], tmp, i*ctx.Mp, ctx.nprime)
		}
	*/
}

func ssaMain(f *ssaContext, z nat) Word {
	K := len(f.Ap)

	// Get convolution order for FFT
	fftOrder := fftOrderK(f.k)
	omega := 2 * f.Mp

	fftForward(f.Ap, K, fftOrder, f.k, omega, f.nprime, 1, f.T)
	fftForward(f.Bp, K, fftOrder, f.k, omega, f.nprime, 1, f.T)

	ssaTermMul(f)

	fftBackwards(f.Ap, K, omega, f.nprime, f.T)

	f.Bp[0] = f.T[f.nprime+1:]
	fermatDiv2Exp(f.Bp[0], f.Ap[0], f.k, f.nprime)
	for i := 1; i < K; i++ {
		f.Bp[i] = f.Ap[i-1]
		fermatDiv2Exp(f.Bp[i], f.Ap[i], f.k+(K-i)*f.Mp, f.nprime)
	}

	// addition of terms in result p
	np1 := f.nprime + 1
	Km1 := K - 1
	pla := f.l*Km1 + np1

	// Reuse A for summation as it's no longer in use
	p := f.B[:pla]
	clear(p)

	tmp := f.T[:np1]
	clear(tmp)

	var carry Word
	for i, sh, lo := Km1, f.l*Km1, f.l*Km1+f.nprime; i >= 0; i, sh, lo = i-1, sh-f.l, lo-f.l {
		t := p[sh:]
		j := (K - i) & (K - 1)

		carry += addVVRagged(t, t, f.Bp[j])

		// if Bp[j] > (i+1)*(2^(2*M))
		tmp[2*f.l] = Word(i + 1)
		if f.Bp[j].cmp(tmp) > 0 {
			carry -= subVW(t[:pla], t[:pla], 1)
			carry -= subVW(p[lo:], p[lo:], 1)
		}
	}

	signedCarry := int(carry)
	switch signedCarry {
	case 0:
		// normal case
	case -1:
		pl := len(z)
		if addVW(p[pla-pl:], p[pla-pl:], 1) != 0 {
			subVW(p[pla-pl-1:], p[pla-pl-1:], 1)
			subVW(p[pla-1:], p[pla-1:], 1)
		}
	case 1:
		// have not encountered this yet
		panic("SignedCarry == 1")
	default:
		// should never happen
		panic("Signed Carry is not zero")
	}

	// set z for normalizes, and return carry
	return fermatNorm(z, p)
}

func ssaTermMul(f *ssaContext) {
	nprime := f.nprime
	if f.recursiveThreshold != 0 && nprime >= f.recursiveThreshold {
		//fmt.Printf("Recur, nprime=%d, 2*len(Ap[0])=%d\n", nprime, 2*len(f.Ap[0]))
		n := nprime
		k := ssaMulParam(n)
		next := &ssaContext{}
		ssaContextInit(next, n, k)
		//fmt.Printf("Old L=%d, New L=%d\n", f.l, next.l)
		if f.nprime != next.n {
			fmt.Println(f.nprime, next.n)
			panic("N is changed")
		}
		if f.nprime%(1<<next.k) != 0 {
			msg := fmt.Sprintf("in term-mul: nprime=%d, next.k=%d", f.nprime, 1<<next.k)
			panic("in term-mul: n/ k != 0" + msg)
		}

		for i := range len(f.Ap) {
			fermatNormalize(f.Ap[i])
			fermatNormalize(f.Bp[i])
			ssaDecompose(next, next.A, next.Ap, f.Ap[i])
			ssaDecompose(next, next.B, next.Bp, f.Bp[i])
			// note!
			//  f.Ap[i][nprime] is the last element (most significant elemnt)
			//  f.Ap[i][:nprime] computes value excluding last element
			f.Ap[i][nprime] = ssaMain(next, f.Ap[i][:nprime])
		}
		return
	}

	// term multiplications (mod 2^Nprime+1)
	//  Ap[i] * Bp[i]
	//

	//fmt.Printf("TermMul %d x %d words\n", len(f.Ap), len(f.Ap[0]))
	tp := f.T[:2*nprime]
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
		if cc = subVV(a[:nprime], tp[:nprime], tp[nprime:]); cc != 0 {
			cc = addVW(a[:nprime], a[:nprime], 1)
		}
		a[nprime] = cc

		/*
			// using tp, subtract first half from second half, and put result in a
			if subVV(a[:nprime], tp[:nprime], tp[nprime:]) != 0 && addVW(a[:nprime], a[:nprime], 1) != 0 {
				a[nprime] = 1
			} else {
				a[nprime] = 0
			}
		*/
	}
}

// calculate FFT of Ap.
// K <- number of terms of current layer.
// ll <- order of fft.
// layer <- the current number of fft layers.
// omega <- t*2*Mp in the current layer.
// n<- length of nprime.
// inc <-  the interval between the terms in the current layer.
//
// let w[x,t]         = 2^(x*t*2*Mp) mod 2^Nprime+1. (Nprime = Mp * K)
// let A              = [a0 a1 a2 ... a(K-1)]
// let  W[K,t]        =
//
//	| w[0,t]  w[1,t]  w[2,t]  w[3,t]    ...  w[K-1,t]         |
//	| w[0,t]  w[2,t]  w[4,t]  w[6,t]    ...  w[2K-2,t]        |
//	| w[0,t]  w[3,t]  w[6,t]  w[9,t]    ...  w[3K-3,t]        |
//	|             ...                                         |
//	|             ...                                         |
//	| w[0,t] w[K-1,t] w[2K-2,t] w[3K-3,t]... w[(K-1)(K-1),t]  |
//
// then FFT(A,K,t) = A * W[K,t]=   [y0   y1   y2   y3  ...  y(K-1)]
//
// let A0 = [a0 a2 a4 ... a(K-2)]
//
//	A1 = [a1 a3 a5 ... a(K-1)]
//
// then yi = (A0 + w[i,t]*A1) * [ w[0,t] w[2i,t] w[4i,t] ... w[(K-2)i,t] ]
//
// w[2i,t] = w[i,2t]
//
//	yi = (A0 + w[i,t]*A1)*[ w[0,2t] w[i,2t] w[2i,2t] ... w[(K-2)i/2,2t] ]
//
// w[i,t] = -w[i+K/2t,t]
//
//	y(K/2t+i) = (A0 - w[i,t])*[ w[0,2t] w[i,2t] w[2i,2t] ... w[(K-2)i/2,2t] ]
//
// Thus FFT(A,K,t)[i]      = FFT(A0,K/2,2t)[i] + w[i,t] * DFT(A1,K/2,2t)[i]
//
//	FFT(A,K,t)[i+K/2t] = FFT(A0,K/2,2t)[i] - w[i,t] * DFT(A1,K/2,2t)[i]
//	i=0,1,2,...,K/2-1. ( called butterfly operation )
//
// Wwe successfully divided the problem into two half-length sub-problems.
//
// fftDirect calculates FFT(Ap,K,1) and arranges the items in fft order in Ap.
func fftForward(Ap []nat, K int, ll [][]int, layer int, omega int, n int, inc int, tp nat) {
	if K == 2 {
		copy(tp[:n+1], Ap[0][:n+1])
		addVV(Ap[0][:n+1], Ap[0][:n+1], Ap[inc][:n+1])
		cy := subVV(Ap[inc][:n+1], tp[:n+1], Ap[inc][:n+1])
		if Ap[0][n] > 1 {
			Ap[0][n] = 1 - subVW(Ap[0][:n], Ap[0][:n], Ap[0][n]-1)
		}
		if cy != 0 {
			Ap[inc][n] = addVW(Ap[inc][:n], Ap[inc][:n], ^Ap[inc][n]+1)
		}
		return
	}

	lk := ll[layer]
	K2 := K >> 1
	fftForward(Ap, K2, ll, layer-1, 2*omega, n, 2*inc, tp)
	fftForward(Ap[inc:], K2, ll, layer-1, 2*omega, n, 2*inc, tp)

	for j, start := 0, 0; j < K2; j, start = j+1, start+2*inc {
		// lk[2*j] is the lower half of the coefficient.
		// This enables DFT to be arranged in order of fft to
		// facilitate the calculation of inverse fft
		fermatMul2Exp(tp, Ap[start+inc], lk[2*j]*omega, n)
		fermatSub(Ap[start+inc], Ap[start], tp, n)
		fermatAdd(Ap[start], Ap[start], tp, n)
	}
}

// fftInverse does the same thing as FFT operation with
//
//	w[x,t] = -(2^(x*t*2*Mp)) mod 2^Nprime+1
//
// K <- number of terms of current layer.
// omega<- t*2*Mp in the current layer.
// n<- length of nprime.
func fftBackwards(Ap []nat, K int, omega int, n int, tp nat) {
	if K == 2 {
		copy(tp[:n+1], Ap[0][:n+1])
		addVV(Ap[0][:n+1], Ap[0][:n+1], Ap[1][:n+1])
		cy := subVV(Ap[1][:n+1], tp[:n+1], Ap[1][:n+1])
		if Ap[0][n] > 1 {
			Ap[0][n] = 1 - subVW(Ap[0][:n], Ap[0][:n], Ap[0][n]-1)
		}
		if cy != 0 {
			Ap[1][n] = addVW(Ap[1][:n], Ap[1][:n], ^Ap[1][n]+1)
		}
		return
	}

	K2 := K >> 1
	fftBackwards(Ap, K2, 2*omega, n, tp)
	fftBackwards(Ap[K2:], K2, 2*omega, n, tp)

	for j := range K2 {
		fermatMul2Exp(tp, Ap[K2+j], j*omega, n)
		fermatSub(Ap[K2+j], Ap[j], tp, n)
		fermatAdd(Ap[j], Ap[j], tp, n)
	}
}
