package big

import (
	"fmt"
	"sync"
)

// A. Schönhage and V. Strassen, "Schnelle Multiplikation großer Zahlen",
//      Computing 7 (1971), pp. 281–292.
//      https://link.springer.com/article/10.1007/BF02242355
//    English translation by Ryan Landay, 2023
//      https://github.com/rlanday/FastMultiplicationOfLargeIntegers
//
// Pierrick Gaudry, Alexander Kruppa,  Paul Zimmermann
//    "A GMP-based Implementation of Schönhage-Strassen’s Large
//       Integer Multiplication Algorithm"
//    ISSAC '07: Proceedings of the 2007 international symposium on
//        Symbolic and algebraic computation
//        https://dl.acm.org/doi/proceedings/10.1145/1277548
//  Available at https://inria.hal.science/inria-00126462v2/document
//
// Wikipedia, entry for "Schönhage–Strassen algorithm"
//   https://en.wikipedia.org/wiki/Sch%C3%B6nhage%E2%80%93Strassen_algorithm
//

// ssaThreshold is minimum number of words needed in the product
// for this algorithm to be used.
var ssaThreshold = 1000

// ssaMul performs the Schönhage-Strassen algorithm computing
//
//	x*y storing result in z
func ssaMul(stk *stack, z, x, y nat) nat {
	n := len(x) + len(y)
	//ssaDump(n)
	k := ssaMulParam(n)
	return ssaMulK(stk, k, z, x, y).norm()
}

// ssaMulParam given product length, returns "k" for the FFT algorithm
func ssaMulParam(n int) (k int) {
	// get best k for fft
	const k0 int = 6
	var lenTable = []int{
		// this will be redone later
		ssaThreshold, 2 * ssaThreshold, 4 * ssaThreshold,
		8 * ssaThreshold, 24 * ssaThreshold, 72 * ssaThreshold,
	}

	var i int
	for i = range len(lenTable) {
		if n < lenTable[i] {
			k = i + k0
			break
		}
	}
	if k == 0 {
		if i == 0 || n < 4*lenTable[i-1] {
			k = i + k0
		} else {
			k = i + k0 + 1
		}
	}

	return
}

func ssaDump(n int) {
	fmt.Println("n initial", n)

	k := ssaMulParam(n)

	n = 1 + ((n - 1) >> uint(k)) // ceil(n/2ᵏ)
	n <<= uint(k)                // ceil(n/2ᵏ) * 2ᵏ

	N := n * _W       // total bits of z, output
	M := N >> uint(k) // divide N to 2ᵏ terms and per part is M bits
	l := n >> uint(k) // per term has l words
	K := 1 << uint(k) // K=2ᵏ

	// M = N / K or  N = M * 2ᵏ
	// l = n / K
	// get prime for fft which is 2^Nprime+1.
	maxLK := max(K, _W)                     // ensure Nprime%_W = 0
	Nprime := (1 + (2*M+k+2)/maxLK) * maxLK // total bits of prime
	nprime := Nprime / _W

	// divide Nprime to 2^k terms.
	// 2^(Mp*K)   mod 2^Nprime+1 = -1
	// 2^(2*Mp*K) mod 2^Nprime+1 = -1
	Mp := Nprime >> uint(k)

	fmt.Println("n", n, "words")
	fmt.Println("k", k)
	fmt.Println("N", N, "bits", "N = M * 2ᵏ")
	fmt.Println("M", M, "bits")
	fmt.Println("l", l, "words")
	fmt.Println("K", K, "2ᵏ")
	fmt.Println("LX", maxLK)
	fmt.Println("NPrime", Nprime, "bits")
	fmt.Println("nprimp", nprime, "words")
	fmt.Println("Mp", Mp, "Nprime = Mp * 2ᵏ")
	fmt.Println("len(A)", K*(nprime+1), "words")
	fmt.Println("tmp", 2*(nprime+1), "words")
}

func ssaMulVector(ap, bp []nat, n int, k int) {
	//fmt.Println("MULVEC GOT n=", n, "k=", k)
	//
	K := 1 << uint(k) // K=2ᵏ

	N := n * _W       // total bits of z, output
	M := N >> uint(k) // divide N to 2ᵏ terms and per part is M bits
	l := n >> uint(k) // per term has l words
	// get prime for fft which is 2^Nprime+1.

	// maxLK is almost always K
	maxLK := max(K, _W) // ensure Nprime%_W = 0

	Nprime := (1 + (2*M+k+2)/maxLK) * maxLK // total bits of prime
	nprime := Nprime / _W
	Mp := Nprime >> uint(k)

	// allocation nat directly instead of "nat(nil).make" which
	// adds 4 extra words at the end.
	A := make(nat, K*(nprime+1))
	B := make(nat, K*(nprime+1))
	T := make(nat, 2*(nprime+1))

	Ap := make([]nat, K)
	Bp := make([]nat, K)

	// Extend x,y to N bits then decompose it to K=2ᵏ terms
	// Each Ap[i] is a slice into A, with nprime+1 words
	// Loop putting l words of x into Ap[i]
	//
	// Same for B, Bp, and y.
	//

	lenvec := (l << k) + 1
	for i := range K {
		Ap[i] = A[i*(nprime+1) : (i+1)*(nprime+1)]
		Bp[i] = B[i*(nprime+1) : (i+1)*(nprime+1)]
		start := i * l

		x := ap[i]
		y := bp[i]

		if start < len(x) {
			end := min(start+l, lenvec)
			copy(Ap[i], x[start:end])
		} // else Ap[i] is all zeros
		if start < len(y) {
			end := min(start+l, lenvec)
			copy(Bp[i], y[start:end])
		} // else Bp[i] is all zeros
	}

	for i := range K {
		ssaMulCore(Ap[i], A, B, Ap, Bp, T, k, l, Mp, nprime)
		//ap[i][n] = carry
	}
}

func fftTermMul(Ap []nat, Bp []nat, nprime int, k int, T nat) {

	if false && nprime >= 300 {
		ssaMulVector(Ap, Bp, nprime, k)
		return
	}
	//	fmt.Println("TERM MUL=", nprime, "k=", k, "AP=", len(Ap))
	// term multiplications (mod 2^Nprime+1)
	//  Ap[i] * Bp[i]
	//
	tp := T[:2*nprime]
	for i := range len(Ap) {
		var cc Word
		a, b := Ap[i], Bp[i]

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
		if subVV(a[:nprime], tp[:nprime], tp[nprime:]) != 0 && addVW(a[:nprime], a[:nprime], 1) != 0 {
			a[nprime] = 1
		} else {
			a[nprime] = 0
		}
	}
}

func ssaMulCore(z nat, A nat, B nat, Ap []nat, Bp []nat, T nat, k int, l int, Mp int, nprime int) nat {
	//n := len(z)

	// Get convolution order for FFT
	fftOrder := fftOrderK(k)

	fftDirect(Ap, len(Ap), fftOrder, k, 2*Mp, nprime, 1, T)
	fftDirect(Bp, len(Ap), fftOrder, k, 2*Mp, nprime, 1, T)

	fftTermMul(Ap, Bp, nprime, k, T)

	/*
		// term multiplications (mod 2^Nprime+1)
		//  Ap[i] * Bp[i]
		//
		tp := T[:2*nprime]
		for i := range K {
			var cc Word
			a, b := Ap[i], Bp[i]

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
			if subVV(a[:nprime], tp[:nprime], tp[nprime:]) != 0 && addVW(a[:nprime], a[:nprime], 1) != 0 {
				a[nprime] = 1
			} else {
				a[nprime] = 0
			}
		}
	*/
	fftInverse(Ap, len(Ap), 2*Mp, nprime, T)

	for i := range len(Ap) {
		fermatDiv2Exp(Bp[i], Ap[i], k, nprime)
	}

	// addition of terms in result p
	np1 := nprime + 1
	Km1 := len(Ap) - 1
	pla := l*Km1 + np1

	// Reuse A for summation as it's no longer in use
	p := A[:pla]
	clear(p)

	for i, sh := Km1, l*Km1; i >= 0; i, sh = i-1, sh-l {
		t := p[sh:]
		if addVV(t[:np1], t[:np1], Bp[i][:np1]) != 0 {
			addVW(t[np1:pla-sh], t[np1:pla-sh], 1)
		}
	}

	// Copy final result to output z
	//
	copy(z, p[:len(z)])
	return z
	//return z.norm()
	//return fermatNormalize(z, len(z)-1)
}

func ssaMulK(stk *stack, k int, z, x, y nat) nat {
	//
	K := 1 << uint(k) // K=2ᵏ

	// compute n such that:
	//
	//  1.  n > xlen+ylen (must be able to hold result)
	//  2.  n is a power of 2 (needed for FFT)
	//  3.  n is divisible by 2ᵏ (a Fermant Number, in Ring ℤ/Fnℤ )
	//
	n := len(x) + len(y)         // number of words in output
	n = 1 + ((n - 1) >> uint(k)) // ceil(n/2ᵏ)
	n <<= uint(k)                // ceil(n/2ᵏ) * 2ᵏ

	N := n * _W // total bits of z, output

	// Compute M and l, such that:
	//
	//    N = M * 2ᵏ bits
	//    n = l * 2ᵏ words
	//
	M := N >> uint(k) // divide N to 2ᵏ terms and per part is M bits
	l := n >> uint(k) // per term has l words

	// get prime for fft which is 2^Nprime+1.

	// maxLK is almost always K
	maxLK := max(K, _W) // ensure Nprime%_W = 0

	// Note: "Nprime" has nothing to do with prime numbers.
	//       Common notation of this value is N' ("n prime")
	Nprime := (1 + (2*M+k+2)/maxLK) * maxLK // total bits of prime
	nprime := Nprime / _W

	// divide Nprime to 2^k terms.
	// 2^(Mp*K)   mod 2^Nprime+1 = -1
	// 2^(2*Mp*K) mod 2^Nprime+1 = -1
	// Nprime = Mp * 2ᵏ
	Mp := Nprime >> uint(k)

	// allocation nat directly instead of "nat(nil).make" which
	// adds 4 extra words at the end.
	A := make(nat, K*(nprime+1))
	B := make(nat, K*(nprime+1))
	T := make(nat, 2*(nprime+1))

	Ap := make([]nat, K)
	Bp := make([]nat, K)

	// Extend x,y to N bits then decompose it to K=2ᵏ terms
	// Each Ap[i] is a slice into A, with nprime+1 words
	// Loop putting l words of x into Ap[i]
	//
	// Same for B, Bp, and y.
	//
	for i := range K {
		Ap[i] = A[i*(nprime+1) : (i+1)*(nprime+1)]
		Bp[i] = B[i*(nprime+1) : (i+1)*(nprime+1)]
		start := i * l
		if start < len(x) {
			end := min(start+l, len(x))
			copy(Ap[i], x[start:end])
		} // else Ap[i] is all zeros
		if start < len(y) {
			end := min(start+l, len(y))
			copy(Bp[i], y[start:end])
		} // else Bp[i] is all zeros
	}

	// when this function is called by nat.mul()
	// z.make(n) does nothing, as it has already been allocated
	z = z.make(n)

	if true {
		return ssaMulCore(z, A, B, Ap, Bp, T, k, l, Mp, nprime)
	} else {
		// Get convolution order for FFT
		fftOrder := fftOrderK(k)
		fftDirect(Ap, K, fftOrder, k, 2*Mp, nprime, 1, T)
		fftDirect(Bp, K, fftOrder, k, 2*Mp, nprime, 1, T)

		// term multiplications (mod 2^Nprime+1)
		//  Ap[i] * Bp[i]
		//
		tp := T[:2*nprime]
		for i := range K {
			var cc Word
			a, b := Ap[i], Bp[i]

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
			if subVV(a[:nprime], tp[:nprime], tp[nprime:]) != 0 && addVW(a[:nprime], a[:nprime], 1) != 0 {
				a[nprime] = 1
			} else {
				a[nprime] = 0
			}
		}

		fftInverse(Ap, K, 2*Mp, nprime, T)

		for i := range K {
			fermatDiv2Exp(Bp[i], Ap[i], k, nprime)
		}

		// addition of terms in result p
		np1 := nprime + 1
		Km1 := K - 1
		pla := l*Km1 + np1

		// Reuse A for summation as it's no longer in use
		p := A[:pla]
		clear(p)

		for i, sh := Km1, l*Km1; i >= 0; i, sh = i-1, sh-l {
			t := p[sh:]
			if addVV(t[:np1], t[:np1], Bp[i][:np1]) != 0 {
				addVW(t[np1:pla-sh], t[np1:pla-sh], 1)
			}
		}

		// Copy final result to output z
		//
		//
		copy(z, p[:len(z)])
		return z
	}
}

var (
	fftOrderAry [][]int
	fftOnce     sync.Once
)

func fftGenerateOrderK(k int) [][]int {
	// get order of terms of fft.
	//   fft[0]: 0 (not used)
	//   fft[1]: 0 1
	//   fft[2]: 0 2 1 3
	//   fft[3]: 0 4 2 6 1 5 3 7
	//
	fftOrder := make([][]int, k+1)
	for i := range fftOrder {
		fftOrder[i] = make([]int, 1<<uint(i))
	}
	for i := 1; i <= k; i++ {
		Kt := 1 << uint(i-1)
		for j := 0; j < Kt; j++ {
			fftOrder[i][j] = fftOrder[i-1][j] * 2
			fftOrder[i][j+Kt] = 1 + fftOrder[i][j]
		}
	}
	return fftOrder
}

func fftOrderK(k int) [][]int {
	const kmax = 20
	fftOnce.Do(func() {
		fftOrderAry = fftGenerateOrderK(kmax)
	})
	if k <= kmax {
		return fftOrderAry
	}

	// possible but unlikely
	return fftGenerateOrderK(k)
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
func fftDirect(Ap []nat, K int, ll [][]int, layer int, omega int, n int, inc int, tp nat) {
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
	fftDirect(Ap, K2, ll, layer-1, 2*omega, n, 2*inc, tp)
	fftDirect(Ap[inc:], K2, ll, layer-1, 2*omega, n, 2*inc, tp)

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
func fftInverse(Ap []nat, K int, omega int, n int, tp nat) {
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
	fftInverse(Ap, K2, 2*omega, n, tp)
	fftInverse(Ap[K2:], K2, 2*omega, n, tp)

	for j := range K2 {
		// 2^x = - 2^(2*Nprime-x) (mod 2^Nprime+1)
		fermatMul2Exp(tp, Ap[K2+j], 2*n*_W-j*omega, n)
		fermatSub(Ap[K2+j], Ap[j], tp, n)
		fermatAdd(Ap[j], Ap[j], tp, n)
	}
}

// Arithmetic in mod 2^64+1
//
// For example, 2^64+1 (2^32+1 on 32-bit systems).
//
// Instead of binary words, let use mod (10^1) +1 or mod 11.
//
// Representing all values under mod 11 (or base 11), requires two decimal digits.
//  0 to 9, and 10.
//
// To add, two numbers mod 11, first add normally in bsae 10.
//
// If the first digit is 0, there is nothing to do.  "02 + 03" = 05.
// If the first digit is 1, there may be an overflow, and we need to normalize the number.
//
// 10 + 03 = 13
//
// To compute 13 mod 11, we would compute 13-11 == 2.
//
// However base (i.e. 11) may (and will) exceed a machine word size, requiring
// doing vector-vector (i.e. subVV).
//
// Instead, "subtract 11" can be done by "subtract 1, then subtract 10".
// Subtraction of 1, is done normally, on ther right most digits
// Subtrcation of 10 is done, by subtracting 1 on the left more digit.
//
// 10+3 = 13
//  3-1 = 2
//  1-1 = 0
// 02 is the answer
//
// 10+10 = 20
// 0 -1 = 9, with borrow
// 2 - 1 - borrow = 0
//
// 09 is the answer.
//

// fermatNormalize normalizes Fermat Numbers, mod 2^n + 1,
// that occur in the process of arithmatic. It does not normalize an artbitrary input.
func fermatNormalize(r nat, n int) nat {
	if r[n] == 0 {
		// already normal
		return r
	}

	// subtract 1

	// case 1 no borrow
	if cc := subVW(r[:n], r[:n], 1); cc == 0 {
		r[n] = 0
		return r
	}

	// case 2, with borrow
	//  when r = 2^Nprime
	clear(r[:n])
	r[n] = 1
	return r
}

func fermatAdd(r, a, b nat, n int) {
	c := a[n] + b[n] + addVV(r[:n], a[:n], b[:n])
	if c > 1 {
		r[n] = 1
		subVW(r[:n+1], a[:n+1], c-1)
	} else {
		r[n] = c
	}
}

func fermatSub(r, a, b nat, n int) {
	c := a[n] - b[n] - subVV(r[:n], a[:n], b[:n])

	// c is an unsigned type, and we did subtraction
	// so either check highbit, or do a cast to signed
	// if int(c) < 0 {
	if (c & (_M ^ (_M >> 1))) != 0 {
		r[n] = 0
		addVW(r[:n+1], r[:n+1], -c)
	} else {
		r[n] = c
	}
}

// fermatDiv2Exp performs a/2ᵈ mod 2ᴺ+1
func fermatDiv2Exp(r nat, a nat, d int, n int) {
	fermatMul2Exp(r, a, 2*n*_W-d, n)
	fermatNormalize(r, n)
}

// fermatMul2Exp performs a*2ᵈ mod 2ᴺ+1, ensure 0 <= d <= 2*N.
func fermatMul2Exp(r nat, a nat, d int, n int) {
	var rd, cc Word
	sh := uint(d % _W)
	m := d / _W

	if m >= n {
		m -= n
		if sh == 0 {
			copy(r[:m+1], a[n-m:n+1])
			rd = r[m]
			copy(r[m:n], a[:n-m])
			// cc remains zero
		} else {
			lshVU(r[:m+1], a[n-m:n+1], sh)
			rd = r[m]
			cc = lshVU(r[m:n], a[:n-m], sh)
		}
		for i, x := range r[m:n] {
			r[m+i] = ^x
		}
		r[n] = 0
		cc++
		addVW(r[:n+1], r[:n+1], cc)
		rd++
		if rd == 0 {
			cc = 1
		} else {
			cc = rd
		}
		t := m
		if rd == 0 {
			t++
		}
		addVW(r[t:n+1], r[t:n+1], cc)
		return
	}

	// m < n

	// TODO SIMPLIFY
	if sh == 0 {
		copy(r[:m+1], a[n-m:n+1])
		for i, x := range r[:m+1] {
			r[i] = ^x
		}
		rd = ^r[m]
		copy(r[m:n], a[:n-m])
		cc = 0
	} else {
		lshVU(r[:m+1], a[n-m:n+1], sh)
		for i, x := range r[:m+1] {
			r[i] = ^x
		}
		rd = ^r[m]
		cc = lshVU(r[m:n], a[:n-m], sh)
	}

	if m != 0 {
		if cc == 0 {
			cc = addVW(r[:n], r[:n], 1)
		} else {
			cc--
		}
		cc = subVW(r[:m], r[:m], cc) + 1
	}
	r[n] = -subVW(r[m:n], r[m:n], cc)
	r[n] -= subVW(r[m:n], r[m:n], rd)
	if r[n]&(_M^(_M>>1)) != 0 { // when r[n]<0
		r[n] = addVW(r[:n], r[:n], 1)
	}
}
