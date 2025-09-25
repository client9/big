package big

var fftThreshold = 2500

// fftMulParam returns parameters n and k for the FFT algorithm
func fftMulParam(xlen, ylen int) int {
	// get best k for fft
	const k0 int = 6
	var lenTable = []int{
		fftThreshold, 2 * fftThreshold, 4 * fftThreshold, 8 * fftThreshold, 24 * fftThreshold, 72 * fftThreshold,
	}

	k := 0
	n := xlen + ylen
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

	return k
}

// References to A. Schönhage and V. Strassen, "Schnelle Multiplikation großer Zahlen", Computing 7 (1971), pp. 281–292.
// and https://en.wikipedia.org/wiki/Sch%C3%B6nhage%E2%80%93Strassen_algorithm#Convolution_theorem
func fftMul(stk *stack, z, x, y nat) nat {

	k := fftMulParam(len(x), len(y))
	return fftMulK(stk, k, z, x, y)
}

// n >= len(x) + len(y)
func fftMulK(stk *stack, k int, z, x, y nat) nat {

	// compute n such that:
	//
	//  1.  n > xlen+ylen (must be able to hold result)
	//  2.  n is a power of 2 (needed for FFT)
	//  3.  n is divisible by 2ᵏ (needed for ℤ/Fnℤ Fermat)

	n := len(x) + len(y)
	n = 1 + ((n - 1) >> uint(k)) // ceil(n/2ᵏ)
	n <<= uint(k)                // ceil(n/2ᵏ) * 2ᵏ

	z = z.make(n)
	N := n * _W       // total bits of z
	M := N >> uint(k) // Divide N to 2ᵏ terms and per part is M bits.
	l := n >> uint(k) // Per term has l words.
	K := 1 << uint(k) // K=2ᵏ

	// get order of terms of fft. fft[1]: 0 1, fft[2]: 0 2 1 3, fft[3]: 0 4 2 6 1 5 3 7, ...
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

	// get prime for fft which is 2^Nprime+1.
	maxLK := max(K, _W)                     // ensure Nprime%_W = 0
	Nprime := (1 + (2*M+k+2)/maxLK) * maxLK // total bits of prime
	nprime := Nprime / _W
	Mp := Nprime >> uint(k) // divide Nprime to 2^k terms. 2^(Mp*K) mod 2^Nprime+1 = -1. 2^(2*Mp*K) mod 2^Nprime+1 = -1.

	A := nat(nil).make(K * (nprime + 1)) // storage for fft
	B := nat(nil).make(K * (nprime + 1))
	Ap := make([]nat, K)
	Bp := make([]nat, K)
	T := nat(nil).make(2*nprime + 2) // temporary storage
	// Extend x,y to N bits then decompose it to 2^k terms
	for i := range K {
		Ap[i] = A[i*(nprime+1) : (i+1)*(nprime+1)]
		Bp[i] = B[i*(nprime+1) : (i+1)*(nprime+1)]
		start := i * l
		if start < len(x) {
			end := min(start+l, len(x))
			copy(Ap[i], x[start:end]) // put l words of x into Ap[i]
		} // else Ap[i] is all zeros
		if start < len(y) {
			end := min(start+l, len(y))
			copy(Bp[i], y[start:end]) // put l words of y into Bp[i]
		} // else Bp[i] is all zeros
	}

	// direct fft
	directFFT(Ap, K, fftOrder, k, 2*Mp, nprime, 1, T)
	directFFT(Bp, K, fftOrder, k, 2*Mp, nprime, 1, T)

	// term multiplications (mod 2^Nprime+1)
	tp := T[:2*nprime]

	for i := range K {
		var cc Word
		a := Ap[i]
		b := Bp[i]
		tp.mul(stk, a[:nprime], b[:nprime])
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

	// inverse fft
	inverseFFT(Ap, K, 2*Mp, nprime, T)

	// division of terms after inverse fft
	for i := range K {
		fermatMul2Exp(Bp[i], Ap[i], 2*Nprime-k, nprime)
		fermatNormalize(Bp[i], nprime)
	}

	// addition of terms in result p
	pla := l*(K-1) + nprime + 1
	p := A[:pla]
	clear(p) // needed
	i := K - 1
	sh := l * i
	for ; i >= 0; i-- {
		t := p[sh:]
		if addVV(t[:nprime+1], t[:nprime+1], Bp[i][:nprime+1]) != 0 {
			addVW(t[nprime+1:pla-sh], t[nprime+1:pla-sh], 1)
		}
		sh -= l
	}
	copy(z[:n], p[:n])
	return z.norm()
}

// calculate FFT of Ap.
// K <- number of terms of current layer.
// ll <- order of fft.
// layer <- the current number of fft layers.
// omega <- t*2*Mp in the current layer.
// n<- length of nprime.
// inc <-  the interval between the terms in the current layer.
//
// let w[x,t] = 2^(x*t*2*Mp) mod 2^Nprime+1. (Nprime = Mp * K)
// let A              =           [a0 a1 a2 ... a(K-1)]
// let  W[K,t]        =                   | w[0,t]  w[1,t]  w[2,t]  w[3,t]    ...  w[K-1,t]             |
//
//	                                | w[0,t]  w[2,t]  w[4,t]  w[6,t]    ...  w[2K-2,t]            |
//	| w[0,t]  w[3,t]  w[6,t]  w[9,t]    ...  w[3K-3,t]            |
//	|             ...                                                     |
//	|             ...                                                     |
//	| w[0,t] w[K-1,t] w[2K-2,t] w[3K-3,t]... w[(K-1)(K-1),t]      |
//
// then FFT(A,K,t) = A * W[K,t]=   [y0   y1   y2   y3  ...  y(K-1)]
//
// let A0 = [a0 a2 a4 ... a(K-2)], A1 = [a1 a3 a5 ... a(K-1)]
// so yi = (A0 + w[i,t]*A1) * [ w[0,t] w[2i,t] w[4i,t] ... w[(K-2)i,t] ].
// Because w[2i,t] = w[i,2t], so yi = (A0 + w[i,t]*A1)*[ w[0,2t] w[i,2t] w[2i,2t] ... w[(K-2)i/2,2t] ].
// Because w[i,t] = -w[i+K/2t,t], so y(K/2t+i) = (A0 - w[i,t])*[ w[0,2t] w[i,2t] w[2i,2t] ... w[(K-2)i/2,2t] ].
//
//	So FFT(A,K,t)[i]      = FFT(A0,K/2,2t)[i] + w[i,t] * DFT(A1,K/2,2t)[i],
//	   FFT(A,K,t)[i+K/2t] = FFT(A0,K/2,2t)[i] - w[i,t] * DFT(A1,K/2,2t)[i].
//	   i=0,1,2,...,K/2-1. ( called butterfly operation )
//
// So we successfully divided the problem into two half-length sub-problems.
//
// directFFT calculates FFT(Ap,K,1) and arranges the items in fft order in Ap.
func directFFT(Ap []nat, K int, ll [][]int, layer int, omega int, n int, inc int, tp nat) {
	if K == 2 {
		copy(tp[:n+1], Ap[0][:n+1])
		// Butterfly operation
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
	directFFT(Ap, K2, ll, layer-1, 2*omega, n, inc*2, tp)
	directFFT(Ap[inc:], K2, ll, layer-1, 2*omega, n, inc*2, tp)

	for j := range K2 {
		// lk[2*j] is the lower half of the coefficient.
		// This enables DFT to be arranged in order of fft to facilitate the calculation of inverse fft
		fermatMul2Exp(tp, Ap[inc], lk[2*j]*omega, n)
		fermatSub(Ap[inc], Ap[0], tp, n)
		fermatAdd(Ap[0], Ap[0], tp, n)

		if j < K2-1 {
			Ap = Ap[2*inc:]
		}
	}
}

// inverseFFT does the same thing as FFT operation, except w[x,t] is -(2^(x*t*2*Mp)) mod 2^Nprime+1 instead.
// K <- number of terms of current layer.
// omega<- t*2*Mp in the current layer.
// n<- length of nprime.
func inverseFFT(Ap []nat, K int, omega int, n int, tp nat) {
	if K == 2 {
		copy(tp[:n+1], Ap[0][:n+1])
		// Butterfly operation
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
	inverseFFT(Ap, K2, 2*omega, n, tp)
	inverseFFT(Ap[K2:], K2, 2*omega, n, tp)

	for j := range K2 {
		fermatMul2Exp(tp, Ap[K2], 2*n*_W-j*omega, n) // 2^x = - 2^(2*Nprime-x) (mod 2^Nprime+1)
		fermatSub(Ap[K2], Ap[0], tp, n)
		fermatAdd(Ap[0], Ap[0], tp, n)
		Ap = Ap[1:]
	}
}

func fermatNormalize(r nat, n int) {
	if r[n] == 0 {
		return
	}
	cc := subVW(r[:n], r[:n], 1)
	if cc != 0 { // only when r = 2^Nprime
		clear(r[:n])
		r[n] = 1
	} else {
		r[n] = 0
	}
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
	if (c & (_M ^ (_M >> 1))) != 0 { // if c<0
		r[n] = 0
		addVW(r[:n+1], r[:n+1], -c)
	} else {
		r[n] = c
	}
}

// r=a*2^d mod 2^(n*_W)+1, ensure 0 <= d <= 2*n*_W.
func fermatMul2Exp(r nat, a nat, d int, n int) {
	sh := uint(d % _W)
	m := d / _W
	var rd, cc Word
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
