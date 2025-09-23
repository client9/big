package big

var fftThreshold = 2500

// shlVU is obsolete,but existing lshVU doesn't handle s==0 case?
func local_shlVU(z, x []Word, s uint) (c Word) {
	if s == 0 {
		copy(z, x)
		return 0
	}
	return lshVU(z, x, s)
}

// References to A. Schönhage and V. Strassen, "Schnelle Multiplikation großer Zahlen", Computing 7 (1971), pp. 281–292.
// and https://en.wikipedia.org/wiki/Sch%C3%B6nhage%E2%80%93Strassen_algorithm#Convolution_theorem
func fftMul(stk *stack, z, x, y nat) nat {
	xl := len(x)
	yl := len(y)
	var n, k, i int
	{
		// get best n,k for fft
		n = xl + yl
		const k0 int = 6
		var lenTable = []int{
			fftThreshold, 2 * fftThreshold, 4 * fftThreshold, 8 * fftThreshold, 24 * fftThreshold, 72 * fftThreshold,
		}
		for i = 0; i < len(lenTable); i++ {
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
		n = 1 + ((n - 1) >> uint(k)) // ceil(n/2^k)
		n <<= uint(k)
	}
	z = z.make(n)
	N := n * _W       // total bits of z
	M := N >> uint(k) // Divide N to 2^k terms and per part is M bits.
	l := n >> uint(k) // Per term has l words.
	K := 1 << uint(k) // K=2^k
	// get order of terms of fft. fft[1]: 0 1, fft[2]: 0 2 1 3, fft[3]: 0 4 2 6 1 5 3 7, ...
	fftOrder := make([][]int, k+1)
	for i = range fftOrder {
		fftOrder[i] = make([]int, 1<<uint(i))
	}
	for i = 1; i <= k; i++ {
		Kt := 1 << uint(i-1)
		for j := 0; j < Kt; j++ {
			fftOrder[i][j] = fftOrder[i-1][j] * 2
			fftOrder[i][j+Kt] = 1 + fftOrder[i][j]
		}
	}
	// get prime for fft which is 2^Nprime+1.
	maxLK := _W // ensure Nprime%_W = 0
	if maxLK < K {
		maxLK = K
	}
	Nprime := (1 + (2*M+k+2)/maxLK) * maxLK // total bits of prime
	nprime := Nprime / _W
	Mp := Nprime >> uint(k) // divide Nprime to 2^k terms. 2^(Mp*K) mod 2^Nprime+1 = -1. 2^(2*Mp*K) mod 2^Nprime+1 = -1.

	A := nat(nil).make(K * (nprime + 1)) // storage for fft
	B := nat(nil).make(K * (nprime + 1))
	Ap := make([]nat, K)
	Bp := make([]nat, K)
	T := nat(nil).make(2*nprime + 2) // temporary storage
	// Extend x,y to N bits then decompose it to 2^k terms
	for i = 0; i < K; i++ {
		Ap[i] = A[i*(nprime+1) : (i+1)*(nprime+1)]
		Bp[i] = B[i*(nprime+1) : (i+1)*(nprime+1)]
		start := i * l
		if start < len(x) {
			end := start + l
			if end > len(x) {
				end = len(x)
			}
			copy(Ap[i], x[start:end]) // put l words of x into Ap[i]
		} // else Ap[i] is all zeros
		if start < len(y) {
			end := start + l
			if end > len(y) {
				end = len(y)
			}
			copy(Bp[i], y[start:end]) // put l words of y into Bp[i]
		} // else Bp[i] is all zeros
	}
	// direct fft
	directFFT(Ap, K, fftOrder, k, 2*Mp, nprime, 1, T)
	directFFT(Bp, K, fftOrder, k, 2*Mp, nprime, 1, T)

	// term multiplications (mod 2^Nprime+1)
	tp := T[:2*nprime]

	t := nat(nil)

	var cc Word
	for i := 0; i < K; i++ {
		a := Ap[i]
		b := Bp[i]
		t = t.mul(stk, a[:nprime], b[:nprime])
		clear(tp)
		copy(tp, t)
		if a[nprime] != 0 {
			cc = addVV(tp[nprime:2*nprime], tp[nprime:2*nprime], b[:nprime])
		} else {
			cc = 0
		}
		if b[nprime] != 0 {
			cc += addVV(tp[nprime:2*nprime], tp[nprime:2*nprime], a[:nprime]) + a[nprime]
		}
		if cc != 0 {
			addVW(tp[:2*nprime], tp[:2*nprime], cc)
		}
		if subVV(a[:nprime], tp[:nprime], tp[nprime:2*nprime]) != 0 && addVW(a[:nprime], a[:nprime], 1) != 0 {
			a[nprime] = 1
		} else {
			a[nprime] = 0
		}
	}

	// inverse fft
	inverseFFT(Ap, K, 2*Mp, nprime, T)

	// division of terms after inverse fft
	for i = 0; i < K; i++ {
		clear(Bp[i])
		mul2ExpMod(Bp[i], Ap[i], 2*Nprime-k, nprime)
		if Bp[i][nprime] != 0 {
			cc = subVW(Bp[i][:nprime], Bp[i][:nprime], 1)
			if cc != 0 { // only when Bp[i] = 2^Nprime
				clear(Bp[i][:nprime])
				Bp[i][nprime] = 1
			} else {
				Bp[i][nprime] = 0
			}
		}
	}

	// addition of terms in result p
	pla := l*(K-1) + nprime + 1
	p := A[:pla]
	clear(p)
	i = K - 1
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
	} else {
		lk := ll[layer]
		K2 := K >> 1
		directFFT(Ap, K2, ll, layer-1, 2*omega, n, inc*2, tp)
		directFFT(Ap[inc:], K2, ll, layer-1, 2*omega, n, inc*2, tp)
		for j := 0; j < K2; j++ {
			// lk[2*j] is the lower half of the coefficient.
			// This enables DFT to be arranged in order of fft to facilitate the calculation of inverse fft
			mul2ExpMod(tp, Ap[inc], lk[2*j]*omega, n)
			// Butterfly operation
			c := Ap[0][n] - tp[n] - subVV(Ap[inc][:n], Ap[0][:n], tp[:n])
			if (c & (_M ^ (_M >> 1))) != 0 { // if c<0
				Ap[inc][n] = 0
				addVW(Ap[inc][:n+1], Ap[inc][:n+1], -c)
			} else {
				Ap[inc][n] = c
			}

			c = Ap[0][n] + tp[n] + addVV(Ap[0][:n], Ap[0][:n], tp[:n])
			if c > 1 {
				Ap[0][n] = 1
				subVW(Ap[0][:n+1], Ap[0][:n+1], c-1)
			} else {
				Ap[0][n] = c
			}
			if j < K2-1 {
				Ap = Ap[2*inc:]
			}
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
	} else {
		K2 := K >> 1
		inverseFFT(Ap, K2, 2*omega, n, tp)
		inverseFFT(Ap[K2:], K2, 2*omega, n, tp)
		for j := 0; j < K2; j++ {
			mul2ExpMod(tp, Ap[K2], 2*n*_W-j*omega, n) // 2^x = - 2^(2*Nprime-x) (mod 2^Nprime+1)
			// Butterfly operation
			c := Ap[0][n] - tp[n] - subVV(Ap[K2][:n], Ap[0][:n], tp[:n])
			if (c & (_M ^ (_M >> 1))) != 0 { // if c<0
				Ap[K2][n] = 0
				addVW(Ap[K2][:n+1], Ap[K2][:n+1], -c)
			} else {
				Ap[K2][n] = c
			}

			c = Ap[0][n] + tp[n] + addVV(Ap[0][:n], Ap[0][:n], tp[:n])
			if c > 1 {
				Ap[0][n] = 1
				subVW(Ap[0][:n+1], Ap[0][:n+1], c-1)
			} else {
				Ap[0][n] = c
			}
			Ap = Ap[1:]
		}
	}
}

// r=a*2^d mod 2^(n*_W)+1, ensure 0 <= d <= 2*n*_W.
func mul2ExpMod(r nat, a nat, d int, n int) {
	sh := uint(d % _W)
	m := d / _W
	var rd, cc Word
	if m >= n {
		m -= n
		local_shlVU(r[:m+1], a[n-m:n+1], sh)
		rd = r[m]
		cc = local_shlVU(r[m:n], a[:n-m], sh)
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
	} else {
		local_shlVU(r[:m+1], a[n-m:n+1], sh)
		for i, x := range r[:m+1] {
			r[i] = ^x
		}
		rd = ^r[m]
		cc = local_shlVU(r[m:n], a[:n-m], sh)
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
}
