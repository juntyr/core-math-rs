/* Correctly-rounded hyperbolic arc-cosine for binary16 value.

Copyright (c) 2025 Paul Zimmermann

This file is ported from the CORE-MATH project
(https://core-math.gitlabpages.inria.fr/).

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

use hexf::hexf32 as h;

use super::utils::snanf16;

/* The following polynomials were generated using Sollya (cf acosh.sollya).
For 0 <= i < 16, P[i] is a degree-5 polynomial approximating
acosh(x)/sqrt(x-1) for 0 <= i < 5, and acosh(x) for 5 <= i < 16.
They were afterwards optimized to decrease the number of exceptions. */

// ([+-]?0x[0-9a-f]+\.[0-9a-f]+p[+-][0-9]+) -> h!("$1")
const P: [[f32; 6]; 16] = [
    [
        h!("0x1.91940ap+0"),
        h!("-0x1.a1ee52p-3"),
        h!("0x1.0dc64ap-4"),
        h!("-0x1.453108p-6"),
        h!("0x1.040a46p-8"),
        h!("-0x1.803198p-12"),
    ], /* [2^0,2^1) */
    [
        h!("0x1.8eb5fap+0"),
        h!("-0x1.6c52a2p-3"),
        h!("0x1.4dddcp-5"),
        h!("-0x1.e9e95cp-8"),
        h!("0x1.b697b2p-11"),
        h!("-0x1.5b2df4p-15"),
    ], /* [2^1,2^2) */
    [
        h!("0x1.854f5ap+0"),
        h!("-0x1.11e128p-3"),
        h!("0x1.362606p-6"),
        h!("-0x1.fb23bp-10"),
        h!("0x1.e2d03p-14"),
        h!("-0x1.8d21eap-19"),
    ], /* [2^2,2^3) */
    [
        h!("0x1.717322p+0"),
        h!("-0x1.608d74p-4"),
        h!("0x1.c3ceb2p-8"),
        h!("-0x1.8859d4p-12"),
        h!("0x1.828e08p-17"),
        h!("-0x1.44eab2p-23"),
    ], /* [2^3,2^4) */
    [
        h!("0x1.52b1b4p+0"),
        h!("-0x1.8e60acp-5"),
        h!("0x1.12e502p-9"),
        h!("-0x1.ee8866p-15"),
        h!("0x1.f10d92p-21"),
        h!("-0x1.a71f3ep-28"),
    ], /* [2^4,2^5) */
    [
        h!("0x1.1b4176p+1"),
        h!("0x1.c58fdp-4"),
        h!("-0x1.3d491cp-9"),
        h!("0x1.244e58p-15"),
        h!("-0x1.2b59e6p-22"),
        h!("0x1.029682p-30"),
    ], /* [2^5,2^6) */
    [
        h!("0x1.744536p+1"),
        h!("0x1.c4d402p-5"),
        h!("-0x1.3c787ap-11"),
        h!("0x1.235b58p-18"),
        h!("-0x1.2a375cp-26"),
        h!("0x1.017ceap-35"),
    ], /* [2^6,2^7) */
    // 4 exceptions
    [
        h!("0x1.cd16eap+1"),
        h!("0x1.c490cp-6"),
        h!("-0x1.3c27dep-13"),
        h!("0x1.22f70cp-21"),
        h!("-0x1.29b90ep-30"),
        h!("0x1.00fd26p-40"),
    ], /* [2^7,2^8) */
    // 1 exception
    [
        h!("0x1.12ec74p+2"),
        h!("0x1.c475e4p-7"),
        h!("-0x1.3c05dap-15"),
        h!("0x1.22cb2cp-24"),
        h!("-0x1.29804ap-34"),
        h!("0x1.00c29ap-45"),
    ], /* [2^8,2^9) */
    // 6 exceptions
    [
        h!("0x1.3f4c7ep+2"),
        h!("0x1.c45ddap-8"),
        h!("-0x1.3be5p-17"),
        h!("0x1.229e72p-27"),
        h!("-0x1.2943fp-38"),
        h!("0x1.008224p-50"),
    ], /* [2^9,2^10) */
    // // 1 exception
    [
        h!("0x1.6baa5p+2"),
        h!("0x1.c45576p-9"),
        h!("-0x1.3bd988p-19"),
        h!("0x1.228f6ep-30"),
        h!("-0x1.29308ep-42"),
        h!("0x1.006ddp-55"),
    ], /* [2^10,2^11) */
    // no exception
    [
        h!("0x1.980882p+2"),
        h!("0x1.c449a8p-10"),
        h!("-0x1.3bc944p-21"),
        h!("0x1.2278ep-33"),
        h!("-0x1.2911bap-46"),
        h!("0x1.004d3p-60"),
    ], /* [2^11,2^12) */
    // no exception
    [
        h!("0x1.c46664p+2"),
        h!("0x1.c44006p-11"),
        h!("-0x1.3bbbdp-23"),
        h!("0x1.22665ep-36"),
        h!("-0x1.28f8a4p-50"),
        h!("0x1.00326p-65"),
    ], /* [2^12,2^13) */
    // 1 exception
    [
        h!("0x1.f0c4d6p+2"),
        h!("0x1.c43246p-12"),
        h!("-0x1.3ba86ap-25"),
        h!("0x1.224b62p-39"),
        h!("-0x1.28d39ep-54"),
        h!("0x1.000a46p-70"),
    ], /* [2^13,2^14) */
    // 1 exception
    [
        h!("0x1.0e912cp+3"),
        h!("0x1.c42b4ep-13"),
        h!("-0x1.3b9eb4p-27"),
        h!("0x1.223e1p-42"),
        h!("-0x1.28c1ap-58"),
        h!("0x1.ffee3cp-76"),
    ], /* [2^14,2^15) */
    // no exception
    [
        h!("0x1.24bfccp+3"),
        h!("0x1.c4263p-14"),
        h!("-0x1.3b97b8p-29"),
        h!("0x1.22347ap-45"),
        h!("-0x1.28b4c4p-62"),
        h!("0x1.ffd2ecp-81"),
    ], /* [2^15,2^16) */ // no exception
];

/// Correctly-rounded hyperbolic arc-cosine for binary16 value.
pub fn cr_acoshf16(x: f16) -> f16 {
    let v = f32::from(x);
    let u = v.to_bits();
    let au = u & 0x7fffffff;

    if au >= 0x7f800000 {
        // NaN or Inf
        // acosh(+Inf) = +Inf, otherwise we get qNaN
        if (u == 0x7f800000) || ((au & 0x7fffff) != 0) {
            // +Inf or NaN
            return x + x;
        }
        return snanf16(); // -Inf
    }

    if (u >= 0x80000000) || (au <= 0x3f800000) {
        // x <= 1
        if u == 0x3f800000 {
            // x = 1
            return 0.0; // for x=1 the code below yields -0
        }
        // for x=1 the code below yields -0
        return snanf16();
    }

    let i = (au >> 23) - 127; // 2^i <= x < 2^(i+1)
    let t = v;
    let tt = t * t;
    let p = P[i as usize];
    let c4 = f32::mul_add(p[5], t, p[4]);
    let mut c2 = f32::mul_add(p[3], t, p[2]);
    c2 = f32::mul_add(c4, tt, c2);
    let c0 = f32::mul_add(p[1], t, p[0]);
    let mut y = f32::mul_add(c2, tt, c0);
    if i < 5 {
        y *= (t - 1.0).sqrt();
    }

    // deal with exceptions
    match (i, u) {
        (1, 0x40422000) => {
            y = h!("0x1.c63ffcp+0"); // x = 0x1.844p+1
        }
        (3, 0x411c0000) => {
            y = h!("0x1.7be006p+1"); // x = 0x1.38p+3
        }
        (3, 0x415f6000) => {
            y = h!("0x1.aa0002p+1"); // x = 0x1.becp+3
        }
        (4, 0x41bfe000) => {
            y = h!("0x1.ef5fecp+1"); // x = 0x1.7fcp+4
        }
        (4, 0x41c04000) => {
            y = h!("0x1.ef9ff6p+1"); // x = 0x1.808p+4
        }
        (4, 0x41dbc000) => {
            y = h!("0x1.006016p+2"); // x = 0x1.b78p+4
        }
        (4, 0x41d3c000) => {
            y = h!("0x1.fc002ap+1"); // x = 0x1.a78p+4
        }
        (5, 0x42374000) => {
            y = h!("0x1.21201cp+2"); // x = 0x1.6e8p+5
        }
        (5, 0x4245c000) => {
            y = h!("0x1.260012p+2"); // x = 0x1.8b8p+5
        }
        (6, 0x42890000) => {
            y = h!("0x1.3ae018p+2"); // x = 0x1.12p+6
        }
        (6, 0x429d6000) => {
            y = h!("0x1.43bff4p+2"); // x = 0x1.3acp+6
        }
        (7, 0x433d2000) => {
            y = h!("0x1.7be006p+2"); // x = 0x1.7a4p+7
        }
        (8, 0x43884000) => {
            y = h!("0x1.934004p+2"); // x = 0x1.108p+8
        }
        (8, 0x43f9a000) => {
            y = h!("0x1.ba000ep+2"); // x = 0x1.f34p+8
        }
        (9, 0x443b0000) => {
            y = h!("0x1.d3e00cp+2"); // x = 0x1.76p+9
        }
        (12, 0x45b78000) => {
            y = h!("0x1.2be008p+3"); // x = 0x1.6fp+12
        }
        (13, 0x4673a000) => {
            y = h!("0x1.4b2008p+3"); // x = 0x1.e74p+13
        }
        _ => (),
    }

    y as f16
}

#[cfg(test)]
mod tests {
    #[test]
    fn exhaustive() {
        for b in 0..=u16::MAX {
            let x = f16::from_bits(b);
            let y1 = super::cr_acoshf16(x);
            let y2 = core_math::acoshf16(x);
            assert_eq!(y1.to_bits(), y2.to_bits(), "acoshf16({x}) = {y1} vs {y2}");
        }
    }
}
