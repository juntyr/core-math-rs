/* Correctly-rounded arc-sine for binary16 value.

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

// the following polynomials were generated using Sollya (cf asin.sollya)

/* degree-5 minimax polynomial for asin(x) over [0,0.25], with relative error
   bounded by 2^-20.817, with coefficients of odd degree only, and degree-1
   coefficient forced to 1
*/
const P0: [f32; 3] = [1.0, h!("0x1.552e8ep-3"), h!("0x1.43696ep-4")];

/* degree-7 minimax polynomial for asin(x) over [0.25,0.5], with relative error
   bounded by 2^-19.202, with coefficients of odd degree only, and degree-1
   coefficient forced to 1
*/
const P1: [f32; 4] = [
    1.0,
    h!("0x1.55a7ep-3"),
    h!("0x1.258e54p-4"),
    h!("0x1.08e44p-4"),
];

/// Correctly-rounded arc-sine for binary16 value.
pub fn cr_asinf16(x: f16) -> f16 {
    let v = f32::from(x);
    let u = v.to_bits();
    let mut au = u & 0x7fffffff;

    const HALF_PI: f32 = h!("0x1.921fb6p+0");

    if au >= 0x3f800000 {
        // NaN, Inf, or |x| >= 1
        if au == 0x3f800000 {
            // |x| = 1
            return if u == 0x3f800000 { HALF_PI } else { -HALF_PI } as f16;
        }
        if (au >> 23) == 0x3ff && ((au & 0x7fffff) != 0) {
            // qNaN or sNaN
            return x + x;
        }
        return snanf16(); // will signal invalid and return sNaN
    }

    let mut t = v;

    /* for |x| < 0.5 we don't reduce t into -t, since the polynomials
    P0 and P1 are valid also on the negative side */

    let reduce = au >= 0x3f000000; // |x| >= 0.5
    if reduce {
        if (u >> 31) != 0 {
            // x < 0
            t = -t;
        }
        // argument reduction: asin(x) = pi/2 - 2*asin(sqrt((1-x)/2))
        t = ((1.0 - t) * 0.5).sqrt();
        au = t.to_bits() & 0x7fffffff;
    }

    // now 0 <= t <= 0.5
    let tt = t * t;

    let mut y;
    if au < 0x3e800000 {
        // 0 <= t < 0.25
        /* For |x| <= 0x1.71p-5, asin(x) rounds to x to nearest,
        we deal with that case separately, so that for x subnormal
        and a power of two, we get an underflow. */
        if (au <= 0x3d388000) && !reduce {
            if au == 0 {
                return x;
            }
            return if au == u {
                v + h!("0x1p-26")
            } else {
                v - h!("0x1p-26")
            } as f16;
        }
        if au == 0x3dd30000 {
            // |x| = 0x1.a6p-4
            return if au == u {
                h!("0x1.a6c00ap-4")
            } else {
                h!("-0x1.a6c00ap-4")
            } as f16;
        }
        if au == 0x3d688000 {
            // |x| = 0x1.d1p-5
            return if au == u {
                h!("0x1.d14004p-5")
            } else {
                h!("-0x1.d14004p-5")
            } as f16;
        }
        if au == 0x3dfa0000 {
            // |x| = 0x1.f4p-4
            return if au == u {
                h!("0x1.f5400ap-4")
            } else {
                h!("-0x1.f5400ap-4")
            } as f16;
        }
        /* Warning for rounding toward -Inf: let P0(t) = t*q(t). If we first
        compute q(t) and then multiply by t, for tiny t and rounding we will
        get q(t)=1, and then t, whereas the correct result is nextbelow(t). */
        let c1 = f32::mul_add(P0[2], tt, P0[1]);
        y = f32::mul_add(c1, tt * t, t);
    } else {
        // 0.25 <= t <= 0.5
        if au == 0x3eb24000 {
            // |x| = 0x1.648p-2
            return if au == u {
                h!("0x1.6c2012p-2")
            } else {
                h!("-0x1.6c2012p-2")
            } as f16;
        }
        if au == 0x3ed96000 {
            // |x| = 0x1.b2cp-2
            return if au == u {
                h!("0x1.c0fffap-2")
            } else {
                h!("-0x1.c0fffap-2")
            } as f16;
        }
        if au == 0x3ef0a000 {
            // |x| = 0x1.e14p-2
            return if au == u {
                h!("0x1.f4fffp-2")
            } else {
                h!("-0x1.f4fffp-2")
            } as f16;
        }
        let c5 = f32::mul_add(P1[3], tt, P1[2]);
        let c1 = f32::mul_add(P1[1], tt, P1[0]);
        y = f32::mul_add(c5, tt * tt, c1);
        y *= t;
    }

    if reduce {
        // argument reconstruction
        y = HALF_PI - 2.0 * y;
        if (u >> 31) != 0 {
            // x < 0
            y = -y;
        }
    }

    y as f16
}

#[cfg(test)]
mod tests {
    #[test]
    fn exhaustive() {
        for b in 0..=u16::MAX {
            let x = f16::from_bits(b);
            let y1 = super::cr_asinf16(x);
            let y2 = core_math::asinf16(x);
            assert_eq!(y1.to_bits(), y2.to_bits(), "asinf16({x}) = {y1} vs {y2}");
        }
    }
}
