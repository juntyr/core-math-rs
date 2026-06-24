/* Correctly-rounded arc-cosine for binary16 value.

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

#![expect(clippy::approx_constant)]

use hexf::hexf32 as h;

use super::utils::snanf16;

// the following polynomials were generated using Sollya (cf acos.sollya)

/* degree-4 minimax polynomial for acos(x) over [0,0.25], with relative error
   bounded by 2^-22.943, manually optimized to reduce the number of exceptions
*/
const P0: [f32; 5] = [
    h!("0x1.921fb4p0"),
    h!("-0x1.fffb44p-1"),
    h!("-0x1.25e6b8p-10"),
    h!("-0x1.3cc114p-3"),
    h!("-0x1.a85b22p-5"),
];

/* degree-4 minimax polynomial for acos(x) over [0.25,0.5], with relative error
   bounded by 2^-20.789, manually optimized to reduce the number of exceptions
*/
const P1: [f32; 5] = [
    h!("0x1.91b678p0"),
    h!("-0x1.f515cap-1"),
    h!("-0x1.bd043ap-4"),
    h!("0x1.7e2d5ap-4"),
    h!("-0x1.190806p-2"),
];

/* degree-4 minimax polynomial for acos(x)/sqrt(1-x) over [0.5,1],
with relative error bounded by 2^-23.583 */
const P2: [f32; 5] = [
    h!("0x1.91fa1cp0"),
    h!("-0x1.ae5c5ep-3"),
    h!("0x1.31640cp-4"),
    h!("-0x1.98038p-6"),
    h!("0x1.251b5p-8"),
];

pub fn cr_acosf16(x: f16) -> f16 {
    let v = f32::from(x);
    let u = v.to_bits();
    let au = u & 0x7fffffff;

    if au >= 0x3f800000 {
        // NaN, Inf, or |x| >= 1
        if au == 0x3f800000 {
            if u == 0x3f800000 {
                return 0.0;
            }
            return h!("0x1.921fb6p+1") as f16; // TODO
        }

        if (au >> 23) == 0x3ff && ((au & 0x7fffff) != 0) {
            // qNaN or sNaN
            return x;
        }
        return snanf16(); // will signal invalid and return sNaN
    }

    let mut t = v;
    let tt = t * t;

    if (u >> 31) != 0 {
        // x < 0
        t = -t;
    }

    let c1;
    let c3;
    let mut y;
    if au < 0x3e800000 {
        // |x| < 0.25
        c1 = f32::mul_add(P0[2], t, P0[1]);
        c3 = f32::mul_add(P0[4], t, P0[3]);
        y = f32::mul_add(c3, tt, c1);
        y = f32::mul_add(y, t, P0[0]);
        if (u >> 31) != 0
        // x < 0
        {
            y = h!("0x1.921fb4p+1") - y;
            /* below we deal with a few exceptions that we were unable to remove
            by tuning the coefficients of p0 */
            if (0x36960000..=0x369a0000).contains(&au) {
                y = h!("0x1.922002p+0"); // x = -0x1.2cp-18 or -0x1.3p-18 or -0x1.34p-18
            }
            if u == 0xbcf80000 {
                y = h!("0x1.99e002p+0"); // x = -0x1.fp-6
            }
            if u == 0xbcfc0000 {
                y = h!("0x1.9a0006p+0"); // -0x1.f8p-6
            }
        }
    } else if au < 0x3f000000 {
        // 0.25 <= |x| < 0.5
        c1 = f32::mul_add(P1[2], t, P1[1]);
        c3 = f32::mul_add(P1[4], t, P1[3]);
        y = f32::mul_add(c3, tt, c1);
        y = f32::mul_add(y, t, P1[0]);
        if (u >> 31) != 0 {
            // x < 0
            y = h!("0x1.921fb4p+1") - y;
        }
    } else {
        // 0.5 <= |x| <= 1
        c1 = f32::mul_add(P2[2], t, P2[1]);
        c3 = f32::mul_add(P2[4], t, P2[3]);
        y = f32::mul_add(c3, tt, c1);
        y = f32::mul_add(y, t, P2[0]);
        y *= (1.0 - t).sqrt();
        if (u >> 31) != 0 {
            // x < 0
            y = h!("0x1.921fb4p+1") - y;
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
            let y1 = super::cr_acosf16(x);
            let y2 = core_math::acosf16(x);
            assert_eq!(y1.to_bits(), y2.to_bits(), "acosf16({x}) = {y1} vs {y2}");
        }
    }
}
