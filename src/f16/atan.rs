/* Correctly-rounded arc-tangent for binary16 value.

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

// the following polynomials were generated using Sollya (cf atan.sollya)

/* Degree-7 minimax polynomial for atan(x) over [0,0.25], with relative error
   bounded by 2^-25.419, with coefficients of odd degree only, and degree-1
   coefficient forced to 1. Coefficients were later optimized to reduce the
   number of exceptions.
*/
const P0: [f32; 4] = [
    h!("0x1.fffffcp-1"),
    h!("-0x1.55546cp-2"),
    h!("0x1.98d0ep-3"),
    h!("-0x1.0c7c54p-3"),
];

/* degree-4 minimax polynomial for atan(x) over [0.25,0.5], with relative error
bounded by 2^-22.573 */
const P1: [f32; 5] = [
    h!("0x1.411612p-14"),
    h!("0x1.ff076cp-1"),
    h!("0x1.1ee64cp-6"),
    h!("-0x1.a6fc96p-2"),
    h!("0x1.81dd3ep-3"),
];

/* degree-4 minimax polynomial for atan(x) over [0.5,0.75], with relative error
bounded by 2^-21.757 */
const P2: [f32; 5] = [
    h!("-0x1.95964cp-8"),
    h!("0x1.0bad76p+0"),
    h!("-0x1.e652e2p-4"),
    h!("-0x1.e70ab4p-3"),
    h!("0x1.a5eb88p-4"),
];

/* degree-4 minimax polynomial for atan(x) over [0.75,1], with relative error
bounded by 2^-23.027 */
const P3: [f32; 5] = [
    h!("-0x1.f75b16p-6"),
    h!("0x1.2d3d1p+0"),
    h!("-0x1.882376p-2"),
    h!("0x1.d18c96p-13"),
    h!("0x1.6aa268p-6"),
];

/// Correctly-rounded arc-tangent for binary16 value.
pub fn cr_atanf16(x: f16) -> f16 {
    let v = f32::from(x);
    let u = v.to_bits();
    let au = u & 0x7fffffff;

    const HALF_PI: f32 = h!("0x1.921fb6p+0");

    if au >= 0x7f800000 {
        // NaN or Inf
        if au == 0x7f800000 {
            // +/-Inf
            return if u == 0x7f800000 { HALF_PI } else { -HALF_PI } as f16;
        }
        return x + x; // will signal invalid for sNaN and return qNaN
    }

    let mut t = v;

    // for x < 0 we use atan(-x) = -atan(x)
    let neg = (u >> 31) != 0;
    const NEG: [f32; 2] = [1.0, -1.0];
    let s = NEG[usize::from(neg)];
    t *= s;

    // now t >= 0

    // for x > 1 we use atan(x) = pi/2 - atan(1/x)
    let reduce = au > 0x3f800000;
    if reduce {
        t = 1.0 / t;
    }

    // now 0 <= t <= 1

    let tt = t * t;
    let mut y;
    if t <= 0.25 {
        // for |x| < 0x1.d14p-6, atan(x) rounds to x to nearest
        if !reduce && (t <= h!("0x1.d14p-6")) {
            if au == 0 {
                return x; // x = 0
            }
            t *= s;
            return f32::mul_add(t, h!("-0x1p-23"), t) as f16;
        }

        // deal with exceptional cases
        if au == 0x3e56a000 {
            // |x| = 0x1.ad4p-3
            return if au == u {
                h!("0x1.a72002p-3")
            } else {
                h!("-0x1.a72002p-3")
            } as f16;
        }
        if au == 0x4115c000 {
            // |x| = 0x1.2b8p+3
            return if au == u {
                h!("0x1.76dffep+0")
            } else {
                h!("-0x1.76dffep+0")
            } as f16;
        }
        if au == 0x42c32000 {
            // |x| = 0x1.864p+6
            return if au == u {
                h!("0x1.8f7ffep+0")
            } else {
                h!("-0x1.8f7ffep+0")
            } as f16;
        }

        let p = P0;
        let c5 = f32::mul_add(p[3], tt, p[2]);
        let c1 = f32::mul_add(p[1], tt, p[0]);
        let c1 = f32::mul_add(c5, tt * tt, c1);
        y = t * c1;
    } else {
        let p = if t <= 0.5 {
            P1
        } else if t <= 0.75 {
            P2
        } else {
            P3
        };
        let c3 = f32::mul_add(p[4], t, p[3]);
        let c2 = f32::mul_add(c3, t, p[2]);
        y = f32::mul_add(p[1], t, p[0]);
        y = f32::mul_add(c2, tt, y);
    }

    if reduce {
        y = HALF_PI - y; // argument reconstruction
    }

    if neg {
        y = -y;
    }

    y as f16
}

#[cfg(test)]
mod tests {
    #[test]
    fn exhaustive() {
        for b in 0..=u16::MAX {
            let x = f16::from_bits(b);
            let y1 = super::cr_atanf16(x);
            let y2 = core_math::atanf16(x);
            assert_eq!(
                y1.to_bits(),
                y2.to_bits(),
                "atanf16({x} @ {b:#04x}) = ({y1} @ {y1b:#04x}) vs ({y2} @ {y2b:#04x})",
                y1b = y1.to_bits(),
                y2b = y2.to_bits()
            );
        }
    }
}
