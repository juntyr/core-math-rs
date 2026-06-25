/* Correctly-rounded cubic root of binary16 value.

Copyright (c) 2025 Maxence Ponsardin.

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

pub fn cr_cbrtf16(x: f16) -> f16 {
    const TM: [u16; 16] = [
        0x00, 0xff, 0xff, 0xff, 0x33, 0x1c, 0x32, 0xff, 0xff, 0xff, 0xff, 0x00, 0xff, 0xff, 0xff,
        0x10,
    ];
    const TE: [u16; 16] = [0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0];
    const TF: [u16; 16] = [
        0x3c00, 0, 0, 0, 0x3d80, 0x3f00, 0x3c80, 0, 0, 0, 0, 0x3e00, 0, 0, 0, 0x3d00,
    ];

    let t = x;
    let mut tu = t.to_bits();
    let debug = tu == 0x7d5c;
    let xf = f32::from(x);
    let xu = xf.to_bits();
    if u32::from(TM[((xu >> 19) % 16) as usize]) == ((xu >> 13) % 64) {
        // exact cases (not supported by cbrtf)
        let expo = ((xu & 0x7fffffff) >> 23) as i32;
        if i32::from(TE[((xu >> 19) % 16) as usize]) == ((expo + 2) % 3) {
            if debug {
                println!("a");
            }
            tu = u32::from(tu & 0x8000)
                .wrapping_add(
                    (((expo - 127 - i32::from(TE[((xu >> 19) % 16) as usize])) / 3) as u32) << 10,
                )
                .wrapping_add(u32::from(TF[((xu >> 19) % 16) as usize])) as u16;
            return f16::from_bits(tu);
        }
    }
    if (tu & 0x03ff) == 0x0151 {
        // only wrong case is 0x1.544pk with k = 1 mod 3
        let expo = i32::from((tu & 0x7fff) >> 10);
        if ((expo % 3) == 1) && (expo < 31) {
            // avoid sNaN and k != 1 mod 3
            tu = ((((expo - 16) / 3 + 15) << 10) + 0x018b + i32::from((tu >> 15) << 15)) as u16;
            if debug {
                println!("b");
            }
            return (f32::from(f16::from_bits(tu)) + h!("0x1p-16") * (f32::from(tu >> 15) - 0.5))
                as f16;
        }
    }
    if debug {
        println!("c");
    }
    if debug {
        println!(
            "{x} {:#x} {xf} {xu} {} {}",
            x.to_bits(),
            xf.cbrt(),
            xf.cbrt() as f16
        );
    }
    // TODO: which cbrt should be called here?
    xf.cbrt() as f16
}

#[cfg(test)]
mod tests {
    #[test]
    fn exhaustive() {
        for b in 0..=u16::MAX {
            let x = f16::from_bits(b);
            let y1 = super::cr_cbrtf16(x);
            let y2 = core_math::cbrtf16(x);
            assert_eq!(
                y1.to_bits(),
                y2.to_bits(),
                "cbrtf16({x} @ {b:#04x}) = ({y1} @ {y1b:#04x}) vs ({y2} @ {y2b:#04x})",
                y1b = y1.to_bits(),
                y2b = y2.to_bits()
            );
        }
    }
}
