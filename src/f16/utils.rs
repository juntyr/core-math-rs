use std::hint::black_box;

pub fn snanf16() -> f16 {
    // signal invalid and return sNaN
    //
    // use black_box to force evaluation at runtime
    // and ensure equivalent behaviour with C
    black_box(0.0) / black_box(0.0)
}
