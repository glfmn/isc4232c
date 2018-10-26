extern crate gnuplot;
extern crate num_complex;

use num_complex::Complex64;

struct Newton {
    x_n: Complex64,
}

impl Newton {
    fn new(x_0: Complex64) -> Self {
        Newton { x_n: x_0 }
    }
}

impl Iterator for Newton {
    type Item = Complex64;

    /// Update a singe time-step, checking the tolerance to see
    /// if convergence has happened
    fn next(&mut self) -> Option<Self::Item> {
        let x_n = self.x_n
                - (self.x_n*self.x_n*self.x_n - 1.)
                  /(3.*self.x_n*self.x_n);
        if ((x_n - self.x_n) / x_n).norm() < 1e-7 {
            None
        } else {
            self.x_n = x_n;
            Some(self.x_n)
        }
    }
}

fn main() {
    // Set the value of the three roots and pre-allocate a vector
    // to store the z values
    let z_1 = Complex64::new(1.,0.);
    let z_2 = Complex64::new(-0.5, -(3f64.sqrt()) / 2.);
    let z_3 = Complex64::new(-0.5, 3f64.sqrt()/2.);
    let mut zs: Vec<u8> = Vec::with_capacity(16000);
    for y in -200..200 {
        for x in -200..200 {
            // Prevent divide-by-zero
            if x == 0 && y == 0 {
                zs.push(0);
                continue
            }
            let mut z = Complex64::new(0.,0.);
            for z_n in Newton::new(Complex64 {
                re: x as f64/100.,
                im: y as f64/100.,
            }) {
                z = z_n;
            }
            if (z - z_1).norm() < 1e-6 {
                zs.push(0);
            } else if (z - z_2).norm() < 1e-6 {
                zs.push(1);
            } else if (z - z_3).norm() < 1e-6 {
                zs.push(2);
            }
        }
    }

    use gnuplot::{Figure, AxesCommon};
    use gnuplot::ContourStyle::*;
    use gnuplot::AutoOption::*;

    let mut fg = Figure::new();
    fg.set_terminal("epscairo", "target/attraction.eps")
        .axes3d()
        .set_title("Basins of Attraction", &[])
        .surface(zs, 400, 400, Some((-2.,-2.,2.,2.)), &[])
        .set_view(0., 0.)
        .set_z_grid(false)
        .set_view_map()
        .show_contours_custom(
            true, false, Cubic(3), Auto, &[0, 1, 2]
        );
    fg.show();
}
