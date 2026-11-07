use std::f64::consts::PI;

/// Gaussian delta approximation:
/// delta_sigma(x) = 1 / (sigma * sqrt(2 * pi)) * exp(-x^2 / (2 * sigma^2))
pub fn gaussian_delta(x: f64, sigma: f64) -> f64 {
    let norm = 1.0 / (sigma * (2.0 * PI).sqrt());
    let exponent = -0.5 * (x / sigma).powi(2);
    if exponent < -50.0 {
        0.0
    } else {
        norm * exponent.exp()
    }
}

/// Numerically integrates y(x) over a uniform grid [x_0, x_n] with step dx using the trapezoidal rule.
pub fn integrate_trapezoidal(y: &[f64], dx: f64) -> f64 {
    if y.len() < 2 {
        return 0.0;
    }
    let n = y.len();
    let mut sum = 0.5 * (y[0] + y[n - 1]);
    for i in 1..(n - 1) {
        sum += y[i];
    }
    sum * dx
}
