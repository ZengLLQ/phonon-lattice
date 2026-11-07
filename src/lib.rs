//! # PhononLattice
//!
//! High-performance crystal lattice dynamics, phonon dispersion, non-analytic corrections,
//! and anharmonic thermal transport in Rust.
//!
//! ## Core Features
//! - **Crystallography & Reciprocal Space:** Unit cells, lattice vectors, metric tensors, fractional-to-Cartesian mappings, and Brillouin zone paths.
//! - **Hermitian Eigensolver:** Self-contained, quadratic-convergent Jacobi complex diagonalization.
//! - **Harmonic Lattice Dynamics:** Interatomic Force Constant (IFC) tensors, exact Acoustic Sum Rule (ASR) invariance, and dynamical matrix construction.
//! - **Polar Non-Analytic Corrections:** Long-range macroscopic electric field coupling and Born effective charge tensors resolving LO-TO splitting.
//! - **Phonon Density of States (DOS):** Total and atom-projected DOS via uniform Monkhorst-Pack reciprocal grids and Gaussian broadening.
//! - **Quantum Harmonic Thermodynamics:** Vibrational zero-point energy, Helmholtz free energy, entropy, and heat capacity $C_v(T)$ reproducing Debye $T^3$ and Dulong-Petit limits.
//! - **Anharmonic Thermal Transport:** Boltzmann Transport Equation (BTE) solver under Single-Mode Relaxation Time Approximation (SMRTA) with boundary, isotope (Tamura), and 3-phonon Umklapp scattering.

#![allow(clippy::needless_range_loop, clippy::approx_constant)]

pub mod crystal;
pub mod dynamics;
pub mod materials;
pub mod math;
pub mod thermodynamics;
pub mod transport;

pub use crystal::{Atom, HighSymmetryPath, ReciprocalLattice, UnitCell};
pub use dynamics::{
    DensityOfStates, DispersionBand, DynamicalMatrixBuilder, InteratomicForceConstants,
    PhononDispersion, PhononDosCalculator,
};
pub use materials::{GaAsMaterial, GrapheneMaterial, NaClMaterial, SiliconMaterial};
pub use math::{
    eigensystem_hermitian, gaussian_delta, integrate_trapezoidal, EigenResult, Matrix3, Vector3,
};
pub use thermodynamics::{
    DebyeModel, DebyeProperties, ThermodynamicCalculator, ThermodynamicState,
};
pub use transport::{
    PhononScatteringModel, ScatteringRates, ThermalConductivityResult, ThermalConductivitySolver,
};
