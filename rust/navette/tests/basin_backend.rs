// SPDX-License-Identifier: LGPL-3.0-or-later
use navette::smatrix::synthesis::optimizer::{OptimizerBackend, available_backends};

#[test]
fn basin_names_and_feature_discovery() {
    for name in ["basin_lm_qr", "basin_trf"] {
        let backend = OptimizerBackend::parse(name).unwrap();
        assert_eq!(backend.as_str(), name);
        assert_eq!(backend.feature(), Some("opt-basin"));
        assert_eq!(backend.is_available(), cfg!(feature = "opt-basin"));
        assert_eq!(available_backends().contains(&name), backend.is_available());
        assert_eq!(backend.bounds_are_native(), name == "basin_trf");
        assert!(!backend.runs_to_max_iterations());
    }
}

#[cfg(not(feature = "opt-basin"))]
#[test]
fn missing_basin_reports_the_build_feature() {
    use navette::smatrix::synthesis::{optimizer::run_optimizer, thick_opt::*};
    for name in ["basin_lm_qr", "basin_trf"] {
        let cfg = LmConfig {
            backend: OptimizerBackend::parse(name).unwrap(),
            ..Default::default()
        };
        let err = run_optimizer(
            &|_, r| {
                *r = vec![1.0];
                Ok(())
            },
            None::<&NoJacobian>,
            &[0.5],
            &[0.0],
            &[1.0],
            &cfg,
        )
        .unwrap_err();
        assert!(err.contains("--features opt-basin"), "{err}");
    }
}

#[cfg(feature = "opt-basin")]
mod enabled {
    use super::*;
    use navette::smatrix::synthesis::{optimizer::run_optimizer, thick_opt::*};
    use std::sync::atomic::{AtomicUsize, Ordering};

    struct Jac<F>(F);
    impl<F: Fn(&[f64], &mut Vec<f64>) -> Result<Option<usize>, String> + Sync> JacobianSource
        for Jac<F>
    {
        fn fill(&self, x: &[f64], j: &mut Vec<f64>) -> Result<Option<usize>, String> {
            (self.0)(x, j)
        }
    }

    fn config(name: &str) -> LmConfig {
        LmConfig {
            backend: OptimizerBackend::parse(name).unwrap(),
            ..Default::default()
        }
    }

    #[test]
    fn interior_solution_cost_and_evaluation_accounting() {
        for name in ["basin_lm_qr", "basin_trf"] {
            for mode in 0..3 {
                let calls = AtomicUsize::new(0);
                let residual = |x: &[f64], r: &mut Vec<f64>| {
                    calls.fetch_add(1, Ordering::Relaxed);
                    *r = vec![x[0] - 1.0, 2.0 * (x[1] + 2.0), 3.0];
                    Ok(())
                };
                let jac = Jac(|_: &[f64], j: &mut Vec<f64>| {
                    if mode == 2 {
                        return Ok(None);
                    }
                    *j = vec![1.0, 0.0, 0.0, 2.0, 0.0, 0.0];
                    Ok(Some(3))
                });
                let result = run_optimizer(
                    &residual,
                    (mode != 0).then_some(&jac),
                    &[0.0, 0.0],
                    &[-5.0; 2],
                    &[5.0; 2],
                    &config(name),
                )
                .unwrap();
                assert!((result.x[0] - 1.0).abs() < 1e-6, "{name}: {result:?}");
                assert!((result.x[1] + 2.0).abs() < 1e-6, "{name}: {result:?}");
                assert!((result.cost - 9.0).abs() < 1e-10);
                assert_eq!(result.evals, calls.load(Ordering::Relaxed));
                assert_eq!(result.analytic_jacobians > 0, mode == 1);
                assert_eq!(result.backend.as_str(), name);
                assert!(result.gain_ratio.is_nan());
                assert!(!matches!(
                    result.termination,
                    LmTermination::MaxIterations | LmTermination::Stalled
                ));
            }
        }
    }

    #[test]
    fn trf_boundary_matrix_and_feasible_difference_probes() {
        // Interior, lower, upper, corner, boundary starts, and a fixed coordinate.
        let cases = [
            ([0.3, 0.7], [0.5, 0.5], [0.0, 0.0], [1.0, 1.0]),
            ([-1.0, 0.7], [0.5, 0.5], [0.0, 0.0], [1.0, 1.0]),
            ([0.3, 2.0], [0.5, 0.5], [0.0, 0.0], [1.0, 1.0]),
            ([-1.0, 2.0], [0.5, 0.5], [0.0, 0.0], [1.0, 1.0]),
            ([0.3, 0.7], [0.0, 1.0], [0.0, 0.0], [1.0, 1.0]),
            ([0.3, 0.7], [0.5, 0.5], [0.5, 0.0], [0.5, 1.0]),
        ];
        for (target, start, lower, upper) in cases {
            for analytic in [false, true] {
                let residual = |x: &[f64], r: &mut Vec<f64>| {
                    assert!(
                        x.iter()
                            .enumerate()
                            .all(|(i, &v)| v >= lower[i] && v <= upper[i])
                    );
                    *r = vec![x[0] - target[0], x[1] - target[1], 0.0];
                    Ok(())
                };
                let jac = Jac(|_: &[f64], j: &mut Vec<f64>| {
                    *j = vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
                    Ok(Some(3))
                });
                let result = run_optimizer(
                    &residual,
                    analytic.then_some(&jac),
                    &start,
                    &lower,
                    &upper,
                    &config("basin_trf"),
                )
                .unwrap_or_else(|e| panic!("target={target:?}, start={start:?}, lower={lower:?}, upper={upper:?}, analytic={analytic}: {e}"));
                let expected: Vec<f64> = (0..2)
                    .map(|i| target[i].clamp(lower[i], upper[i]))
                    .collect();
                for (&got, &want) in result.x.iter().zip(&expected) {
                    assert!((got - want).abs() < 1e-6, "{target:?}: {result:?}");
                }
                let cost = (expected[0] - target[0]).powi(2) + (expected[1] - target[1]).powi(2);
                assert!((result.cost - cost).abs() < 1e-8);
            }
        }
    }

    #[test]
    fn trf_all_fixed_needs_only_one_residual() {
        let jac = Jac(
            |_: &[f64], _: &mut Vec<f64>| -> Result<Option<usize>, String> {
                panic!("fixed problem needs no Jacobian")
            },
        );
        let result = run_optimizer(
            &|_, r| {
                *r = vec![2.0];
                Ok(())
            },
            Some(&jac),
            &[0.5],
            &[0.5],
            &[0.5],
            &config("basin_trf"),
        )
        .unwrap();
        assert_eq!(result.x, [0.5]);
        assert_eq!(result.cost, 4.0);
        assert_eq!(result.evals, 1);
        assert_eq!(format!("{:?}", result.termination), "Converged");
    }

    #[test]
    fn zero_columns_and_singular_jtj_do_not_require_unique_parameters() {
        for name in ["basin_lm_qr", "basin_trf"] {
            for analytic in [false, true] {
                let jac = Jac(|_: &[f64], j: &mut Vec<f64>| {
                    *j = vec![1.0, 1.0, 0.0, 2.0, 2.0, 0.0];
                    Ok(Some(2))
                });
                let result = run_optimizer(
                    &|x, r| {
                        *r = vec![x[0] + x[1] - 3.0, 2.0 * (x[0] + x[1] - 3.0)];
                        Ok(())
                    },
                    analytic.then_some(&jac),
                    &[0.5; 3],
                    &[-5.0; 3],
                    &[5.0; 3],
                    &config(name),
                )
                .unwrap();
                assert!(
                    (result.x[0] + result.x[1] - 3.0).abs() < 1e-6,
                    "{name}: {result:?}"
                );
                assert!((result.x[2] - 0.5).abs() < 1e-10);
                assert!(result.cost < 1e-10);
            }
        }
    }

    #[test]
    fn active_bounds_with_singular_jtj_and_a_zero_column() {
        for name in ["basin_lm_qr", "basin_trf"] {
            for analytic in [false, true] {
                let jac = Jac(|_: &[f64], j: &mut Vec<f64>| {
                    *j = vec![1.0, 1.0, 0.0, 2.0, 2.0, 0.0];
                    Ok(Some(2))
                });
                let result = run_optimizer(
                    &|x, r| {
                        *r = vec![x[0] + x[1] - 4.0, 2.0 * (x[0] + x[1] - 4.0)];
                        Ok(())
                    },
                    analytic.then_some(&jac),
                    &[0.5; 3],
                    &[0.0; 3],
                    &[1.0; 3],
                    &config(name),
                )
                .unwrap();
                assert!((result.x[0] - 1.0).abs() < 1e-6, "{name}: {result:?}");
                assert!((result.x[1] - 1.0).abs() < 1e-6, "{name}: {result:?}");
                assert!((result.x[2] - 0.5).abs() < 1e-10);
                assert!((result.cost - 20.0).abs() < 1e-8);
            }
        }
    }

    #[test]
    fn nonzero_floor_returns_an_accurate_stalled_result() {
        let residual = |x: &[f64], r: &mut Vec<f64>| {
            *r = vec![x[0] - 0.3, x[1] - 0.7, 2.0];
            Ok(())
        };
        let jac = Jac(|_: &[f64], j: &mut Vec<f64>| {
            *j = vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
            Ok(Some(3))
        });
        for analytic in [false, true] {
            let result = run_optimizer(
                &residual,
                analytic.then_some(&jac),
                &[0.0, 1.0],
                &[0.0; 2],
                &[1.0; 2],
                &config("basin_trf"),
            )
            .unwrap();
            assert_eq!(result.termination, LmTermination::Stalled, "{result:?}");
            assert!((result.cost - 4.0).abs() < 1e-12);
            assert!((result.x[0] - 0.3).abs() < 1e-7);
            assert!((result.x[1] - 0.7).abs() < 1e-7);
            assert_eq!(result.analytic_jacobians > 0, analytic);
        }
    }

    #[test]
    fn native_gradient_convergence_is_reported() {
        for name in ["basin_lm_qr", "basin_trf"] {
            for analytic in [false, true] {
                let jac = Jac(|_: &[f64], j: &mut Vec<f64>| {
                    *j = vec![1.0, 0.0];
                    Ok(Some(2))
                });
                let result = run_optimizer(
                    &|x, r| {
                        *r = vec![x[0], 2.0];
                        Ok(())
                    },
                    analytic.then_some(&jac),
                    &[0.0],
                    &[-1.0],
                    &[1.0],
                    &config(name),
                )
                .unwrap();
                assert_eq!(
                    result.termination,
                    LmTermination::Gradient,
                    "{name}: {result:?}"
                );
                assert_eq!(result.cost, 4.0);
                assert_eq!(result.x, [0.0]);
            }
        }
    }

    #[test]
    fn lm_native_cost_and_step_convergence_are_reported() {
        for (ftol, xtol, expected) in [
            (1.0, 0.0, LmTermination::Cost),
            (0.0, 1.0, LmTermination::Step),
        ] {
            let mut cfg = config("basin_lm_qr");
            cfg.gtol = 0.0;
            cfg.ftol = ftol;
            cfg.xtol = xtol;
            let result = run_optimizer(
                &|x, r| {
                    *r = vec![x[0] - 0.3, 2.0];
                    Ok(())
                },
                None::<&NoJacobian>,
                &[0.2],
                &[0.0],
                &[1.0],
                &cfg,
            )
            .unwrap();
            assert_eq!(result.termination, expected, "{result:?}");
        }
    }

    #[test]
    fn budgets_count_probes_and_stop_at_iteration_boundaries() {
        for name in ["basin_lm_qr", "basin_trf"] {
            let residual = |x: &[f64], r: &mut Vec<f64>| {
                *r = vec![x[0].exp() - 10.0];
                Ok(())
            };
            let mut cfg = config(name);
            cfg.max_evals = 1;
            let result =
                run_optimizer(&residual, None::<&NoJacobian>, &[0.1], &[0.0], &[5.0], &cfg)
                    .unwrap();
            assert_eq!(result.termination, LmTermination::MaxIterations);
            assert_eq!(result.iterations, 0);
            assert!(result.evals > cfg.max_evals);
            assert!((result.cost - (result.x[0].exp() - 10.0).powi(2)).abs() < 1e-12);
            cfg.max_evals = 100_000;
            cfg.max_iterations = 1;
            let result =
                run_optimizer(&residual, None::<&NoJacobian>, &[0.1], &[0.0], &[5.0], &cfg)
                    .unwrap();
            assert_eq!(result.termination, LmTermination::MaxIterations);
            assert_eq!(result.iterations, 1);
        }
    }

    #[test]
    fn invalid_inputs_are_rejected_before_callbacks() {
        let residual = |_: &[f64], _: &mut Vec<f64>| -> Result<(), String> {
            panic!("invalid inputs must not reach the problem")
        };
        for name in ["basin_lm_qr", "basin_trf"] {
            for (x, lower, upper) in [
                (vec![], vec![], vec![]),
                (vec![0.5], vec![], vec![1.0]),
                (vec![f64::NAN], vec![0.0], vec![1.0]),
                (vec![0.5], vec![1.0], vec![0.0]),
                (vec![0.5], vec![0.0], vec![f64::INFINITY]),
            ] {
                let err = run_optimizer(
                    &residual,
                    None::<&NoJacobian>,
                    &x,
                    &lower,
                    &upper,
                    &config(name),
                )
                .unwrap_err();
                assert!(err.starts_with(name), "{err}");
            }
            for tolerance in [f64::NAN, -1.0, f64::INFINITY] {
                let mut cfg = config(name);
                cfg.ftol = tolerance;
                assert!(
                    run_optimizer(&residual, None::<&NoJacobian>, &[0.5], &[0.0], &[1.0], &cfg)
                        .unwrap_err()
                        .contains("ftol")
                );
            }
        }
        assert!(
            run_optimizer(
                &residual,
                None::<&NoJacobian>,
                &[2.0],
                &[0.0],
                &[1.0],
                &config("basin_trf")
            )
            .unwrap_err()
            .contains("outside")
        );
    }

    #[test]
    fn changing_residual_lengths_and_nonfinite_derivatives_are_errors() {
        for name in ["basin_lm_qr", "basin_trf"] {
            let calls = AtomicUsize::new(0);
            let residual = |x: &[f64], r: &mut Vec<f64>| {
                *r = vec![x[0]; 1 + usize::from(calls.fetch_add(1, Ordering::Relaxed) > 0)];
                Ok(())
            };
            assert!(
                run_optimizer(
                    &residual,
                    None::<&NoJacobian>,
                    &[0.5],
                    &[0.0],
                    &[1.0],
                    &config(name)
                )
                .unwrap_err()
                .contains("residual length")
            );
            let jac = Jac(|_: &[f64], j: &mut Vec<f64>| {
                *j = vec![f64::NAN];
                Ok(Some(1))
            });
            assert!(
                run_optimizer(
                    &|x, r| {
                        *r = vec![x[0]];
                        Ok(())
                    },
                    Some(&jac),
                    &[0.5],
                    &[0.0],
                    &[1.0],
                    &config(name)
                )
                .unwrap_err()
                .contains("solver failed")
            );
        }
    }

    #[test]
    fn callback_errors_and_invalid_shapes_are_errors_not_panics() {
        for name in ["basin_lm_qr", "basin_trf"] {
            let cfg = config(name);
            let err = run_optimizer(
                &|_, _| Err("sentinel residual".into()),
                None::<&NoJacobian>,
                &[0.5],
                &[0.0],
                &[1.0],
                &cfg,
            )
            .unwrap_err();
            assert!(
                err.contains(name) && err.contains("sentinel residual"),
                "{err}"
            );
            for rows in [Some(2), Some(1)] {
                let jac = Jac(move |_: &[f64], j: &mut Vec<f64>| {
                    j.clear();
                    Ok(rows)
                });
                let err = run_optimizer(
                    &|x, r| {
                        *r = vec![x[0]];
                        Ok(())
                    },
                    Some(&jac),
                    &[0.5],
                    &[0.0],
                    &[1.0],
                    &cfg,
                )
                .unwrap_err();
                assert!(err.contains("Jacobian"), "{err}");
            }
            let err = run_optimizer(
                &|_, r| {
                    r.clear();
                    Ok(())
                },
                None::<&NoJacobian>,
                &[0.5],
                &[0.0],
                &[1.0],
                &cfg,
            )
            .unwrap_err();
            assert!(err.contains("empty residual"), "{err}");
            let jac = Jac(|_: &[f64], _: &mut Vec<f64>| Err("sentinel Jacobian".into()));
            let err = run_optimizer(
                &|x, r| {
                    *r = vec![x[0]];
                    Ok(())
                },
                Some(&jac),
                &[0.5],
                &[0.0],
                &[1.0],
                &cfg,
            )
            .unwrap_err();
            assert!(err.contains("sentinel Jacobian"), "{err}");
        }
    }
}
