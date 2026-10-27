# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Regression cases for the standalone SciPy comparison harness."""

import numpy as np
import pytest

from validation.review.lm_check import thickness_deviation


def test_retained_negligible_films_do_not_change_the_parameter_count():
    actual = np.array([4e-6, 302.8989])
    expected = np.array([1.02720976e-7, 302.89896])
    assert thickness_deviation(actual, expected) == pytest.approx(6e-5)


def test_the_ci_basin_fd_result_is_within_the_existing_thickness_tolerance():
    actual = np.array([294.207866, 1.24931187e-4])
    expected = np.array([294.207836, 3.57098603e-5])
    assert thickness_deviation(actual, expected, zero_tol=1e-4) < 1e-4


@pytest.mark.parametrize("actual,expected", [
    ([302.89896], [1e-7, 302.89896]),
    ([1e-7, 302.89896], [302.89896]),
    ([], [1e-8]),
    ([], []),
])
def test_actual_removal_of_negligible_films_preserves_parity(actual, expected):
    assert thickness_deviation(np.array(actual), np.array(expected)) == 0.0


@pytest.mark.parametrize("actual,expected", [
    ([100.0], [100.0, 1.0]),
    ([0.0, 100.0], [100.0, 0.0]),
    ([294.207866, 2e-4], [294.207836, 3.57098603e-5]),
    ([float("nan")], []),
    ([float("inf")], [float("inf")]),
])
def test_disagreement_and_invalid_results_still_fail(actual, expected):
    assert thickness_deviation(
        np.array(actual), np.array(expected), zero_tol=1e-4) > 1e-4
