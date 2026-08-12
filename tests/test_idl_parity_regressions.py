import os

import numpy as np
from astropy.table import Table

import banyan_sigma
from banyan_sigma import membership_probability


def _models():
    path = os.path.join(
        os.path.dirname(banyan_sigma.__file__),
        "data",
        "banyan_sigma_parameters.fits",
    )
    models = Table.read(path)
    names = np.array([str(value).strip() for value in models["NAME"]])
    return models[np.isin(names, ["CAR", "FIELD"])]


def _models_with_mixture():
    path = os.path.join(
        os.path.dirname(banyan_sigma.__file__),
        "data",
        "banyan_sigma_parameters.fits",
    )
    models = Table.read(path)
    names = np.array([str(value).strip() for value in models["NAME"]])
    return models[np.isin(names, ["ABDMG", "FIELD"])]


def _calculate(**mode):
    return membership_probability(
        ra=[56.75],
        dec=[24.12],
        pmra=[42.3],
        pmdec=[-63.4],
        epmra=[0.2],
        epmdec=[0.2],
        custom_models=_models(),
        restrained_distance_range=np.array([[23.93], [75.0]]),
        use_plx=False,
        **mode,
    )


def test_measured_distance_does_not_overwrite_optimal_output():
    output = _calculate(
        use_rv=False,
        use_dist=True,
        dist=[39.5],
        edist=[0.6],
    )

    assert output[("CAR", "D_OPT")].iloc[0] != 39.5
    np.testing.assert_allclose(
        output[("CAR", "D_OPT")].iloc[0],
        42.30869969044697,
        atol=2e-6,
        rtol=0,
    )
    np.testing.assert_allclose(
        output[("CAR", "ED_OPT")].iloc[0],
        0.56699385,
        atol=2e-8,
        rtol=0,
    )


def test_measured_rv_does_not_overwrite_optimal_output():
    output = _calculate(
        use_rv=True,
        use_dist=False,
        rv=[15.2],
        erv=[0.4],
    )

    assert output[("CAR", "RV_OPT")].iloc[0] != 15.2
    np.testing.assert_allclose(
        output[("CAR", "RV_OPT")].iloc[0],
        13.28004508,
        atol=2e-7,
        rtol=0,
    )
    np.testing.assert_allclose(
        output[("CAR", "ERV_OPT")].iloc[0],
        0.32693565,
        atol=1e-8,
        rtol=0,
    )


def test_uvw_uncertainties_match_idl_analytical_propagation():
    output = _calculate(use_rv=False, use_dist=False)

    expected = np.array([0.51357009, 0.62149400, 0.27196058])
    actual = np.array(
        [
            output[("CAR", "EU")].iloc[0],
            output[("CAR", "EV")].iloc[0],
            output[("CAR", "EW")].iloc[0],
        ]
    )
    np.testing.assert_allclose(actual, expected, atol=2e-8, rtol=0)


def test_mixture_metrics_use_component_likelihood_weights():
    output = membership_probability(
        ra=[56.75],
        dec=[24.12],
        pmra=[42.3],
        pmdec=[-63.4],
        epmra=[0.2],
        epmdec=[0.2],
        custom_models=_models_with_mixture(),
        restrained_distance_range=np.array([[23.93], [75.0]]),
        use_rv=False,
        use_dist=False,
        use_plx=False,
    )

    np.testing.assert_allclose(
        output[("ABDMG", "D_OPT")].iloc[0],
        81.32006444,
        atol=2e-6,
        rtol=0,
    )
    np.testing.assert_allclose(
        output[("ABDMG", "EV")].iloc[0],
        0.52473393,
        atol=2e-8,
        rtol=0,
    )
