import pandas as pd
import pytest

from banyan_sigma import cli


def _targets():
    return pd.DataFrame(
        {
            "NAME": ["target-a", "target-b"],
            "RA": [10.0, 20.0],
            "DEC": [-10.0, -20.0],
            "PMRA": [100.0, 200.0],
            "PMDEC": [-100.0, -200.0],
            "EPMRA": [1.0, 2.0],
            "EPMDEC": [1.5, 2.5],
        }
    )


def test_association_results_uses_strict_absolute_probability_threshold():
    targets = _targets()
    output = pd.DataFrame(
        {
            ("ALL", "ABDMG"): [0.20, 0.01],
            ("ALL", "BPMG"): [0.03, 0.80],
            ("ALL", "FIELD"): [0.77, 0.19],
        }
    )

    results = cli.association_results(targets, output)

    assert results[["INPUT_ROW", "ASSOCIATION"]].to_dict("records") == [
        {"INPUT_ROW": 1, "ASSOCIATION": "ABDMG"},
        {"INPUT_ROW": 1, "ASSOCIATION": "BPMG"},
        {"INPUT_ROW": 2, "ASSOCIATION": "BPMG"},
    ]
    assert results["PROBABILITY"].tolist() == [0.20, 0.03, 0.80]
    assert results["NAME"].tolist() == ["target-a", "target-a", "target-b"]
    assert "FIELD" not in results["ASSOCIATION"].tolist()


def test_association_results_writes_headers_when_no_association_passes():
    targets = _targets().iloc[:1]
    output = pd.DataFrame(
        {
            ("ALL", "ABDMG"): [0.01],
            ("ALL", "FIELD"): [0.99],
        }
    )

    results = cli.association_results(targets, output)

    assert results.empty
    assert results.columns.tolist() == [
        "INPUT_ROW",
        *targets.columns,
        "ASSOCIATION",
        "PROBABILITY",
    ]


def test_classify_csv_detects_optional_columns_and_writes_results(tmp_path, monkeypatch):
    targets = _targets().iloc[:1].assign(PLX=25.0, EPLX=0.2)
    input_path = tmp_path / "targets.csv"
    output_path = tmp_path / "results.csv"
    targets.to_csv(input_path, index=False)
    call = {}

    def fake_membership_probability(**kwargs):
        call.update(kwargs)
        return pd.DataFrame(
            {
                ("ALL", "ABDMG"): [0.02],
                ("ALL", "FIELD"): [0.98],
            }
        )

    monkeypatch.setattr(cli, "membership_probability", fake_membership_probability)

    results = cli.classify_csv(input_path, output_path)

    assert call["column_names"]["PLX"] == "PLX"
    assert call["use_plx"] is True
    assert call["use_rv"] is False
    assert call["use_dist"] is False
    assert results["ASSOCIATION"].tolist() == ["ABDMG"]
    pd.testing.assert_frame_equal(pd.read_csv(output_path), results)


def test_column_names_rejects_incomplete_optional_measurement_pair():
    targets = _targets().assign(RV=5.0)

    with pytest.raises(ValueError, match="RV requires companion column ERV"):
        cli._column_names(targets)


def test_main_rejects_probability_outside_unit_interval():
    with pytest.raises(SystemExit):
        cli.main(["input.csv", "output.csv", "--min-probability", "1.1"])
