"""Unit tests for the module survey script in `tools/survey_module.py`."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOL_PATH = REPO_ROOT / "tools" / "survey_module.py"


def _load_surveyor():
    """*import the surveyor by path, since `tools/` is not an installed package*"""
    spec = importlib.util.spec_from_file_location("survey_module", TOOL_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


surveyor = _load_surveyor()


def _static_survey(tmp_path: Path, source: str):
    modulePath = tmp_path / "sample.py"
    modulePath.write_text(source)
    return surveyor.static_survey(modulePath)


def _hit_kinds(survey) -> list[str]:
    return [hit.kind for hit in survey.hits]


# LINE COUNT


def test_line_count_is_the_number_of_physical_lines(tmp_path: Path) -> None:
    # ARRANGE
    source = "import os\n\nvalue = 1\n"

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert survey.lineCount == 3


def test_a_file_without_a_trailing_newline_still_counts_its_last_line(tmp_path: Path) -> None:
    # ARRANGE
    source = "import os\nvalue = 1"

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert survey.lineCount == 2


# FUNCTION LENGTHS


def test_a_function_over_fifty_lines_is_reported_with_its_name_length_and_start_line(tmp_path: Path) -> None:
    # ARRANGE
    body = "\n".join(f"    value = {index}" for index in range(60))
    source = f"def produce_product(self):\n{body}\n"

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert [(function.name, function.length, function.line) for function in survey.functions] == [
        ("produce_product", 61, 1)
    ]


def test_a_function_of_exactly_fifty_lines_is_not_reported(tmp_path: Path) -> None:
    # ARRANGE
    body = "\n".join(f"    value = {index}" for index in range(49))
    source = f"def short(self):\n{body}\n"

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert survey.functions == []


def test_a_decorated_function_is_measured_from_its_def_line(tmp_path: Path) -> None:
    # ARRANGE
    body = "\n".join(f"    value = {index}" for index in range(60))
    source = f"@property\ndef produce_product(self):\n{body}\n"

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert survey.functions[0].line == 2
    assert survey.functions[0].length == 61


def test_a_method_is_reported_by_its_qualified_name(tmp_path: Path) -> None:
    # ARRANGE
    body = "\n".join(f"        value = {index}" for index in range(60))
    source = f"class soxs_mbias(object):\n    def produce_product(self):\n{body}\n"

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert survey.functions[0].name == "soxs_mbias.produce_product"


def test_functions_over_two_hundred_lines_are_separable_from_the_rest(tmp_path: Path) -> None:
    # ARRANGE
    longBody = "\n".join(f"    value = {index}" for index in range(250))
    shortBody = "\n".join(f"    value = {index}" for index in range(60))
    source = f"def huge(self):\n{longBody}\n\n\ndef large(self):\n{shortBody}\n"

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert [function.name for function in survey.functions if function.length > surveyor.UNSPLIT_REASON_LINES] == [
        "huge"
    ]


# DUPLICATION HITS


def test_a_naive_utcnow_call_is_a_duplication_hit(tmp_path: Path) -> None:
    # ARRANGE
    source = 'utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")\n'

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert _hit_kinds(survey) == ["naive-utcnow"]
    assert survey.hits[0].line == 1


def test_a_timezone_aware_now_call_is_not_a_duplication_hit(tmp_path: Path) -> None:
    # ARRANGE
    source = "stamp = datetime.now(timezone.utc)\n"

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert _hit_kinds(survey) == []


def test_a_qc_row_builder_is_a_duplication_hit(tmp_path: Path) -> None:
    # ARRANGE
    source = (
        "self.qc = pd.concat(\n"
        "    [\n"
        "        self.qc,\n"
        '        pd.DataFrame([{"qc_name": "STRUCTX", "qc_value": 1.0}]),\n'
        "    ],\n"
        "    ignore_index=True,\n"
        ")\n"
    )

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert _hit_kinds(survey) == ["qc-row-builder"]


def test_a_products_row_builder_is_its_own_duplication_kind(tmp_path: Path) -> None:
    # ARRANGE
    source = (
        "self.products = pd.concat(\n"
        "    [\n"
        "        self.products,\n"
        '        pd.DataFrame([{"product_label": "MBIAS"}]),\n'
        "    ],\n"
        "    ignore_index=True,\n"
        ")\n"
    )

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert _hit_kinds(survey) == ["products-row-builder"]


def test_the_legacy_series_to_frame_row_builder_is_still_a_duplication_hit(tmp_path: Path) -> None:
    # ARRANGE
    source = (
        "self.qc = pd.concat(\n"
        "    [\n"
        "        self.qc,\n"
        '        pd.Series({"qc_name": "STRUCTX", "qc_value": 1.0}).to_frame().T,\n'
        "    ],\n"
        "    ignore_index=True,\n"
        ")\n"
    )

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert _hit_kinds(survey) == ["qc-row-builder"]


def test_a_series_built_from_a_variable_is_not_a_row_builder(tmp_path: Path) -> None:
    # ARRANGE
    source = "self.qc = pd.concat([self.qc, pd.Series(qcValues).to_frame().T], ignore_index=True)\n"

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert _hit_kinds(survey) == []


def test_a_whole_table_concatenation_is_not_a_row_builder(tmp_path: Path) -> None:
    # ARRANGE
    source = "self.qc = pd.concat([self.qc, recipeQc], ignore_index=True)\n"

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert _hit_kinds(survey) == []


def test_a_multi_row_dataframe_concatenation_is_not_a_row_builder(tmp_path: Path) -> None:
    # ARRANGE
    source = 'self.qc = pd.concat([self.qc, pd.DataFrame([{"a": 1}, {"a": 2}])], ignore_index=True)\n'

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert _hit_kinds(survey) == []


def test_a_savefig_call_is_a_duplication_hit(tmp_path: Path) -> None:
    # ARRANGE
    source = 'plt.savefig(filePath, dpi=720)\n'

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert _hit_kinds(survey) == ["savefig"]


def test_duplication_hits_are_reported_in_source_order(tmp_path: Path) -> None:
    # ARRANGE
    source = (
        "plt.savefig(filePath)\n"
        "utcnow = datetime.utcnow()\n"
        'self.products = pd.concat([self.products, pd.DataFrame([{"a": 1}])], ignore_index=True)\n'
    )

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert _hit_kinds(survey) == ["savefig", "naive-utcnow", "products-row-builder"]


# BARE EXCEPTS


def test_a_bare_except_is_counted_with_its_line(tmp_path: Path) -> None:
    # ARRANGE
    source = "try:\n    run()\nexcept:\n    pass\n"

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert survey.bareExcepts == [3]


def test_a_typed_except_is_not_counted(tmp_path: Path) -> None:
    # ARRANGE
    source = "try:\n    run()\nexcept ValueError:\n    pass\n"

    # ACT
    survey = _static_survey(tmp_path, source)

    # ASSERT
    assert survey.bareExcepts == []


# TOOL ERRORS


def test_an_unparseable_module_raises_a_survey_error(tmp_path: Path) -> None:
    # ARRANGE
    modulePath = tmp_path / "broken.py"
    modulePath.write_text("def broken(:\n")

    # ACT / ASSERT
    with pytest.raises(surveyor.SurveyError):
        surveyor.static_survey(modulePath)


def test_a_missing_module_raises_a_survey_error(tmp_path: Path) -> None:
    # ACT / ASSERT
    with pytest.raises(surveyor.SurveyError):
        surveyor.static_survey(tmp_path / "absent.py")


# RUFF FINDINGS


def test_ruff_findings_are_counted_by_code(tmp_path: Path) -> None:
    # ARRANGE
    ruffJson = (
        '[{"code": "E501", "location": {"row": 10}},'
        ' {"code": "E501", "location": {"row": 12}},'
        ' {"code": "B007", "location": {"row": 20}}]'
    )

    # ACT
    counts = surveyor.parse_ruff_counts(ruffJson)

    # ASSERT
    assert counts == {"E501": 2, "B007": 1}


def test_an_empty_ruff_report_counts_nothing() -> None:
    # ACT
    counts = surveyor.parse_ruff_counts("[]")

    # ASSERT
    assert counts == {}


def test_an_unreadable_ruff_report_raises_a_survey_error() -> None:
    # ACT / ASSERT
    with pytest.raises(surveyor.SurveyError):
        surveyor.parse_ruff_counts("not json")


@pytest.mark.parametrize("ruffJson", ["null", '{"code": "E501"}', '["E501"]', '[{"code": 501}]'])
def test_well_formed_json_of_the_wrong_shape_raises_a_survey_error(ruffJson: str) -> None:
    # ACT / ASSERT
    with pytest.raises(surveyor.SurveyError):
        surveyor.parse_ruff_counts(ruffJson)


def test_a_finding_without_a_rule_code_is_counted_as_a_syntax_error() -> None:
    # ACT
    counts = surveyor.parse_ruff_counts('[{"code": null, "location": {"row": 1}}]')

    # ASSERT
    assert counts == {"syntax-error": 1}


# COVERAGE


def _coverage_report(tmp_path: Path, body: str) -> Path:
    reportPath = tmp_path / "coverage.json"
    reportPath.write_text(body)
    return reportPath


def test_coverage_is_read_from_the_surveyed_modules_own_entry(tmp_path: Path) -> None:
    # ARRANGE
    reportPath = _coverage_report(
        tmp_path,
        '{"files": {"soxspipe/recipes/soxs_mbias.py": {"summary": {"percent_covered": 90.43}},'
        ' "soxspipe/commonutils/toolkit.py": {"summary": {"percent_covered": 12.0}}},'
        ' "totals": {"percent_covered": 78.61}}',
    )

    # ACT
    percent = surveyor.read_coverage_percent(reportPath, "soxspipe/recipes/soxs_mbias.py")

    # ASSERT
    assert percent == pytest.approx(90.43)


def test_a_module_missing_from_the_coverage_report_raises_rather_than_reading_zero(tmp_path: Path) -> None:
    # ARRANGE
    reportPath = _coverage_report(tmp_path, '{"files": {}, "totals": {"percent_covered": 78.61}}')

    # ACT / ASSERT
    with pytest.raises(surveyor.SurveyError):
        surveyor.read_coverage_percent(reportPath, "soxspipe/recipes/soxs_mbias.py")


def test_an_unreadable_coverage_report_raises_a_survey_error(tmp_path: Path) -> None:
    # ARRANGE
    reportPath = _coverage_report(tmp_path, "not json")

    # ACT / ASSERT
    with pytest.raises(surveyor.SurveyError):
        surveyor.read_coverage_percent(reportPath, "soxspipe/recipes/soxs_mbias.py")


# REPORT RENDERING


def _survey(**overrides):
    defaults = {
        "path": Path("soxspipe/recipes/soxs_mbias.py"),
        "lineCount": 481,
        "functions": [surveyor.FunctionLength("soxs_mbias.produce_product", 149, 131)],
        "hits": [surveyor.Hit("naive-utcnow", 247), surveyor.Hit("qc-row-builder", 329)],
        "bareExcepts": [],
        "coveragePercent": 90.43,
        "ruffCounts": {"E501": 2},
        "docstringFindings": 3,
        "driftedFunctions": 1,
    }
    defaults.update(overrides)
    return surveyor.Survey(**defaults)


def test_the_report_states_the_line_count_and_coverage() -> None:
    # ACT
    report = surveyor.render_report(_survey())

    # ASSERT
    assert "- 481 lines." in report
    assert "- Coverage 90.4%." in report


def test_the_over_fifty_line_list_stays_in_source_order_when_a_function_exceeds_two_hundred() -> None:
    # ARRANGE
    functions = [
        surveyor.FunctionLength("recipe.first", 10, 60),
        surveyor.FunctionLength("recipe.second", 100, 300),
        surveyor.FunctionLength("recipe.third", 500, 60),
    ]

    # ACT
    report = surveyor.render_report(_survey(functions=functions))

    # ASSERT
    overFifty = next(line for line in report.splitlines() if line.startswith("- Functions over 50 lines:"))
    assert overFifty.index("recipe.first") < overFifty.index("recipe.second") < overFifty.index("recipe.third")


def test_the_report_names_every_oversized_function_with_its_length_and_line() -> None:
    # ACT
    report = surveyor.render_report(_survey())

    # ASSERT
    assert "- Functions over 50 lines: 1 — `soxs_mbias.produce_product` (131, L149)." in report
    assert "- Functions over 200 lines: none." in report


def test_the_report_groups_duplication_hits_by_kind_with_their_lines() -> None:
    # ACT
    report = surveyor.render_report(_survey())

    # ASSERT
    expected = (
        "- Duplication hits: 2 — 1 naive `utcnow` timestamp site (L247),"
        " 1 QC-table `pd.concat` row builder (L329)."
    )
    assert expected in report


def test_the_report_says_none_rather_than_zero_when_there_is_nothing_to_report() -> None:
    # ACT
    empty = _survey(hits=[], functions=[], ruffCounts={}, docstringFindings=0, driftedFunctions=0)
    report = surveyor.render_report(empty)

    # ASSERT
    assert "- Functions over 50 lines: none." in report
    assert "- Duplication hits: none." in report
    assert "- Ruff findings: none." in report
    assert "- Docstring drift: none." in report


def test_the_report_states_the_bare_except_and_savefig_counts() -> None:
    # ACT
    report = surveyor.render_report(_survey(bareExcepts=[100, 200]))

    # ASSERT
    assert "- Bare excepts: 2 (L100, L200). `plt.savefig` calls: none." in report


def test_the_report_marks_coverage_as_not_measured_when_it_was_skipped() -> None:
    # ACT
    report = surveyor.render_report(_survey(coveragePercent=None))

    # ASSERT
    assert "- Coverage not measured." in report


def test_the_report_is_markdown_bullets_only_so_it_pastes_into_a_ticket() -> None:
    # ACT
    report = surveyor.render_report(_survey())

    # ASSERT
    assert all(line.startswith("- ") for line in report.splitlines() if line.strip())


# DOCSTRING DRIFT


def test_docstring_drift_counts_the_findings_and_the_drifted_functions(tmp_path: Path) -> None:
    # ARRANGE
    modulePath = tmp_path / "drifted.py"
    modulePath.write_text(
        '''
def prepare_frame(frame, save=False):
    """*prepare a frame*

    **Key Arguments:**

    - ``frame`` -- the frame to prepare
    """
    return frame
'''
    )

    # ACT
    findings, drifted = surveyor.measure_docstring_drift(modulePath, REPO_ROOT)

    # ASSERT
    assert drifted == 1
    assert findings >= 1


def test_a_module_the_docstring_checker_cannot_parse_raises_a_survey_error(tmp_path: Path) -> None:
    # ARRANGE
    modulePath = tmp_path / "broken.py"
    modulePath.write_text("def broken(:\n")

    # ACT / ASSERT
    with pytest.raises(surveyor.SurveyError):
        surveyor.measure_docstring_drift(modulePath, REPO_ROOT)


def test_an_absent_docstring_checker_raises_a_survey_error(tmp_path: Path) -> None:
    # ARRANGE
    modulePath = tmp_path / "sample.py"
    modulePath.write_text("value = 1\n")

    # ACT / ASSERT
    with pytest.raises(surveyor.SurveyError):
        surveyor.measure_docstring_drift(modulePath, tmp_path)


# PATHS AND COMMANDS


def test_a_module_path_is_reported_relative_to_the_repository_root() -> None:
    # ACT
    modulePath = surveyor._repository_relative(REPO_ROOT / "soxspipe" / "recipes" / "soxs_mbias.py", REPO_ROOT)

    # ASSERT
    assert modulePath == "soxspipe/recipes/soxs_mbias.py"


def test_a_module_outside_the_repository_raises_a_survey_error(tmp_path: Path) -> None:
    # ACT / ASSERT
    with pytest.raises(surveyor.SurveyError):
        surveyor._repository_relative(tmp_path / "stray.py", REPO_ROOT)


def test_a_command_that_cannot_be_run_raises_a_survey_error(tmp_path: Path) -> None:
    # ACT / ASSERT
    with pytest.raises(surveyor.SurveyError):
        surveyor._run(["a-command-that-does-not-exist"], tmp_path)


def test_a_command_failing_with_an_unexpected_status_raises_a_survey_error() -> None:
    # ACT / ASSERT
    # THE REPOSITORY ROOT, NOT A TEMPORARY DIRECTORY, SO A SUBPROCESS MEASURED BY
    # pytest-cov FINDS THE BRANCH COVERAGE SETTING ITS PARENT RUN IS USING
    with pytest.raises(surveyor.SurveyError):
        surveyor._run([sys.executable, "-c", "import sys; sys.exit(3)"], REPO_ROOT)


def test_a_command_failing_with_an_allowed_status_returns_its_output() -> None:
    # ACT
    output = surveyor._run([sys.executable, "-c", "print('done'); raise SystemExit(1)"], REPO_ROOT, (0, 1))

    # ASSERT
    assert output.strip() == "done"


def test_the_ruff_count_runs_ruff_over_the_module(monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    calls = []

    def _fake_run(command, workingDirectory, allowedStatuses=(0,)):
        calls.append((command, allowedStatuses))
        return '[{"code": "E501", "location": {"row": 1}}]'

    monkeypatch.setattr(surveyor, "_run", _fake_run)

    # ACT
    counts = surveyor.count_ruff_findings("soxspipe/recipes/soxs_mbias.py", REPO_ROOT)

    # ASSERT
    command, allowedStatuses = calls[0]
    assert counts == {"E501": 1}
    assert command[:2] == ["ruff", "check"]
    assert command[-2:] == ["--", "soxspipe/recipes/soxs_mbias.py"]
    # RUFF EXITS 1 WHENEVER IT FINDS SOMETHING, WHICH IS THE NORMAL CASE FOR A SURVEY
    assert allowedStatuses == (0, 1)


def test_coverage_measures_the_package_not_the_single_module(monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    calls = []

    def _fake_run(command, workingDirectory, allowedStatuses=(0,)):
        calls.append((command, allowedStatuses))
        reportPath = Path(command[-3].split(":", 1)[1])
        reportPath.write_text('{"files": {"soxspipe/recipes/soxs_mbias.py": {"summary": {"percent_covered": 90.1}}}}')
        return ""

    monkeypatch.setattr(surveyor, "_run", _fake_run)

    # ACT
    percent = surveyor.measure_coverage("soxspipe/recipes/soxs_mbias.py", REPO_ROOT)

    # ASSERT
    command, allowedStatuses = calls[0]
    assert percent == pytest.approx(90.1)
    assert f"--cov={surveyor.COVERAGE_PACKAGE}" in command
    # A RED SUITE STILL MEASURES COVERAGE, SO PYTEST EXITING 1 MUST NOT ABORT THE SURVEY
    assert allowedStatuses == (0, 1)


# COMMAND LINE


def test_the_command_line_prints_the_report_and_exits_clean(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture
) -> None:
    # ARRANGE
    monkeypatch.chdir(REPO_ROOT)
    monkeypatch.setattr(surveyor, "count_ruff_findings", lambda modulePath, repoRoot: {"E501": 4})

    # ACT
    status = surveyor.main(["soxspipe/recipes/soxs_mbias.py", "--no-coverage"])

    # ASSERT
    report = capsys.readouterr().out
    assert status == surveyor.EXIT_CLEAN
    assert "- Coverage not measured." in report
    assert "- Ruff findings: 4 — E501 ×4." in report


def test_the_command_line_measures_coverage_unless_it_is_skipped(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture
) -> None:
    # ARRANGE
    monkeypatch.chdir(REPO_ROOT)
    monkeypatch.setattr(surveyor, "count_ruff_findings", lambda modulePath, repoRoot: {})
    monkeypatch.setattr(surveyor, "measure_coverage", lambda modulePath, repoRoot: 90.13)

    # ACT
    status = surveyor.main(["soxspipe/recipes/soxs_mbias.py"])

    # ASSERT
    assert status == surveyor.EXIT_CLEAN
    assert "- Coverage 90.1%." in capsys.readouterr().out


def test_the_command_line_exits_two_when_the_survey_cannot_be_produced(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture
) -> None:
    # ARRANGE
    monkeypatch.chdir(REPO_ROOT)

    # ACT
    status = surveyor.main(["soxspipe/recipes/no_such_module.py", "--no-coverage"])

    # ASSERT
    assert status == surveyor.EXIT_TOOL_ERROR
    assert "module survey:" in capsys.readouterr().err
