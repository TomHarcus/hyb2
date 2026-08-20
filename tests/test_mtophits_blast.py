import subprocess

from hyb2.chimera.mtophits_blast import deduplicate_by_second_fragement_start

GROUPED = (
    "READ1\tsubjA\tc3\tc4\tc5\tc6\tc7\tc8\tc9\tc10\t1e-10\n"
    "READ1\tsubjC\tc3\tc4\tc5\tc6\tc7\tc8\tc9\tc10\t1e-10\n"
    "READ1\tsubjE\tc3\tc4\tc5\tc6\tc7\tc8\tc9\tc10\t1e-02\n"
    "READ2\tsubjB\tc3\tc4\tc5\tc6\tc7\tc8\tc9\tc10\t1e-05\n"
    "READ2\tsubjD\tc3\tc4\tc5\tc6\tc7\tc8\tc9\tc10\t1e-20\n"
)

INTERLEAVED = (
    "READ1\tsubjA\tcol3\tcol4\tcol5\tcol6\tcol7\tcol8\tcol9\tcol10\t1e-10\n"
    "READ2\tsubjB\tcol3\tcol4\tcol5\tcol6\tcol7\tcol8\tcol9\tcol10\t1e-05\n"
    "READ1\tsubjC\tcol3\tcol4\tcol5\tcol6\tcol7\tcol8\tcol9\tcol10\t1e-10\n"
    "READ2\tsubjD\tcol3\tcol4\tcol5\tcol6\tcol7\tcol8\tcol9\tcol10\t1e-20\n"
)


def _awk_reference(text: str) -> str:
    result = subprocess.run(
        ["awk", "-f", "legacy_bin/mtophits_blast"],
        input=text,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def test_matches_legacy_awk_output_grouped():
    assert deduplicate_by_second_fragement_start(GROUPED) == _awk_reference(GROUPED)


def test_keeps_ties_and_drops_worse_evalue_within_group():
    result = deduplicate_by_second_fragement_start(GROUPED)
    assert result.count("READ1") == 2
    assert "subjE" not in result


def test_matches_legacy_awk_shared_state_quirk_on_interleaved_ids():
    result = deduplicate_by_second_fragement_start(INTERLEAVED)
    assert result == _awk_reference(INTERLEAVED)
    assert result.count("\n") == 2
