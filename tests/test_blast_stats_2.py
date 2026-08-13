import subprocess

from hyb2.legacy_leaves.blast_stats_2 import blast_stats_2

MULTI_BIOTYPE = (
    "READ1_5\tHuman|X|GENE_A_typeA\t.\t.\t.\t.\t1\t10\t1\t10\t1e-10\t.\n"
    "READ2_3\tHuman|X|GENE_B_typeB\t.\t.\t.\t.\t1\t10\t1\t10\t1e-10\t.\n"
    "READ3_2\tHuman|X|GENE_C_typeA\t.\t.\t.\t.\t1\t10\t1\t10\t1e-10\t.\n"
    "READ4_7\tHuman|X|GENE_D_typeC\t.\t.\t.\t.\t1\t10\t1\t10\t1e-10\t.\n"
    "READ5_1\tHuman|X|GENE_E_typeB\t.\t.\t.\t.\t1\t10\t1\t10\t1e-10\t.\n"
    "READ6_9\tHuman|X|GENE_F_typeC\t.\t.\t.\t.\t1\t10\t1\t10\t1e-10\t.\n"
    "READ7_4\tHuman|X|GENE_G_typeD\t.\t.\t.\t.\t1\t10\t1\t10\t1e-10\t.\n"
)


def _awk_reference(text: str) -> str:
    result = subprocess.run(
        ["awk", "-f", "bin/blast_stats_2"],
        input=text,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def test_matches_legacy_awk_output_on_real_fixture():
    with open("fixtures/legacy_run/direct_viewpoint/test_Zika_virusRNA.blast") as f:
        blast_text = f.read()

    assert blast_stats_2(blast_text) == _awk_reference(blast_text)


def test_matches_legacy_awk_output_multi_biotype_with_ties():
    assert blast_stats_2(MULTI_BIOTYPE) == _awk_reference(MULTI_BIOTYPE)


def test_sorted_descending_by_collapsed_count():
    result = blast_stats_2(MULTI_BIOTYPE)
    rows = result.splitlines()[1:]
    counts = [int(row.split("\t")[1]) for row in rows]
    assert counts == sorted(counts, reverse=True)


def test_header_printed_on_empty_input():
    assert blast_stats_2("") == "biotype\tcollapsed\tall\n"
