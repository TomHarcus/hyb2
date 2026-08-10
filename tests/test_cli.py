from hyb2.cli import build_parser, main


def test_no_args_prints_help(capsys):
    assert main([]) == 0
    assert "usage" in capsys.readouterr().out.lower()


def test_parses_legacy_flags():
    args = build_parser().parse_args(
        ["-i", "test.sam", "-d", "ref.fasta", "-o", "test", "-a", "Zika_virusRNA", "-x", "1001", "-l", "500"]
    )
    assert args.in_file == "test.sam"
    assert args.db_1 == "ref.fasta"
    assert args.gene_1 == "Zika_virusRNA"
    assert args.x_coord == 1001
    assert args.length == 500
    assert args.fold == "cplfold"  # default until couplefold is integrated (Phase 3)
