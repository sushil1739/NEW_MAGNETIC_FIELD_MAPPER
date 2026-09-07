from pathlib import Path

from tms_mapper.edep import read_mapper_header, validate_mapper_file


def test_validate_edep_order(tmp_path: Path):
    p = tmp_path / "Mapper.txt"
    p.write_text(
        "# shape 2 2 2\n"
        "# comment\n"
        "0 0 0 10 20 30\n"
        "0 0 0 1 0 0 1\n"
        "0 0 30 1 0 0 1\n"
        "0 20 0 1 0 0 1\n"
        "0 20 30 1 0 0 1\n"
        "10 0 0 1 0 0 1\n"
        "10 0 30 1 0 0 1\n"
        "10 20 0 1 0 0 1\n"
        "10 20 30 1 0 0 1\n"
    )
    h = read_mapper_header(p)
    assert h.spacing == (10.0, 20.0, 30.0)
    result = validate_mapper_file(p)
    assert result["shape"] == (2, 2, 2)
    assert result["rows"] == 8
