from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from auto_nest import parse_args

def test_parse_basic():
    args = parse_args(['-s','100x100','part.dxf'])
    assert args.sheet == '100x100'
    assert args.files == ['part.dxf']
