from pathlib import Path
import sys
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from auto_nest import parse_args

def test_parse_basic():
    args = parse_args(['-s','100x100','part.dxf'])
    assert args.sheet == '100x100'
    assert args.files == ['part.dxf']

def test_parse_flags():
    argv = ['-s','50x50','--pop','10','--gen','5','--polish','1',
            '--fill-gaps','file1.dxf']
    args = parse_args(argv)
    assert args.sheet == '50x50'
    assert args.pop == 10
    assert args.gen == 5
    assert args.polish == 1
    assert args.fill_gaps is True

from config_util import load_config


def test_invalid_arg():
    with pytest.raises(SystemExit):
        parse_args(['--unknown'])


def test_load_config_missing(tmp_path):
    cfg = load_config(tmp_path / 'nope.yaml')
    assert cfg == {}


def test_load_config_ignored(tmp_path):
    cfg_path = tmp_path / 'cfg.yaml'
    cfg_path.write_text('sheet: 50x50\niterations: 2')
    cfg = load_config(cfg_path)
    assert cfg == {}


