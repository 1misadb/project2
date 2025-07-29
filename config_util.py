import yaml
from pathlib import Path

DEFAULT_CONFIG = Path('config.yaml')

def load_config(path=None):
    path = Path(path or DEFAULT_CONFIG)
    if path.exists():
        with open(path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f) or {}
    return {}
