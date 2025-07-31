"""Placeholder utilities for configuration handling.

This project no longer reads configuration files. ``load_config`` is kept
for compatibility but always returns an empty dictionary.
"""


def load_config(_=None):
    """Return an empty configuration dict.

    Parameters are ignored because loading from files is not supported.
    """
    return {}
