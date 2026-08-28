"""Package-wide logging.

Sapphire historically reported via ``print`` and by appending to ``Sapphire_Info.txt`` /
``Sapphire_Errors.log`` in the run directory. Those files are kept (they are part of the
user-facing contract) but console output now goes through :mod:`logging`, so callers can
silence or redirect it: ``logging.getLogger("Sapphire").setLevel(logging.WARNING)``.
"""
import logging

_ROOT = "Sapphire"


def get_logger(name: str | None = None) -> logging.Logger:
    """Return a child of the ``Sapphire`` logger. Attaches a plain stderr handler once."""
    root = logging.getLogger(_ROOT)
    if not root.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(name)s: %(message)s"))
        root.addHandler(handler)
        root.setLevel(logging.INFO)
    if name and name != _ROOT:
        return root.getChild(name.removeprefix(_ROOT + "."))
    return root
