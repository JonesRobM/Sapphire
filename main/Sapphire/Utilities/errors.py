"""One place for Sapphire's "log it and carry on" behaviour.

Historically every module caught ``Exception`` and appended to ``Sapphire_Errors.log``; nothing
ever raised, so breakage (e.g. the numpy-2 ``trapz`` removal) went unnoticed for years.
:func:`report` keeps the log-file contract but also emits a ``logging`` warning and, when
strict mode is on (``Process(strict=True)`` or :func:`set_strict`), re-raises.
"""
from Sapphire.Utilities.log import get_logger

log = get_logger("Sapphire.errors")
_STRICT = False


def set_strict(flag: bool) -> None:
    global _STRICT
    _STRICT = bool(flag)


def is_strict() -> bool:
    return _STRICT


def report(exc, message, base_dir=None, strict=None):
    """Record *exc* with context *message*; re-raise if strict.

    ``message`` may contain one ``%s`` for the exception text. ``base_dir`` (a run directory)
    selects the ``Sapphire_Errors.log`` to append to; ``None`` logs to the console only.
    """
    text = message % exc if "%s" in message else f"{message}: {exc}"
    if base_dir:
        try:
            with open(str(base_dir) + "Sapphire_Errors.log", "a") as f:
                f.write("\n" + text.strip() + "\n")
        except OSError:  # pragma: no cover - unwritable base_dir
            pass
    log.warning(text.strip().replace("\n", " "))
    if _STRICT if strict is None else strict:
        raise exc
