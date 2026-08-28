"""Every module in the package must import in a fresh interpreter (optional extras excepted)."""
import importlib
import pkgutil

import pytest

import Sapphire

OPTIONAL = {"Sapphire.Light.ClassSpec": "pyGDM2"}

MODULES = sorted(m.name for m in pkgutil.walk_packages(Sapphire.__path__, "Sapphire."))


@pytest.mark.parametrize("name", MODULES)
def test_module_imports(name):
    if name in OPTIONAL:
        pytest.importorskip(OPTIONAL[name])
    importlib.import_module(name)


def test_version():
    assert Sapphire.__version__
