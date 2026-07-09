"""NValchemi GPU-accelerated backend for MLIP batched operations."""

from src.utils.mlips.nvalchemi.nvalchemi_utils import (
    NVALCHEMI_AVAILABLE,
    check_nvalchemi_available,
    atoms_to_atomic_data,
    atomic_data_to_atoms,
    extract_batch_results,
)

__all__ = [
    "NVALCHEMI_AVAILABLE",
    "check_nvalchemi_available",
    "atoms_to_atomic_data",
    "atomic_data_to_atoms",
    "extract_batch_results",
    "M3GNetWrapper",
    "CHGNetWrapper",
]


# Lazy imports to avoid requiring nvalchemi/matgl at import time
def __getattr__(name: str):  # noqa: N807
    if name in ("M3GNetWrapper", "CHGNetWrapper"):
        from src.utils.mlips.nvalchemi.matgl_wrappers import (
            M3GNetWrapper,
            CHGNetWrapper,
        )

        globals()["M3GNetWrapper"] = M3GNetWrapper
        globals()["CHGNetWrapper"] = CHGNetWrapper
        return globals()[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
