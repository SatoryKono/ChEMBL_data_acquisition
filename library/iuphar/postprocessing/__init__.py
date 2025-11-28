"""Post-processing routines for IUPHAR data."""

from .iuphar import IUPHARPostProcessingError, process_iuphar_targets

__all__ = ["IUPHARPostProcessingError", "process_iuphar_targets"]
