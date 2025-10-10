"""I/O helpers for CSV-based data pipelines."""

from importlib import import_module, reload

_metadata = import_module(f"{__name__}.metadata")
_paths = import_module(f"{__name__}.paths")
_writers = import_module(f"{__name__}.writers")
_readers = reload(import_module(f"{__name__}.readers"))

pa = _readers.pa
locale = _readers.locale
read_csv = _readers.read_csv
read_ids = _readers.read_ids
CsvReadError = _readers.CsvReadError
write_csv = _writers.write_csv
default_output_path = _paths.default_output_path
write_meta_yaml = _metadata.write_meta_yaml

__all__ = [
    "CsvReadError",
    "default_output_path",
    "pa",
    "locale",
    "read_csv",
    "read_ids",
    "write_csv",
    "write_meta_yaml",
]
