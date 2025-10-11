"""I/O helpers for CSV-based data pipelines.

The package consolidates CSV reader and writer utilities together with
high-level helpers such as :func:`save_standard_outputs` which persist
canonical pipeline artefacts deterministically.
"""

from importlib import import_module, reload

_metadata = import_module(f"{__name__}.metadata")
_paths = import_module(f"{__name__}.paths")
_writers = import_module(f"{__name__}.writers")
_output_writer = import_module(f"{__name__}.output_writer")
_readers = reload(import_module(f"{__name__}.readers"))

pa = _readers.pa
locale = _readers.locale
read_csv = _readers.read_csv
read_ids = _readers.read_ids
CsvReadError = _readers.CsvReadError
write_csv = _writers.write_csv
save_standard_outputs = _output_writer.save_standard_outputs
StandardOutputArtifacts = _output_writer.StandardOutputArtifacts
default_output_path = _paths.default_output_path
derive_output_labels = _paths.derive_output_labels
write_meta_yaml = _metadata.write_meta_yaml

__all__ = [
    "CsvReadError",
    "derive_output_labels",
    "default_output_path",
    "pa",
    "locale",
    "read_csv",
    "read_ids",
    "write_csv",
    "save_standard_outputs",
    "StandardOutputArtifacts",
    "write_meta_yaml",
]
