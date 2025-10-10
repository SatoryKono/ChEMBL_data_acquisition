"""Unified entry points for all first-party console scripts."""

from .activity import ActivityPipelineCLI, build_parser
from .activity import main as activity_main
from .commands import (
    check_determinism_main,
    chunk_io_main,
    csv_utils_main,
    get_activities_main,
    get_activity_data_main,
    get_assay_data_main,
    get_data_main,
    get_document_data_main,
    get_document_type_main,
    get_input_initialisation_main,
    get_target_data_main,
    get_testitem_data_main,
    make_md_summary_main,
    mapper_main,
    table_quality_main,
)

__all__ = [
    "ActivityPipelineCLI",
    "activity_main",
    "build_parser",
    "check_determinism_main",
    "chunk_io_main",
    "csv_utils_main",
    "get_activities_main",
    "get_activity_data_main",
    "get_assay_data_main",
    "get_data_main",
    "get_document_data_main",
    "get_document_type_main",
    "get_input_initialisation_main",
    "get_target_data_main",
    "get_testitem_data_main",
    "make_md_summary_main",
    "mapper_main",
    "table_quality_main",
]
