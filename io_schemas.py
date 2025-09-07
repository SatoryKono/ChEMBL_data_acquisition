"""Pydantic data models for document classification I/O."""

from __future__ import annotations

from typing import List, Optional

from pydantic import BaseModel, Field  # type: ignore[import-not-found]


class DocumentRecord(BaseModel):
    """Normalized input record for classification."""

    document_id: Optional[str] = None
    title: str
    abstract: str
    doi: str
    pubmed_publicationtype: List[str] = Field(default_factory=list)
    scholar_publicationtypes: List[str] = Field(default_factory=list)
    openalex_publicationtypes: List[str] = Field(default_factory=list)
    crossref_type: Optional[str] = None
    openalex_typecrossref: Optional[str] = None
    pubmed_mesh_descriptors: List[str] = Field(default_factory=list)
    pubmed_mesh_qualifiers: List[str] = Field(default_factory=list)
    openalex_meshdescriptors: List[str] = Field(default_factory=list)
    openalex_meshqualifiers: List[str] = Field(default_factory=list)
    pubmed_chemicallist: List[str] = Field(default_factory=list)
    optional_experiment_kind: Optional[str] = None


class ClassificationResult(BaseModel):
    """Result of classification with detailed explanation."""

    document_id: Optional[str] = None
    final_class: str
    S_in_vivo: int
    S_in_vitro: int
    pt_comparative: bool
    pt_binding: bool
    vivo_hits: List[str]
    vitro_hits: List[str]
    conflicts: List[str]
    confidence: int
    source_priority_used: str
    xrf_is_generic: bool
    xrf_is_specific: bool
    explain_short: str
