"""Factories and helpers for client construction."""

from .client_factory import ClientFactory
from .helpers import build_normalizer, build_transport
from .semantic_scholar_client_factory_impl import SemanticScholarClientFactoryImpl
from .openalex_client_factory_impl import OpenAlexClientFactoryImpl

__all__ = [
    "ClientFactory",
    "OpenAlexClientFactoryImpl",
    "build_normalizer",
    "build_transport",
    "SemanticScholarClientFactoryImpl",
]
