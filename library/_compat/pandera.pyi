from types import ModuleType
from typing import Any

class _PanderaModule(ModuleType):
    Column: Any
    Check: Any
    DataFrameSchema: Any
    DataFrameModel: Any


def load_pandera_pandas() -> _PanderaModule: ...

pa: _PanderaModule
Column: Any
Check: Any
DataFrameSchema: Any
DataFrameModel: Any
