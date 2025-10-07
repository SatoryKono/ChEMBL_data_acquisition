# Настройка окружения

## Инструменты

- Python 3.11 или 3.12.
- `pip` + `requirements-lock.txt` для воспроизводимых установок.
- `pre-commit` (форматирование `black`, линт `ruff`).
- Опционально `.env` с переменными `CHEMBL_DA_*`.

## Установка

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements-lock.txt
pre-commit install
```

Для работы только с документацией можно установить `requirements-dev.txt`.

## Основные каталоги

| Путь | Описание |
|------|----------|
| `library/` | Код пакета. |
| `scripts/` | CLI-обёртки над пайплайнами. |
| `config/` | YAML и словари; локальные правки — в `config/config.local.yaml`. |
| `tests/` | Тесты (`unit`, `integration`, `e2e`). |
| `data/input` | Тестовые CSV. |
| `reports/` | Итоговые отчёты pytest. |

## Форматирование

```bash
make lint
make typecheck
make format
```

`pre-commit` запускает те же проверки перед коммитом.

## Переменные окружения

| Переменная | Назначение |
|------------|------------|
| `CHEMBL_DA_BASE_PATH` | Базовый путь для плейсхолдера в `config/config.yaml`. |
| `CHEMBL_DA_*` | Алиасы конфигурации (см. `library/config/models.py`). |
| `HTTP_PROXY` / `HTTPS_PROXY` | Прокси в корпоративных сетях. |

Чувствительные значения храните в `.env` (подхватывается `python-dotenv`).
