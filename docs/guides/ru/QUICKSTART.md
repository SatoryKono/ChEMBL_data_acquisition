# Быстрый старт

Этот чек-лист помогает перейти от чистого клона репозитория к проверенному запуску пайплайна. Каждый шаг повторяет настройки разработки, чтобы аналитики и инженеры работали в единой среде.

## 1. Подготовьте окружение

1. **Совместите версию Python.** Установите Python `3.11.12` (см. `.python-version`) и обновите `pip`, `setuptools`, `wheel`.

   ```bash
   python -m pip install --upgrade pip setuptools wheel
   ```

2. **Создайте и активируйте виртуальное окружение.** Изолированная среда не загрязняет системный интерпретатор.

   ```bash
   python -m venv .venv
   source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
   ```

3. **Установите зависимости из lock-файла и dev-набор.** Lock-файл синхронизирует локальное окружение с CI, а editable-установка открывает консольные скрипты из документации.

   ```bash
   pip install -r requirements-lock.txt
   pip install -e .[dev]
   ```

4. **Включите автоматизацию репозитория.** Pre-commit-хуки гарантируют запуск форматирования, линтинга, mypy и тестов перед коммитом.

   ```bash
   pre-commit install
   ```

## 2. Выполните минимальный пайплайн

1. **Создайте рабочий каталог для выгрузки.**

   ```bash
   mkdir -p tmp/quickstart
   ```

2. **Запустите пайплайн таргетов на встроенном smoke-входе.** Команда выгружает несколько таргетов ChEMBL, сохраняет «сырой» и очищенный файлы и демонстрирует поэтапное логирование.

   ```bash
   get-target-data chembl \
     --input tests/data/chembl_targets_min.csv \
     --column target_chembl_id \
     --limit 5 \
     --raw-out tmp/quickstart/targets.raw.csv \
     --final-out tmp/quickstart/targets.final.csv
   ```

   *Добавьте `--log-level DEBUG` для подробной диагностики или `--print-config`, чтобы посмотреть итоговую конфигурацию до запуска.*

## 3. Проверьте установку

Выполните проверки качества, повторяющие CI. Предполагается, что виртуальное окружение всё ещё активно.

```bash
make lint
make test
make smoke
```

`make lint` объединяет Ruff, Black (в режиме проверки) и MyPy. `make test` запускает полный набор PyTest, а `make smoke` повторяет сетевые smoke-тесты с локальными фикстурами (пайплайн test-item по умолчанию пропускается).
