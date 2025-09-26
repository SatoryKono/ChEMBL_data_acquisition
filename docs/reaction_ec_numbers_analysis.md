# Анализ заполнения `reaction_ec_numbers`

## 1. Источник и обработка в `scripts/get_target_data.py`
- Колонка `reaction_ec_numbers` приходит в объединённый датафрейм из выгрузки UniProt (`fetch_uniprot`) и сохраняется в таблице `combined_df` вместе с колонкой `ec_numbers` из UniProt.【F:scripts/get_target_data.py†L880-L907】
- При подготовке данных для обогащения IUPHAR значения `reaction_ec_numbers` используются только для расчёта агрегированной колонки `ec_number`, после чего исходные столбцы `ec_numbers` и `reaction_ec_numbers` целиком удаляются вызовом `DataFrame.drop(...)`. Это означает, что дальнейшие шаги конвейера больше не видят `reaction_ec_numbers`.【F:scripts/get_target_data.py†L904-L913】

## 2. Постобработка и финализация таргетов
- На этапе финализации `target_postprocessing.align_target_columns` ожидает наличие колонки `reaction_ec_numbers`. Если столбца нет, хелпер `_series_or_default` создаёт заполнение из дефолтного значения `"-"`, что и приводит к пустым данным в готовом `output/target.csv`.【F:library/target_postprocessing.py†L160-L170】【F:library/target_postprocessing.py†L268-L287】

## 3. Формирование `reaction_ec_numbers` в `output/target_uniprot.csv`
- Выгрузка UniProt строится функцией `collect_info`, которая извлекает каталитические реакции и EC-номера (через `extract_activity`) и добавляет их в результирующий словарь перед записью в CSV.【F:library/uniprot_library.py†L1129-L1185】
- Список колонок для UniProt-выхода включает `reaction_ec_numbers`, поэтому значение без изменений оказывается в `output/target_uniprot.csv`.【F:library/uniprot_library.py†L88-L140】【F:library/uniprot_library.py†L1239-L1265】

## 4. Причина расхождения
- В конвейере `all` колонка удаляется из объединённых данных перед постобработкой, поэтому финальный CSV получает дефолтные `"-"` вместо реальных значений.
- В чистой UniProt-выгрузке колонка не удаляется и записывается как есть.

## 5. Рекомендации по исправлению
1. **Не удалять колонку при сборке комбинированных данных.** Убрать `reaction_ec_numbers` из списка на удаление или сохранить его под отдельным именем, чтобы финализация увидела исходные данные.
2. **Опционально:** если колонка `ec_number` нужна одновременно с исходным списком реакций, можно оставить оба столбца или пересчитать `ec_number` на более позднем шаге (например, в `target_postprocessing`).
3. После изменений убедиться, что `TargetsSchema` допускает наличие обеих колонок (она уже объявлена как необязательная) и что в `target_postprocessing.finalise_targets` данные больше не заменяются на дефолтные значения.

Внесение этих правок позволит `reaction_ec_numbers` из UniProt пройти весь конвейер и сохраниться в итоговом `output/target.csv`.
