
    # Выходные данные пайплайна целей: постобработка isoform

    Дополнение к [описанию выходов](./ru/OUTPUT.md), описывающее этап
    `isoform.*`, выполняемый поверх CSV `output.target_YYYYMMDD.csv`. Эквиваличный
    материал на английском приведён в
    [`OUTPUT_TARGETS_EN.md`](./OUTPUT_TARGETS_EN.md).

    ## Постобработка `isoform.*`

    ### Этапы
    1. **Загрузка входного CSV.** Файл читается с перебором кодировок
       `utf-8`, `utf-8-sig`, `cp1252`. Из таблицы выбираются только столбцы
       `isoform_synonyms`, `isoform_names`, `isoform_ids`, `uniprot_id_primary`,
       `target_chembl_id`.
    2. **Нормализация строк.** `isoform_synonyms` и `isoform_names`
       переводятся в нижний регистр и разбиваются по символу `|` с обрезкой
       пробелов и удалением пустых элементов. `isoform_ids` разбиваются по `|`
       без изменения регистра.
    3. **Выравнивание записей.** Списки имён, идентификаторов и синонимов
       выравниваются по индексу, недостающие элементы заполняются `null`.
    4. **Токенизация синонимов.** Каждый синоним разбивается по `":"`, пустые
       части удаляются. Для каждого токена строятся варианты
       `[token, token без "pde", token без "pld"]` с устранением дубликатов.
    5. **Формирование таблиц.** Создаются две таблицы: `names` — из
       нормализованных `isoform_names`, и `synonyms` — из токенов. После
       переименования столбца `tokens` в `name` обе таблицы объединяются.
    6. **Удаление дублей.** Последовательно выполняются три шага:
       `distinct(id, name, target_chembl_id, uniprot_id_primary)`, стабильная
       сортировка по `(uniprot_id_primary, id)` и `distinct(id, target_chembl_id, name)`,
       затем финальный `distinct(id, name)`.
    7. **Экспорт.** Результат записывается в файл
       `isoform.output.<basename>.csv` рядом с исходным CSV (или по указанному
       пути) в кодировке UTF-8 без BOM.

    ### Инварианты
    - Столбцы результата: ровно `id`, `uniprot_id_primary`, `target_chembl_id`,
      `name` в указанном порядке.
    - `name` всегда в нижнем регистре; `id` сохраняет регистр входного
      `isoform_ids` (заполняется `null`, если идентификатор отсутствовал).
    - Дубликаты удаляются строго в порядке, задающем стабильный и
      воспроизводимый вывод.
    - Поля `"", "n/a", "none"` отфильтровываются перед объединением таблиц.

    ### Пример (вход → выход)

    Входной фрагмент `output.target_20250228.csv`:

    ```csv
    target_chembl_id,uniprot_id_primary,isoform_ids,isoform_names,isoform_synonyms
    CHEMBL3587,P12345,"ENSP0001|ENSP0002","Alpha|Beta","PDE3A:Alpha|Alpha variant|Beta"
    CHEMBL1250,P67890,"ENSP0003","Gamma","Gamma isoform 1|PLD1A"
    CHEMBL3135,Q11111,"","",""
    CHEMBL2205,Q22222,"ENSP0004","None","None"
    ```

    Результат `isoform.output.target_20250228.csv`:

    ```csv
    id,uniprot_id_primary,target_chembl_id,name
    ,P12345,CHEMBL3587,beta
    ENSP0001,P12345,CHEMBL3587,alpha
    ENSP0001,P12345,CHEMBL3587,pde3a
    ENSP0001,P12345,CHEMBL3587,3a
    ENSP0002,P12345,CHEMBL3587,beta
    ENSP0002,P12345,CHEMBL3587,alpha variant
    ,P67890,CHEMBL1250,pld1a
    ,P67890,CHEMBL1250,1a
    ENSP0003,P67890,CHEMBL1250,gamma
    ENSP0003,P67890,CHEMBL1250,gamma isoform 1
    ```

    ### Детерминизм
    - Сортировка выполняется стабильным `mergesort`, что обеспечивает одинаковый
      порядок строк при повторных запусках.
    - Операции не зависят от локали, даты или случайных чисел; результат
      воспроизводим при неизменном входе.
    - Кодировка вывода — UTF-8 без BOM, перевод строки `
`.

    ### Совместимость
    - Функция `library.postprocessing.target.process_targets` доступна начиная с
      версии `chembl-data-acquisition` 0.1.3.
    - Требуемая версия Python — 3.11+, библиотека `pandas` — 2.3.x.
    - Ожидается, что входной файл соответствует схеме таргетов, описанной в
      [`docs/ru/OUTPUT.md`](./ru/OUTPUT.md).

    ### Пример запуска

    ```bash
    python - <<'PY'
    from library.postprocessing import target
    target.process_targets("data/output/output.target_20250228.csv")
    PY
    ```

    В результате появится файл `data/output/isoform.output.target_20250228.csv`
    с постобработанными именами изоформ.
