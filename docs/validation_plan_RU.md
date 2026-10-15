# OpenChipCool v0.4 — Qu & Mudawar (2002): проверяемая схема валидации

**На дату формирования данного пакета: реальный OpenFOAM CHT НЕ ЗАПУЩЕН, сетки НЕ ПРОВЕРЕНЫ `checkMesh`, численные поля CFD НЕ ПОЛУЧЕНЫ, экспериментальная валидация НЕ ДОСТИГНУТА.** Это важное инженерное ограничение, а не незаметная стадия работы.

## 1. Первичные источники

1. Qu, W. & Mudawar, I. (2002), *Experimental and numerical study of pressure drop and heat transfer in a single-phase micro-channel heat sink*, IJHMT **45** (12), 2549–2565. DOI https://doi.org/10.1016/S0017-9310(01)00337-4 ; [авторский PDF](https://engineering.purdue.edu/mudawar/files/articles-all/2002/2002_1.pdf). Fig. 2–3 — установка/схема ячейки; Table 1 — размеры; Fig. 4 — Δp(Re) между пленумами; Fig. 5 — рост температуры воды; Fig. 7 — температуры *меди* в точках термопар.
2. IAPWS (2011), *Revised Supplementary Release on Properties of Liquid Water at 0.1 MPa*, [official release](https://iapws.org/technical-guidance/release/LiquidWater). Eq. (7), Table 5 — вязкость; Eq. (8), Table 6 — теплопроводность; Table 8 — контрольные значения для верификации программной реализации. В v0.4 используется **табличная интерполяция** этих корреляций в CHT словаре OpenFOAM `transport tabulated` (исполняемая совместимость не проверена). ρ и cp по-прежнему постоянны.
3. Celik, I. B., Ghia, U., Roache, P. J., Freitas, C. J., Coleman, H. & Raad, P. E. (2008), *Procedure for Estimation and Reporting of Uncertainty Due to Discretization in CFD Applications*, J. Fluids Engineering **130**, 078001, DOI https://doi.org/10.1115/1.2960953 . Используется методология GCI, но **не** заявление о достигнутой сходимости.
4. OpenFOAM.com v2312 documentation: [steady conjugate heat transfer](https://doc.openfoam.com/2312/tools/processing/solvers/rtm/heat-transfer/), [surfaceFieldValue](https://doc.openfoam.com/2312/tools/post-processing/function-objects/field/surfaceFieldValue/), [tabulated fluid transport](https://doc.openfoam.com/2312/tools/processing/models/thermophysical/transport/rtm/tabulated/).

## 2. Принципиальные отличия настоящего эксперимента от текущей вычислительной геометрии

| Элемент | Эксперимент / статья | v0.4 core CFD |
|---|---|---|
| 21 микроканал, 231×713 мкм, шаг 467 мкм | Есть, Table 1 | Геометрия подготовлена |
| Длина 44.764 мм, ширина ≈10 мм | Есть, Table 1 и Fig. 2 | Геометрия подготовлена |
| Тепловой поток 100/200 W/cm² | По площади *верхней плоскости* | Нижняя граница меди: равная общая мощность |
| Inlet water 15°C | Эксперимент | Задан |
| Верхняя поликарбонатная крышка, 700 мкм | Table 1, численная модель статьи | **Не смоделирована**: адиабатическая стенка |
| 3 поперечных глубоких прорези по 1.6 мм и воздух | Fig. 2–3 | **Не смоделированы**, локальные TC температуры смещены |
| 4 термопары (x = 5, 16.588, 28.176, 39.764 мм) | Table 1, Fig. 7 | Пробы x подготовлены, поперечное положение и отверстия требуют проверки |
| Входной/выходной глубокий и мелкий пленумы, давление между отводами | Fig. 2; Eq. (18)–(22) | **Не смоделированы**: необходимо отдельно установить гидравлические сопротивления |
| Температурная зависимость вязкости | Использовалась в статье | Словарь `tabulated` по IAPWS **подготовлен**, ещё не скомпилирован/протестирован |
| Итерационная / сеточная сходимость | В статье сравнение собственных CFD сеток | **Нет результатов солвера**; известен лишь расчёт числа ячеек |

**Физический смысл:** совмещение измеренного *полного* Δp из Fig. 4 с расчётом потерь только в микроканалах некорректно. Сравнение Fig. 7 с моделью без глубоких прорезей также не может служить полным подтверждением CHT-модели установки. Нельзя калибровать неизвестные потери пленумов по Fig. 4 и далее именовать те же самые точки *независимой валидацией*.

## 3. Корректное сравнение рис. 4, 5, 7

- **Fig. 4:** `Δp_total = Δp_core + Δp_c1 + Δp_c2 + Δp_e2 + Δp_e1`. В `src/openchipcool/cfd/validation.py::pressure_budget` реализованы Eq. (19)–(22) Qu & Mudawar. Для этой оценки **обязательны независимые** площади/скорости обоих пленумов, определённые из чертежа или отдельного расчёта. `Δp_core` берётся из area-averaged `p` OpenFOAM при x=0 и x=L. Реальное измерение идёт между *глубокими* пленумами.
- **Fig. 5:** сравнить `T_out − T_in` в Kelvin с экспериментальными точками; `T_out` должен быть **массово-взвешенным** по выходному патчу. Проверить `Q_water = mdot Cp (T_out − T_in)` и отличать энергетическую проверку от CHT-валидации.
- **Fig. 7:** сравнить `T_tc,i − T_in` вдоль продольной линии термопар; **не** сравнивать температура поверхности канала с термопарой внутри меди. Нужно отдельно восстановить отверстия и прорези, подтвердить положение центра спая TC в Fig. 2–3. В v0.4 заданы лишь приближённые координаты с z=3.175 мм, место сверить до использования.
- Для Fig. 4 готов `qu_mudawar_2002_pressure_approx.csv`, для Fig. 5 — `qu_mudawar_2002_tout_fig5_visual_estimate.csv`, для Fig. 7 — `qu_mudawar_2002_tc_fig7_visual_estimate.csv`. Это **ручные приближённые считывания с печатных графиков, НЕ данные авторов**; имеют заявленные в CSV условные интервалы ошибки оцифровки. Нужна повторная цифровая оцифровка высокого разрешения и два независимых чтения.

## 4. Контроль качества расчёта и критерии (предложения, НЕ достигнуты)

1. `blockMesh` + `splitMeshRegions -cellZonesOnly -overwrite` + `checkMesh -region water/copper -allTopology -allGeometry`; сохранённые полные логи, без ошибок/отрицательных объёмов.
2. Запуск OpenFOAM.com `chtMultiRegionSimpleFoam`, зафиксировать версию (`foamVersion`), рабочие словари после адаптации, время CPU, трассу residuals, невязку расхода, невязку энергии, карту давлений и температур; сам факт напечатанного `End` ещё не есть сходимость.
3. Три согласованные 3D-сетки, r=2 по каждой координате. Для каждого рабочей точки: `Δp_core`, `T_out`, `T_tc1...T_tc4`, баланс массы и энергии. Если величина монотонна и наблюдается асимптотический режим, оценить p и GCI (Celik et al.). Если немонотонна, **не выдавать GCI без оговорок**.
4. Предлагаемые пороги: mass imbalance <0.1%, energy imbalance <1%, GCI целевых величин <5%, средний абсолютный разброс с точками Fig. 7 <5 K и ошибки Δp_total <10%. Пороги **только внутренние цели проекта**; должны быть пересмотрены с учётом погрешности измерений и чтения графиков.
5. Валидировать **независимо**: некоторые режимы для настройки модельных параметров, другие удержать для hold-out. Если геометрия пленумов недоступна, не сообщать full-facility-validation, а явно обозначить scope = core CHT with component-level pressure budget.
6. Переходить к 36-мм концепции **только после** указанных проверок; отдельно оценить влияние четырех неоднородных тепловых карт и давления на распределение расхода по каналам. В перенесённой геометрии необходимо заново рассчитать Δp(flow), h_eq(flow), R_th(flow) и *их неопределённость*. Эксперимент Qu–Mudawar **не является прямой экспериментальной валидацией нового устройства**.

## 5. Воспроизводимый запуск на установленном OpenFOAM.com

```bash
# Войти в окружение своей конкретной установленной версии OpenFOAM.com
# например: source /usr/lib/openfoam/openfoam2312/etc/bashrc
python -m pip install -e '.[dev]'
pytest -q
python scripts/build_validation_matrix.py --output cfd_cases/validation_matrix
# Требуются команды: blockMesh, splitMeshRegions, checkMesh, chtMultiRegionSimpleFoam
python scripts/run_openfoam_matrix.py cfd_cases/validation_matrix/qm_q100_Re890_coarse --execute
# Только после успеха, корректной интерпретации pressure taps и независимой оценки пленумов:
python scripts/extract_openfoam_result.py cfd_cases/validation_matrix/qm_q100_Re890_coarse \
    --plenum-correction-pa 1234 \
    --plenum-justification path/to_independent_plenum_calculation.md \
    --out cfd_exports/q100_Re890_coarse.csv
# 1234 Pa ВЫШЕ — только формат ввода, не численный результат! 
python scripts/validate_all_figures.py --cfd-export cfd_exports/combined_validated_exports.csv
```

**Примечания к исполнению:** инструмент `extract_openfoam_result.py` требует реальные `postProcessing`-файлы и полный лог; ненулевые показатели контроля проверяются. Оценка плёнумов невозможна без независимо обоснованных скоростей, поэтому нулевое/фиктивное значение **не разрешает** заявлять полную валидацию. Скрипт сравнений тоже не превращает диаграммы в сертификат экспериментального подтверждения. Набор `--execute` имеет строгие nonzero exit при отсутствии инструмента или при ошибке.

## 6. Фактический статус

- Python static geometry, zone face bookkeeping, water properties, matrix generation — **unit tests PASS**.
- CFD mesh and regional checkMesh — **NOT RUN**.
- Real CHT output and pressure/temperature fields — **NOT RUN**.
- GCI, сравнение с Fig 4/5/7 — **NOT RUN**.
- Коллекторы и глубокие слоты — **NOT MODELED**.
- GitHub — **NOT PUBLISHED**, потому что текущий инструмент подключённого GitHub не позволяет создать новый репозиторий.
