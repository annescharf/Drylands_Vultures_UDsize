# Shared agent brief: Drylands vulture space-use paper

**Read `AIMS.md` before this brief.** If a task here seems to drift away from the questions in `AIMS.md`, stop and open an issue instead of continuing.

This brief is for **two agents working in the same shared folder and GitHub repo**:
- **[A] Anne's agent**: all processing and preparatory analysis, up to and including the final tables for modelling (this includes the residency classification, exposure extraction and their diagnostics).
- **[M] Martina's agent**: the final statistical models that produce the results reported in the paper, and writing support.

Most of the work before the final models is [A]'s. [M] mainly reviews at Checkpoint 1 and takes over once the modelling tables are delivered.

Every task below is tagged with its owner. Read the whole brief, including the other agent's tasks, because your inputs and outputs depend on each other.

---

## 0. Context

**Project.** Global analysis of how vulture space use responds to vegetation productivity (NDVI) across aridity gradients. Monthly utilization distributions (UDs, dynamic Brownian bridge movement model, dBBMM) for GPS-tracked vultures, 2003–2025, ~13 species, ~2,400 individuals, five realms. Analysis in R (Ubuntu).

**People and roles.**
- **Dennis Kim**: lead author.
- **Anne and Martina**: co-senior authors.
- **Anne** owns all processing and preparatory analysis: raw trajectories and NDVI rasters → annotated UDs → residency classification → exposure tables → final modelling tables.
- **Martina** owns the final statistical models and the writing.
- **Threshold and design decisions** (e.g. overlap threshold, `k_settle`, `D_region`, drought cut-offs, destination-choice design) are **decided jointly by Anne and Martina**. [A] proposes them with evidence (plots, sensitivity tables), Anne and Martina agree, then [A] implements. [A] never fixes such a threshold on its own; log the agreed value and date in DECISIONS.md.

**Existing work.** Martina will put her June 2026 model scripts (`drylands_vulture_udsize_change`) in the repo under `legacy/models/`. Anne will put her existing processing scripts under `legacy/processing/`. The humans do this, not the agents. Agents treat `legacy/` as read-only context and never edit it. If something expected is missing from `legacy/`, open an issue for Anne or Martina rather than looking for it elsewhere.

**Paper framing (decided).** *Adjust, relocate, or persist.* Under changing productivity, a vulture can:
- **Adjust**: change UD size while staying in the same area.
- **Relocate**: leave and settle elsewhere.
- **Persist**: stay without detectable adjustment.

How the strategies map onto the questions (see `AIMS.md`):
- **Adjust** → Q1 (Task 3).
- **Relocate** → Q2 (what triggers leaving, Task 4) and Q3 (does leaving reduce the anomalies a bird experiences, Task 5).
- **Persist** → tested directly by the persist test (Task 8, keep-or-drop), or inferred from the absence of both an adjustment response (Q1) and a leaving response (Q2).

Do NOT use the word "suffer". We have no fitness data. "Persist" means "no detectable response"; it cannot separate "did not need to respond" from "could not respond".

**Test analyses (keep or drop).** Task 7 (long-term residency) and Task 8 (persist test) are tests. Each has explicit keep criteria. Anne and Martina decide after seeing the results; if a test fails its criteria, it is dropped from the paper.

**Core hypotheses.**
- **Resource abundance, short timescale.** Drought raises herbivore and livestock mortality, which creates carcass pulses. Predicted response: UD *contracts* during and after drought (lagged).
- **Resource depletion, long timescale.** Sustained aridity lowers herbivore populations and disperses carcasses. Predicted response: UD *expands* with chronic dryness and long-term drying.

**Out of scope. Do not analyse, do not mention in outputs:**
- Livestock density (all livestock models). This is future work.
- Soaring, thermals, uplift, flight energetics.
- Grouping species by traits or context.
- Full event-based drought analysis (separate paper). Drought enters only as a covariate (Task 2).

**Already done. Do not redo:**
- Mean NDVI vs weighted-mean NDVI (wmNDVI) within UD: no meaningful difference. Report as a robustness check only.

---

## 1. How the two agents work together

### Repo and shared folder layout

```
repo/                          (GitHub; code, docs, small tables only)
  README.md
  AGENTS.md                    (copy of this section: coordination rules)
  BRIEF.md                     (this brief)
  DECISIONS.md                 (shared, append-only decision log)
  data_contracts/              (one file per handoff table: schema, keys, units, CRS, version)
  legacy/                      (added by Anne and Martina; read-only for agents)
    models/                    Martina's June 2026 model scripts
    processing/                Anne's existing processing scripts
  processing/                  [A] scripts, numbered P00_, P01_, ...
  analysis/                    [M] scripts, numbered M00_, M01_, ...
  reports/                     checkpoint and final reports

shared_folder/                 (large data; NOT in git)
  raw/                         raw trajectories, NDVI, drought data (read-only for [M])
  derived/                     [A] writes: annotated UDs, footprints, exposure tables
  analysis_outputs/            [M] writes: model objects, figures, result tables
```

Add `.gitignore` rules so nothing from `shared_folder/` and no `.rds` larger than ~20 MB is committed.

### Rules
1. **Write only in your own area.** [A] writes in `processing/` and `shared_folder/derived/`. [M] writes in `analysis/` and `shared_folder/analysis_outputs/`. Both may append to `DECISIONS.md` and write in `reports/`.
2. **Handoffs go through data contracts.** Every table [A] delivers to [M] has a file in `data_contracts/` giving:
   - file path and version;
   - one row per what (unit of observation);
   - key columns;
   - every column with type, units and meaning;
   - CRS for spatial objects;
   - known caveats.
   Never overwrite a delivered table. Write a new version (`_v2`) and update the contract.
3. **Requests go through GitHub issues.** Label them `request:processing` (for [A]) or `request:analysis` (for [M]). State what is needed, why, and which task it blocks. Close an issue only when the deliverable and its contract are in place.
4. **Start every session** by pulling, reading `DECISIONS.md`, and checking open issues assigned to you (`gh issue list`).
5. **Branches and PRs.** Work on branches `anne/<topic>` or `martina/<topic>`. Open a PR to `main`. Anne or Martina merges; agents do not merge their own PRs.
6. **DECISIONS.md entries** carry date, author (`[A]` or `[M]`), task number, the decision and the reason.
7. **Checkpoints.** At each CHECKPOINT, stop, write the report in `reports/`, open an issue tagging both Anne and Martina, and wait.
8. **Report numbers honestly.** Flag weak or unstable results. Do not tune specifications until something is significant.

### Technical rules (both agents)
- Use `mgcv::bam()`. For the individual random effect `s(animal_id, bs="re")`: if edf ≈ 0 with `discrete=TRUE`, refit with `discrete=FALSE`. Keep the term if ΔAIC > 10 favours it.
- Do not include both species as a fixed effect and `(1|species)` / `s(species, bs="re")`. They are not identifiable together.
- In the NDVI trend model, the spatial smooth `te(longitude, latitude, month)` was settled at k = 45 after `k.check`. It is near-collinear (concurvity) with a species random intercept, so the standalone species intercept was dropped.
- Dryland classification:
  - The 4-category version (Non-dryland, Dry subhumid, Semiarid, Arid) is primary for long-term models.
  - Binary was used where species × category cells got sparse.
  - Three species occur in only one dryland category (*Gyps himalayensis*, *Torgos tracheliotus*, *Sarcoramphus papa*). They are non-identifiable in species × dryland smooths.
- Hemisphere-specific month smooths control for inverted southern seasonality.
- Plotting smooths with gratia: call `smooth_estimates()` **without** `select` and filter by the `.smooth` label afterwards. Alternatively, build `select` from `sapply(model$smooth, function(s) s$label)`, never from `names(model$coefficients)`.
- RAM is limited:
  - Save derived objects with `saveRDS()`, never `save.image()`.
  - Free large objects after extracting what you need.
  - Run long jobs detached (`nohup Rscript ... &`) with logs.
- Every script starts with a header block: purpose, inputs, outputs, decisions.

---

## 2. Task overview and dependencies

| Task | Owner | Needs | Delivers |
|---|---|---|---|
| T0a Data inventory | [A] | — | data products, UD settings, sample flow raw → UDs |
| T0b Model inventory | [M] | — | legacy model scripts, sample flow UDs → model datasets |
| T1a–b UD similarity tables | [A] | T0a | pairwise overlap table, monthly centroid/area table |
| T1c–g Residency states | [A]; thresholds agreed by Anne + Martina | T1a–b | stay blocks, regions, state table, validation |
| T1h Footprints | [A] | T1c–g | stay-block and region footprint polygons |
| **CHECKPOINT 1** | [A] writes; Anne + Martina decide thresholds | T0–T1 | agreed thresholds |
| T2a–d NDVI and drought exposure | [A] | T1h | exposure table |
| T2e Exposure diagnostics | [A] | T2a–d | diagnostics report |
| T3–T8 modelling tables | [A] | T1, T2 | one model-ready table per task, with contract |
| **HANDOFF** | [A] → [M] | all modelling tables | — |
| T3 Q1 Adjust models | [M] | T3 table | models, hypothesis table |
| T4a Q2 departure hazard models | [M] | T4 table | models |
| T4b Q2 destination choice | [A] builds used/available table (design agreed jointly); [M] fits model | T1, T2 | models |
| T5 Q3 Buffer | [A] builds experienced/foregone table; [M] runs tests | T1, T2 | results |
| T6 Age subset | [A] joins age and builds table; [M] runs models | T1, age table | results |
| T7 Long-term residency test (keep or drop) | [A] builds animal-year table; [M] runs checks | T0, T1, T2 | robustness table, keep/drop recommendation |
| T8 Persist test (keep or drop) | [A] builds drought-episode table (thresholds agreed jointly); [M] runs model | T1, T2 | results, keep/drop recommendation |
| **CHECKPOINT 2 / T9** Final report | [M]; [A] writes processing sections | all | final report |

---

## Task 0 — Inventory and sample-flow table

### T0a [A] Data inventory
1. Read the scripts in `legacy/processing/`. Identify and document in `reports/T0a_data_inventory.md`:
   - NDVI product: sensor, spatial and temporal resolution, years covered, QA/cloud handling;
   - movement trajectory format, fix rates, and how trajectories were split into months;
   - **dBBMM settings**: package and function calls (`move`/`move2`), raster resolution, margin and window sizes, location error, how monthly UDs were cut, how isopleths were derived (e.g. `getVolumeUD` → 50% and 95% contours);
   - **UD grid check**: do the monthly UD rasters of one individual already share a common grid (same origin, resolution and CRS), or does each month have its own extent and grid? Record the answer. Raster overlap in T1b depends on it;
   - which UD products are stored: rasters, isopleth polygons (50%, 95%?), summary tables;
   - how wmNDVI and meanNDVI were annotated to UDs.
2. Sample flow from raw data to annotated UDs: individuals, species and months at each step, with reasons for drops.
3. Sampling over time: individuals and UD-months per calendar year × dryland category × species, and per study/project ID. Heatmap.
4. Tracking duration per individual by species and dryland category.

**Deliverables:** `T0a_data_inventory.md`, `sample_flow_processing.csv`, `sampling_by_year.png`, `tracking_duration.csv`.

### T0b [M] Model inventory
1. Read all scripts in `legacy/models/`. Write `reports/T0b_model_inventory.md`: each script's purpose, inputs, outputs and models fitted; which objects feed which.
2. Record the residency definition currently used: distance threshold (50/100 km), anchor point, minimum duration.
3. Sample flow from annotated UDs to each model's final dataset (resident 50 km / 100 km → each model).
   - Earlier work reported inconsistent numbers (17 spp/2,420 ind; 13 spp/2,409 ind; 11 spp/28,051 rows in the change model; 164 ind in the long-term model). Reconcile them, together with T0a, and explain every drop.

**Deliverables:** `T0b_model_inventory.md`, `sample_flow_models.csv`.

---

## Task 1 — Redefine residency as a behavioural state from monthly UDs

**Goal:** replace "within N km of the tagging site" with a state classification based on each individual's own sequence of monthly UDs. This keeps seasonal migrants and relocating birds in the data instead of discarding them.

**Key constraint: two spatial scales.** Migratory birds return to the same *region* each year (e.g. Africa ↔ Europe) but often not to the same *site* within it. A single overlap threshold would wrongly label each new wintering site as a relocation or a transient range. Classification is therefore done at two scales:
- **Site** (fine scale): overlapping monthly UDs.
- **Region** (coarse scale): groups of stay blocks that are close together, used to detect migration and true relocation.

### T1a [A] UD inputs (dBBMM)
- Use the dBBMM monthly UD products identified in T0a.
- **Grid check (from T0a):** if one individual's monthly UD rasters already share a common grid, use them as they are. Only if they do not, resample them onto a common per-individual equal-area grid before computing raster overlap, and log the resolution in DECISIONS.md. Confirm with Anne before resampling.
- If only polygons are stored, use the 50% and 95% isopleth polygons.

### T1b [A] Site-scale similarity tables
- For every individual, compute similarity between all pairs of its monthly UDs.
- **Primary metric:** overlap of the **50% core isopleths**, as Jaccard (intersection / union) and intersection / smaller area. Do not use 95% isopleths as primary; at vulture range sizes they overlap almost always. Provide 95% as a sensitivity column.
- If UD rasters are usable on a common grid, also compute Bhattacharyya's coefficient on the UD probability surfaces, at least for a subset, to check agreement with the polygon metric.
- Also deliver a monthly table: UD centroid (lon/lat), 50% and 95% areas, equivalent radius `sqrt(area/pi)`.
- Work in a local azimuthal equal-area projection centred on each individual's data.

**Deliverables:** `derived/ud_pair_overlap_v1.rds` (individual, month_i, month_j, jaccard50, overlap50_min, jaccard95, BC optional), `derived/ud_month_geometry_v1.rds`, with data contracts.

**From T1c onwards, all of Task 1 is [A]'s.** Every threshold below is a starting value for exploration. [A] runs the sensitivity versions, presents them at Checkpoint 1, and implements the values Anne and Martina agree on.

### T1c [A] Stay blocks (cumulative stays longer than one month)
- A **stay block** is a run of ≥ 2 consecutive months whose UDs are linked at site scale: 50% overlap ≥ threshold with the previous month, or normalised centroid distance (centroid distance / mean equivalent radius) < 1. Allow single-month gaps (missing data) inside a block.
- Start with overlap threshold 0.3; sensitivity at 0.2 and 0.5.
- For each block, store: start and end month; duration; block centroid (UD-area-weighted mean of monthly centroids); calendar months covered.
- Months not in any block are **unlinked months**: transit, exploration, or single-month stays.

### T1d [A] Regions: group stay blocks
- Per individual, cluster stay-block centroids by distance (hierarchical, single or average linkage) at a coarse threshold `D_region`.
- **Choose `D_region` from the data, not a priori.** Plot the distribution of distances between consecutive block centroids, pooled and by species. It is expected to be bimodal: local site shifts (tens to a few hundred km) vs migratory or relocation jumps (hundreds to thousands of km). Set `D_region` at the trough. Report sensitivity at ±50%. If a species shows no clear bimodality (e.g. partially migratory species with irregular movements), report it and propose a species-specific threshold.
- Each cluster of blocks is a **region**. Label each block, and each month within it, with region ID and block ID.

### T1e [A] Two-level classification

**Region level:**

| State | Rule |
|---|---|
| **Resident** | all stay blocks in one region over the tracking period (or over the window considered) |
| **Seasonal migrant** | ≥ 2 regions occupied alternately, each re-occupied in overlapping calendar months (± 1 month tolerance) in ≥ 2 years. Exact sites within each region may change between years |
| **Regional relocation** | stays in region A, then settles in a new region B (stay blocks totalling ≥ `k_settle` months), and does not return to A within the following 12 months, outside the predictable seasonal pattern. Right-censor if tracking ends earlier |
| **Exploratory / nomadic** | many short blocks or unlinked months across regions never revisited; no stable regional pattern |

**Site level (within a region):**

| State | Rule |
|---|---|
| **Site-faithful** | revisits the same site (block monthly UDs overlap at site scale) across years within a region |
| **Site switch** | returns to the region but uses a non-overlapping site. This is NOT a relocation; record it separately as interannual site shift |

**Duration thresholds:** `k_settle` (cumulative months to count a new region as settled) at 3, 6 and 12 months.
- 3 months captures one seasonal stay.
- 6 months is the primary candidate.
- 12 months is a sensitivity bound only, because it removes migrants' seasonal ranges from "settled".

Report how the state proportions change with `k_settle` and `D_region`, by species.

**Note:** interannual site switching within a region is potentially a response variable in its own right (does a migrant change wintering site after a dry year?). Flag the number of switches by species.

### T1f [A] Individual-level behavioural summary
For each individual: proportion of months in each region-level state; number of regions; number of regional relocations; number of site switches; number of unlinked/exploratory months in the first 24 months of tracking vs later. This last measure is a proxy for the immature exploratory phase when age is unknown.

### T1g [A] Validation
- Plot ~20 individuals, stratified across:
  - resident species;
  - long-distance migrants such as *Neophron percnopterus*, including cases where wintering sites change between years;
  - partially migratory species such as *Gyps fulvus*;
  - individuals with many unlinked months.
- For each, show monthly UD outlines coloured by region and block, with the state sequence below.
- Compare the new classification with the old 50/100 km residency flag in a confusion table.

**Deliverables [A]:** `derived/stay_blocks_v1.rds`, `derived/ud_month_states_v1.rds` (one row per individual-month: block ID, region ID, region-level state, site-level state), `individual_behaviour_summary.csv`, `block_jump_distances.png`, validation figures, confusion table, with data contracts. Produce these for each sensitivity setting until thresholds are agreed; after Checkpoint 1, keep only the agreed version as `_v1` final.

### T1h [A] Footprints
From the agreed block and region tables, build:
- **stay-block footprint**: union of the monthly 50% UDs in each block;
- **region footprint**: union of all block footprints in each region.

**Deliverables:** `derived/block_footprints_v1.gpkg`, `derived/region_footprints_v1.gpkg`, with data contracts.

### CHECKPOINT 1 ([A] writes; Anne and Martina decide)
[A] writes `reports/CHECKPOINT1_report.md`. [M] adds the T0b model inventory and reviews the rest. Include:
- the reconciled sample-flow table (T0a + T0b);
- the UD grid check result and what was done about it;
- state proportions by species, `k_settle` and `D_region`;
- the block jump-distance histogram and the chosen `D_region`;
- validation plots;
- proposed values for `k_settle`, the site overlap threshold and `D_region`, with reasons, framed as options for Anne and Martina to decide;
- anything in the data that breaks the classification.

Open an issue tagging Anne and Martina. Wait for their joint decision, log it in DECISIONS.md, then build T1h with the agreed values.

---

## Task 2 — Rebuild NDVI exposure and add an independent drought index

**Why.** The current `ndvi_anomaly` is the residual of a GAM baseline (spatial-seasonal surface + dryland trends + species and individual effects). It therefore mixes four things:
1. real weather anomalies;
2. NDVI noise (cloud, compositing);
3. baseline misfit;
4. changes caused by the bird itself: when the UD moves, wmNDVI changes because the bird samples different pixels.

Component 4 makes "UD change ~ anomaly change" partly circular.

### T2a [A] Pixel-level standardized NDVI anomaly (environmental, not model-based)
- From raw monthly NDVI, compute per-pixel, per-calendar-month climatology (mean and SD) over the full NDVI record. Record the baseline period in DECISIONS.md.
- Pixel anomaly: `z = (NDVI − clim_mean) / clim_sd`. Mask pixels with very low clim_sd (bare desert) to avoid inflated z. Choose and log the threshold.

### T2b [A] Fixed-footprint exposure (removes the bird's own movement)
For each individual-month, extract the UD-weighted (or area-weighted for footprints) mean pixel anomaly over:
- **(i) the current month's UD**: experienced conditions;
- **(ii) the previous month's UD**: "what happened to the place you were". [M] uses this as the predictor in UD *change* models;
- **(iii) the current stay-block footprint** (T1h). [M] uses this for the departure hazard (T4) and for counterfactual exposure (T5). Also extract over the region footprint as a sensitivity check for migrants.

### T2c [A] Independent drought index
- **Preferred:** SPEI at 3- and 6-month timescales (SPEI-3, SPEI-6), and SPEI-12 for context. Compute from TerraClimate monthly precipitation and PET with the `SPEI` R package, or use SPEIbase (0.5°) if its coverage reaches 2025.
  - Check temporal coverage. Tracking runs to 2025.
  - Coarse resolution (0.5°) is acceptable relative to typical vulture UD sizes. Document the median UD area against the pixel area.
- **Alternative or companion:** TerraClimate monthly PDSI and climatic water deficit (`def`).
- If there is no network access, list the exact files and variables needed in an issue for Anne and continue with other tasks.
- Extract to the same footprints (i)–(iii).

### T2d [A] Drought flags
- Drought month: SPEI-6 ≤ −1 (moderate). Severe: ≤ −1.5. These are standard SPEI category cut-offs.
- **Vegetation-confirmed drought** flag: SPEI-6 ≤ −1 AND pixel NDVI anomaly < 0 over the footprint.
- Lagged versions (1, 2, 3, 6 months).

**Deliverables [A]:** `derived/exposure_month_v1.rds` (one row per individual-month; columns for each footprint × variable × lag), with data contract.

### T2e [A] Exposure diagnostics
- Correlation between the old residual anomaly, the new pixel anomaly, cumulative 3- and 6-month NDVI anomaly, and SPEI-3/6.
- How much of the old anomaly is explained by UD displacement? Regress old anomaly change on centroid displacement.
- Report whether the cumulative NDVI anomaly identifies the same months as SPEI-defined drought (agreement table).

**Deliverables [A]:** `reports/T2e_anomaly_diagnostics.md`. The drought cut-offs in T2d are confirmed by Anne and Martina after this report.

---

## Modelling tables and HANDOFF [A] → [M]

[A] builds one model-ready table per analysis, each with a data contract, in `shared_folder/derived/model_tables/`:

| Table | For | Unit | Contents |
|---|---|---|---|
| `mt_Q1_udsize_v1` | T3a | individual-month in a stay block | log UD size, the four NDVI components (new pixel anomaly, positive/negative parts), SPEI-3/6 lags, dryland (4-cat and binary), hemisphere, month, fixes per month, range type, time since tagging, species, animal_id, study ID |
| `mt_Q1_udchange_v1` | T3b | consecutive month pair within a stay block | log UD ratio, anomaly change over previous month's UD, seasonal change, all timescale versions (lags, rolling, cumulative, SPEI), same controls |
| `mt_Q2_hazard_v1` | T4a | individual-month in a stay block | departure indicators (relocation, site switch, any; seasonal-migration departures flagged), footprint anomaly lags, SPEI-6, chronic NDVI, block duration so far, controls |
| `mt_Q2_destination_v1` | T4b | relocation × candidate (1 used + 20 available) | stratum ID, used flag, NDVI anomaly, chronic NDVI, SPEI at each candidate |
| `mt_Q3_buffer_v1` | T5 | individual-month after leaving a stay block | experienced and foregone anomaly, move type, dryland, controls |
| `mt_age_v1` | T6 | individual-month (age subset) | age, age class, age at tagging, plus T4a columns |
| `mt_residency_year_v1` | T7 | animal-year | months in stay blocks, months in main region, months observed, year, time since tagging, deployment year, **long-term NDVI trend position at the bird's main region**, chronic NDVI, dryland, study ID |
| `mt_persist_v1` | T8 | individual × drought episode | episode ID, start, duration, severity (min/mean SPEI-6), outcome (adjusted / relocated / persisted), UD log ratio vs the bird's own baseline, chronic NDVI, dryland, hemisphere, month, time since tagging, species, animal_id, study ID |

Before handing off, [A] checks each table: no duplicated keys, no unexplained NAs, row counts match the sample-flow table. Then [A] opens one issue, "Modelling tables ready", listing the tables and contracts. [M] starts T3–T8 only after that issue.

---

## Task 3 [M] — Q1 "Adjust": UD size and UD change within stay blocks

**Data:** `mt_Q1_udsize_v1` and `mt_Q1_udchange_v1` from [A]. These contain months inside stay blocks only (T1c), from residents *and* from seasonal migrants' stay blocks in each region. Exclude unlinked months, because a travelling month's UD is not a home range.

### 3a. UD size (level) model
- Response: `log(UD_size)`.
- NDVI components (keep the existing four-way decomposition logic, but use the new environmental anomaly):
  - chronic baseline (individual mean);
  - long-term trend position;
  - seasonal deviation;
  - **pixel anomaly** (T2a), replacing the GAM residual.
- Include `× dryland` interactions (4-category if cells allow, else binary; log which).
- Controls:
  - hemisphere-specific cyclic month smooth;
  - sampling effort (fixes per month);
  - range type (resident's region vs migrant's seasonal region);
  - time since tagging.
- **Species structure:** a global smooth plus species deviations (Pedersen et al. 2019 "GS" model), e.g. `s(x) + s(x, species, bs = "fs", m = 1)`. Also `s(animal_id, bs = "re")`.
- **Drought asymmetry:** fit the anomaly as separate negative and positive parts, or test an asymmetric smooth. Separately, fit SPEI-6 (lag 0–3) in place of the NDVI anomaly. Compare by AIC on identical rows.

### 3b. UD change (reactivity) model
- Response: `log(UD_t / UD_{t−1})`, using only pairs of consecutive months within the same stay block.
- Predictor: change in pixel anomaly over the **previous** month's UD footprint (T2b-ii). Also the seasonal-change component.
- Timescale comparison on identical rows: concurrent, lag 2, 3 and 6, rolling 3- and 6-month means, cumulative 3 and 6, and SPEI-3/SPEI-6. Present as an AIC table.
- Same GS species structure.

### 3c. Hypothesis evaluation table
For each timescale component, report the sign and CI of the global effect, overall and by dryland category, against these predictions:
- abundance: drought → contraction (short, lagged);
- depletion: chronic or long-term drying → expansion.

State plainly which predictions are supported, contradicted, or unresolved.

### 3d. Species outliers
List species whose deviation smooths differ strongly from the global smooth. Include n individuals, n regions/studies and n months for each. In particular, check whether the *Torgos tracheliotus* positive effect rests on one study or region and how many individuals. Do not write ecological explanations; just report support.

**Deliverables [M]:** `analysis/M03_Q1_adjust_models.R`, model summaries, `k.check` and concurvity output, AIC tables, figures (global + species smooths), `Q1_hypothesis_table.csv`.

---

## Task 4 — Q2 "Relocate": what triggers leaving, and where do birds go?

### T4a [M] Discrete-time departure hazard
- Data: `mt_Q2_hazard_v1` from [A].
- Unit: individual-month within a stay block (T1c).
- Response: 1 if the stay block ends this month; 0 otherwise. Fit separate models for:
  - (a) **regional relocation** only (T1e);
  - (b) **site switch** within a region (migrants changing site between years, residents shifting site);
  - (c) any departure.
  Exclude predictable seasonal migration departures from (a) and (c), or model them separately, so that migration timing does not dominate.
- Predictors:
  - pixel NDVI anomaly over the stay-block footprint (T2b-iii), lagged 1–3 months;
  - SPEI-6 over the footprint;
  - chronic NDVI;
  - dryland category;
  - block duration so far (`s()` smooth);
  - hemisphere month smooth;
  - time since tagging.
- Random effects: species and individual.
- Use binomial `bam()` with cloglog link. Report the number of events. If there are fewer than ~10 events per parameter, simplify the model and say so.

### T4b Destination choice ([A] builds table, [M] fits model)
- Design (agreed by Anne and Martina before [A] builds the table):
  - origin: stay-block footprint during the 3 months before departure;
  - destination: new block footprint during the 3 months after settling;
  - 20 random available destinations per relocation, at the same displacement distance (random bearings, restricted to land).
- [A] generates the available destinations as polygons (origin footprint shifted to each random location) and extracts NDVI anomaly, chronic NDVI and SPEI for origin, destination and available footprints.
- [M] fits a conditional logistic regression (used vs available) with `survival::clogit`.
- Report: do birds move toward greener or less drought-affected areas, or just away?

**Deliverables:** [A] `mt_Q2_destination_v1` with contract; [M] `analysis/M04_Q2_relocation.R`, event counts by species and dryland, model outputs, figures.

---

## Task 5 — Q3 "Relocate: does it pay off?": does moving reduce the anomalies a bird experiences?

- [A] For every month after a bird leaves a stay block (regional relocation, site switch, or exploratory months), extract:
  - **experienced** anomaly: over the current UD;
  - **foregone** anomaly: over the stay-block footprint it left, for the same month.
- Do not treat a migrant's predictable seasonal move as "leaving". For migrants, compare the site used this year with the site(s) used in the same region in previous years (site switch as a buffering move).
- [M] tests whether experienced anomalies are less negative and less variable than foregone ones (paired comparison within individual; mixed model with species and individual random effects).
- [M] contrasts with persisting individuals: within-individual variance of experienced anomaly in residents vs movers, **always benchmarked against the foregone exposure**, not raw.
- Report by dryland category.

**Deliverables:** [A] `mt_Q3_buffer_v1` with contract; [M] `analysis/M05_Q3_buffering.R`, results table, figure.

---

## Task 6 — Age subset (Europe) ([A] steps 1–2, [M] steps 3–4)

An age table for a subset of European individuals will be provided by a collaborator. [A] checks for it in `shared_folder/raw/`; if absent, prepare the code and skip.

1. [A] Join by individual ID and build `mt_age_v1`. Report coverage:
   - individuals, species and months with known age;
   - age at tagging distribution (nestlings vs adults). For nestlings, age = time since tagging;
   - **months per dryland category (4-cat)**. Much of Iberia is dryland (dry subhumid and semiarid), so do not assume this subset is non-dryland.
2. [A] **Stop and check with Anne and Martina before steps 3–4.** Known likely limits to confirm against the coverage table: few months in the Arid category, and a European context with supplementary feeding, so the subset does not stand in for African or Asian arid systems.
3. [M] Validate the T1f proxy: does the exploratory phase coincide with immature ages? Use species-appropriate adult age (e.g. ~5 years for *Gyps fulvus* / *Aegypius monachus*; check the literature value for each species present and log it).
4. [M] Refit the T4a departure hazard on this subset with age (or age class) added. Report how much the NDVI and drought effects change when age is controlled. Where dryland coverage allows, test whether the age effect differs between dryland and non-dryland months.

**Deliverables:** [A] `mt_age_v1` with contract and coverage table; [M] `analysis/M06_age_subset.R`, outputs, short note.

---

## Task 7 — Long-term residency test: keep or drop ([A] builds `mt_residency_year_v1`, [M] runs steps 1–6)

If kept, this becomes the long-timescale part of Q2 (relocate).

Earlier result: monthly residency probability declined strongly over 20 years in Semiarid (β ≈ −0.19) and Arid (β ≈ −0.41) categories, and was flat in Dry subhumid. This came from a binomial GLMM of resident months per animal-year, using the old km-based residency.

Before this result can be used:
1. Refit using the new states (proportion of stay-block months per animal-year, and proportion in the individual's main region).
2. Add **time since tagging** (and age where known) to separate calendar trend from individuals drifting away from capture sites and maturing. Include deployment year (cohort) and calendar year together only if identifiable; explain the age-period-cohort limitation.
3. **Leave-one-study-out:** refit dropping each study/project in turn. Report whether the arid/semiarid decline holds.
4. **Sampling imbalance:** weight animal-years equally, and refit restricted to years with ≥ N individuals per dryland category. Choose N from the T0a heatmap.
5. **Test against productivity.** Refit with the long-term NDVI trend position at the bird's main region as the predictor (in place of, or alongside, calendar year), with the same controls. A decline over calendar years alone does not answer the depletion hypothesis.
6. Conclude in one paragraph: robust, fragile, or artefact.

**Keep criteria (all must hold):** the decline survives steps 2–4 in direction and remains clearly different from zero, AND residency is related to the long-term NDVI trend in step 5. Otherwise recommend dropping it. Anne and Martina decide.

**Deliverables [M]:** `analysis/M07_residency_trend_checks.R`, table of estimates across all checks, keep/drop recommendation.

---

## Task 8 — Persist test: keep or drop ([A] builds `mt_persist_v1`, [M] runs the model)

**Goal:** test "persist" directly by making the strategy itself the response: during a drought episode, does the bird adjust, relocate, or neither?

1. [A] **Define drought episodes** within stay blocks: consecutive months with SPEI-6 ≤ the agreed drought cut-off (T2d), with a minimum duration. Propose options for the minimum duration (e.g. 1, 2, 3 months) with the number of episodes each produces by dryland category and species.
2. [A] **Classify each episode's outcome**, within the episode and up to the agreed lag (e.g. 3 months) after it:
   - **Relocated:** the stay block ends (regional relocation or site switch, excluding predictable seasonal migration).
   - **Adjusted:** did not leave, and UD size changed beyond a threshold relative to the bird's own baseline (e.g. |log UD ratio| larger than 1 SD of its non-drought month-to-month changes).
   - **Persisted:** neither.
   Propose options for the adjustment threshold and lag, with how the outcome proportions change under each. Anne and Martina decide.
3. [M] **Model the outcome** (three categories) as a function of drought severity, dryland category and chronic NDVI, with species and individual random effects. Options: `mgcv::gam(family = multinom(K = 2))` with random effects as `s(..., bs = "re")`, or `brms` with a categorical family if `mgcv` cannot handle the random effects. Report the number of episodes per outcome × dryland category.

**Keep criteria (all must hold):**
- enough episodes in each outcome category to estimate the model (as a guide, ≥ 30 per outcome overall, and dryland categories with too few episodes merged or dropped);
- outcome proportions are reasonably stable across the threshold options (no change in which strategy dominates in each dryland category);
- the result adds something beyond Q1 and Q2 (e.g. strategy proportions differ between dryland categories or with drought severity).

Otherwise recommend dropping it; "persist" is then inferred from the absence of Q1 and Q2 effects. Anne and Martina decide.

**Deliverables:** [A] `mt_persist_v1` with contract and threshold sensitivity table; [M] `analysis/M08_persist_test.R`, results, keep/drop recommendation.

---

## Task 9 — Final report

### CHECKPOINT 2
[M] writes `reports/FINAL_report.md`; [A] writes the processing and classification sections and checks the methods. Include:
- one paragraph per question (Q1 adjust; Q2 and Q3 relocate; persist, from Task 8 or from the absence of Q1 and Q2 effects) with the key numbers;
- the hypothesis table from 3c;
- robustness summary for Task 6;
- keep/drop recommendations for Tasks 7 and 8, with the evidence for each;
- list of results that are weak, unstable or contradictory;
- a list of figures, with a proposed main/supplement split;
- methods paragraph drafts for: UD processing and annotation ([A]); residency classification ([A]); NDVI/SPEI exposure ([A]); models ([M]).

Do not write introduction or discussion text. Open an issue tagging Anne and Martina.
