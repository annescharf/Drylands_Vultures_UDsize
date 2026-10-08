# AIMS: what this paper is for

Keep this short document as the reference point. Every analysis must answer one of the questions below. If it doesn't, it doesn't go in this paper.

Confidence tags:
- **[Certain]**: shown by results we already have.
- **[Likely]**: strong expectation from current results or the literature.
- **[Guessing]**: plausible, but no direct evidence yet.

---

## Overall aim

**Test whether vultures respond to changes in vegetation productivity by adjusting their space use, relocating, or persisting, and whether the response follows the resource-abundance (short-term drought → carcass pulse → smaller ranges) or resource-depletion (long-term drying → fewer, dispersed carcasses → larger ranges) hypothesis, across aridity levels.**

What is new compared with the literature [Likely]:
- Earlier vulture tracking studies are single-region and single-species, and explain range size mainly by season, sex, age and breeding. We test a continuous environmental driver (productivity) across ~13 species and five realms with one consistent method.
- Single-region studies found support for either abundance or depletion, never both. We test both, separated by timescale.
- Vultures respond to vegetation only indirectly, through carrion. This tests resource-tracking ideas, which were developed mostly for herbivores, one trophic step removed.

---

## How the three strategies map onto the questions

| Strategy | Tested by |
|---|---|
| **Adjust** (change UD size, stay in the same area) | Q1 |
| **Relocate** (leave and settle elsewhere) | Q2 (what triggers leaving, where birds go) and Q3 (does leaving reduce the anomalies a bird experiences). If the long-term residency test is kept, it adds the long-timescale part of Q2 |
| **Persist** (stay without detectable adjustment) | Directly by the persist test, if kept. Otherwise inferred from the absence of both an adjustment response (Q1) and a leaving response (Q2) |

"Persist" means "no detectable response". Without fitness data it cannot separate "did not need to respond" from "could not respond", and the paper says so.

---

## Questions

### Q1 — Adjust: does UD size follow productivity, and on which timescale?
- **Abundance prediction:** drought months (negative anomaly, lagged) → UD contracts.
- **Depletion prediction:** chronically dry sites and long-term drying → UD expands.
- **Tasks:** T2 (new NDVI anomaly, SPEI drought index), T3 (UD size and UD change models).
- **Likely result:**
  - Chronically drier sites have larger UDs. **[Certain]** in previous models (dryland UDs ~34% larger; drier still within drylands). **[Likely]** to hold in the new models. Supports depletion at the long timescale.
  - Short-term drought response. Direction **[Guessing]**. In earlier models, drier anomalies went with *larger* UDs overall, which goes against abundance. A contraction in drylands appeared only in the livestock model, now dropped. The new environmental anomaly removes the part caused by the bird's own movement, so the short-term effect may get **weaker [Likely]**.
  - The response is lagged by ~3 months **[Certain]** for raw NDVI. Whether this holds for the new anomaly or SPEI is **[Guessing]**.
  - Species differ a lot **[Certain]**. A global average effect will exist but explain little of month-to-month change **[Likely]** (deviance explained ~8% before).

### Q2 — Relocate (1): do poor conditions make vultures leave, and where do they go?
- **Tasks:** T1 (residency states: stay blocks, regions, relocations, site switches), T4a (departure hazard), T4b (destination choice).
- **Likely result:**
  - Regional relocations will be rare **[Likely]**. Site switches within a region, especially among migrants, will be more common **[Guessing]**.
  - Negative anomalies or drought raise the chance of leaving: **[Guessing]**, probably a weak effect.
  - In birds of known age, much of the leaving will be explained by being young **[Likely]**.
  - Birds move toward greener or less dry areas: **[Guessing]**.

### Q3 — Relocate (2): does leaving pay off, i.e. does moving reduce the anomalies a bird experiences?
- **Tasks:** T2b (fixed-footprint exposure), T5 (experienced vs foregone anomaly).
- **Likely result:**
  - Birds that leave experience less negative anomalies than at the place they left: **[Guessing]**. If it holds, it is the most original result of the paper and the strongest case for Ecology Letters.
  - If no buffering is found, Q3 becomes a short null result supporting Q2.

---

## Test analyses: keep or drop

These are run, then Anne and Martina decide. If a test fails its keep criteria, it is dropped from the paper.

### Long-term residency test (T7) → would extend Q2 (relocate)
- **Earlier result [Certain]:** residency declined over 20 years in semiarid and arid sites under the old km-based definition.
- **Keep if both hold:** (1) the decline survives the controls for time since tagging, leave-one-study-out and sampling balance; (2) residency is related to the long-term NDVI trend at the bird's site, not only to calendar year.
- **Likely result:** the decline gets weaker once time since tagging and sampling imbalance are controlled **[Guessing]**. Whether it links to the long-term NDVI trend is **[Guessing]**.

### Persist test (T8) → tests persist directly
- **What it does:** for each drought episode a bird experiences in a stay block, classifies the outcome as adjusted, relocated or persisted, and models which outcome happens depending on drought severity, dryland category and chronic NDVI.
- **Keep if all hold:** enough episodes in each outcome; outcome proportions stable across the threshold options; the result adds something beyond Q1 and Q2 (e.g. strategy proportions differ between dryland categories or with drought severity).
- **Likely result:** persist is the most common outcome **[Guessing]**, because relocations are rare and month-to-month UD responses are noisy. Whether the proportions differ with aridity is **[Guessing]**.

---

## How tasks map to questions

| Task | Serves | Needed for the paper? |
|---|---|---|
| T0 Inventory, sample flow | all | Yes. Reviewers will need the sample-flow table |
| T1 Residency states | Q1 (which months count), Q2, Q3, tests | Yes. Everything depends on it |
| T2 NDVI anomaly + SPEI | Q1, Q2, Q3, tests | Yes |
| T3 UD size and change models | Q1 (adjust) | Yes (core result) |
| T4 Departure and destination | Q2 (relocate) | Yes, but shrink to descriptive if events are too few |
| T5 Experienced vs foregone | Q3 (relocate) | Yes, quick to run |
| T6 Age subset | Robustness for Q2 | Only as a robustness check |
| T7 Long-term residency test | Q2 (relocate, long timescale) | Keep or drop |
| T8 Persist test | Persist | Keep or drop |

**Minimum paper:** Q1 plus either Q2 or Q3. **Full paper:** Q1, Q2 and Q3, plus T7 and T8 if they pass their keep criteria.

---

## Rabbit-hole rules (for people and agents)

1. **No new questions.** Anything that doesn't answer Q1, Q2, Q3 or the two listed tests goes to a "future work" list in `DECISIONS.md`, not into the analysis.
2. **No new data sources** beyond what is listed: tracking data, UDs, NDVI, SPEI or TerraClimate drought, and the European age subset.
3. **Out of scope:** livestock, soaring and thermals, grouping species by traits, full event-based drought analysis.
4. **No single-species stories.** Report species deviations and their sample support; don't build explanations for individual species.
5. **Stop-and-ask triggers.** Stop and open an issue when:
   - a result contradicts a [Certain] item above;
   - a task needs more than one extra round of model specification to "work";
   - a sensitivity check changes a conclusion's direction.
6. **Thresholds are decided by Anne and Martina.** Agents propose options, each with the reasons for it and the likely consequences of choosing it (e.g. how many individuals or months change state, which results are sensitive to it). Agents never fix a threshold on their own.
7. **Null results are results.** A clean null for Q2 or Q3 is reported as such; don't add predictors until something becomes significant.
8. **Tests stay tests.** If T7 or T8 fails its keep criteria, drop it. Don't rework it until it passes.
