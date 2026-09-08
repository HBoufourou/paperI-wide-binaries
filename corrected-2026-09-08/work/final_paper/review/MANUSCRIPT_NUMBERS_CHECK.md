# Manuscript numerical editorial check

Status: PASS; 204 checks, 0 discrepancies.

Manuscript SHA-256: `78d7c4f10c9069c6d69f508204644e67d399bd3c7f8b42456cc3893bd38e5377`.
Checker SHA-256: `2eebe477f2f8d7e7ffc654548eff6bb63d760d5746611d27690905ab9eb2d3c2`.

Editorial agreement of the abstract, all numerical results in section5, conclusion, and relevant design/old-new-count accounting with saved CSV/JSON. No production import or scientific reconstruction.

No manuscript text, frozen source, protocol, score or decision was edited. The companion JSON records every comparison and all input hashes.

The two result tables were parsed and all80 numerical cells compared. Counts and transitions require exact integer agreement; displayed percentages/ranges are checked at their stated decimal precision. The four-decimal mean uses0.00005 percentage-point tolerance; three-decimal interval endpoints use0.0005 points. The effective size and maximum weight use half a last displayed decimal. Text anchors normalize whitespace only, to tolerate the concurrent editorial spacing correction. The first anchor-only failure is preserved in MANUSCRIPT_NUMBERS_CHECK_pre_spacing_normalization.json/.md; it contained no numerical discrepancy.

The declared501 combined full66-envelope passes equal257+244 across realizations. Its classification Wilson counts are266+265=531. Both retention and classification labels match their CSV columns. The diagnostic contrast has451 correct gains minus440 losses=11 net, but148 wrong gains minus2 losses=146 additional wrong views. Full34 envelope loses973−123=850 correct views net and removes106−6=100 wrong views net; conditional loses111−70=41 correct and removes57 wrong. All are paired catalogue views, not independent astrophysical systems.

Historical count reconciliation: the earlier independent population proof covers441024 rows =440000 scientific parents in22 libraries +1024 trajectory fixtures. The final experiment uses880000 scientific parents =640000 fresh +240000 from12 historical inputs; it does not add all22 earlier libraries or include the old fixtures.

The abstract and conclusion retain the distinction between observed95%/1% classification, marginal Wilson gates and99% true-model retention. No numerical discrepancy was found in the quoted extrema, identities, ranges, signed means or transitions. A physical or inference PASS is not inferred from this editorial check.

## Discrepancies

None.

## Hashed inputs

- `MANUSCRIPT_REVIEW_VERSION.md`: `78d7c4f10c9069c6d69f508204644e67d399bd3c7f8b42456cc3893bd38e5377`
- `stats\results\metrics.csv`: `160619d0bb1769249d97fda28bb39f9d57432ea638f61f5999db353654f1ff5b`
- `stats\results\paired.csv`: `d07d7ba6aad16333537c27ff699230ba5e85474ed78694250643054dee0fb0b6`
- `stats\results\method_summary.csv`: `9efbc6e10dcd2f6a902055b3edacadecf3997a77645320fb15e0515f1b222b66`
- `stats\results\information_summary.csv`: `66869743f2ebddbd2f63593e505d2a69047e95852ab51c227f3211addb55e177`
- `stats\results\envelope_summary.csv`: `74fbdc30258f6a7d4329352f3a84a247fba3c65751d9680ca1882a3c206bd4ac`
- `stats\results\envelope_effects.csv`: `dfd458d642b17fe313bf6ac508cf6621d69c48f1bec235ac78b21ed0b1cc6c6a`
- `stats\results\realization_summary.csv`: `4f0f74b4f7d43e4f526bbab3b9ef011e5b363e229d5264d3db4dbe569c383256`
- `stats\results\library_support.csv`: `c9c6b09b8c38abb5c8778bb472685727383619427a46115e1159550d478b1cd5`
- `stats\results\calibration\manifest.json`: `3fa529bdaaca7e03ea9d0f62fef383c8280bc08fa89db109d08b76f8ae785235`
- `stats\results\test\manifest.json`: `688fe86ac232adc9335737ef37b3562e94694ed44fe88892a3a7192a0e96fdf5`
- `stats\results\EXECUTION_PROTOCOL.json`: `a31ae2e50feb7af250d9a8d9f8c61713005fabde04d096d1f570417c4db53685`
- `review/physical_banks_verification.json`: `006f16b21c9685858652b5d4af0e8d7a9279838c9439f607e5fadefd62b5ac8a`
- `review/delivery_api_verification.json`: `848be529d6729a2084f65e086d039e1ea00189ee3ae50ae83bc29fa05e3a774b`
- `../dr4_chain/review/population_verification.json`: `41db1f71cd144fd28022e03fff1adeba6cc65f8f1eb754235b9565d60a76c56c`
- `dr4_chain\decision\results\metrics.csv`: `e54fb60b79ee37033fb1d3742c50c6f97c613f8e8bf6cf084caa5bb24fa3075c`
