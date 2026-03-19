# Bio MCP — Human TODOs

Things that need your action (not automatable by agents).

## Required

- [ ] **Register PharmVar API key** — Go to [pharmvar.org](https://www.pharmvar.org), create account, generate API key in Account Settings. Then set `PHARMVAR_API_KEY` env var in your shell profile or `.mcp.json` env block. Without this, `drugs_star_alleles` returns auth errors.

## Recommended

- [ ] **Restart Claude Code** (or `/mcp` reset) — The MCP server was reinstalled but Claude Code caches the tool list. Restart the session or run `/mcp` to pick up the new `mcp__bio__*` tools (population, panels, phenotype, expression, gwas, literature, composite).

- [ ] **Bump version to 0.4.0** — The expansion added 13 APIs and 6 domains. Consider bumping in `pyproject.toml` when you're satisfied with testing.

## Optional / Future

- [ ] **Review composite tool output** — Try `composite_variant_context("1-11796321-G-A", gene_symbol="MTHFR")` and `composite_gene_dossier("BRCA1")` in a session. Check if the provenance structure and concurrent dispatch work as expected.

- [ ] **Evaluate Tier 1 usage** after 10+ review sessions — Are gnomAD, PanelApp, HPO/Monarch the most-used new tools? If so, expand to Tier 2 (GTEx, HGNC, GWAS). If not, focus on improving Tier 1.

- [ ] **Set up CI canary tests** — The 24 new smoke tests hit live APIs. Consider a weekly GitHub Actions run to catch API breakage early (contract tests per model-review recommendation).
