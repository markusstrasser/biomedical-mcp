"""MyGene.info REST API client — gene annotations aggregator."""

import hashlib
import logging

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

MYGENE_BASE = "https://mygene.info/v3"

DEFAULT_FIELDS = (
    "symbol,name,entrezgene,ensembl.gene,uniprot,pathway.kegg,"
    "go,genomic_pos,alias,type_of_gene"
)


def _is_retryable(exc: BaseException) -> bool:
    if isinstance(exc, httpx.HTTPStatusError):
        return exc.response.status_code in (429, 500, 502, 503, 504)
    return isinstance(exc, (httpx.ConnectError, httpx.ReadTimeout))


class MyGene:
    def __init__(self, cache: Cache):
        self.cache = cache
        self.client = httpx.Client(
            base_url=MYGENE_BASE,
            timeout=30,
            headers={"User-Agent": "biomedical-mcp/0.1"},
        )

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=2, max=15), retry=retry_if_exception(_is_retryable))
    def _get(self, path: str, params: dict | None = None) -> dict:
        resp = self.client.get(path, params=params)
        resp.raise_for_status()
        return resp.json()

    def _cached(self, prefix: str, path: str, params: dict | None = None, max_age_days: int = 30) -> dict:
        raw = f"{path}:{sorted(params.items()) if params else ''}"
        key = f"mygene:{prefix}:{hashlib.md5(raw.encode()).hexdigest()}"
        cached = self.cache.get(key, max_age_days=max_age_days)
        if cached is not None:
            return cached
        result = self._get(path, params)
        self.cache.set(key, result)
        return result

    def gene_info(self, gene_id: str) -> dict | None:
        """Get gene info by symbol, Entrez ID, or Ensembl ID.

        Detects ID type: ENSG* -> Ensembl, numeric -> Entrez, else -> symbol query.
        """
        # Resolve to Entrez ID if needed
        if gene_id.startswith("ENSG"):
            data = self._cached("query", "/query", {"q": f"ensembl.gene:{gene_id}", "species": "human", "fields": DEFAULT_FIELDS, "size": "1"})
            hits = data.get("hits", [])
            if not hits:
                return None
            return self._extract_gene(hits[0])
        elif gene_id.isdigit():
            data = self._cached("gene", f"/gene/{gene_id}", {"fields": DEFAULT_FIELDS})
            return self._extract_gene(data)
        else:
            # Symbol lookup
            data = self._cached("query", "/query", {"q": f"symbol:{gene_id}", "species": "human", "fields": DEFAULT_FIELDS, "size": "1"})
            hits = data.get("hits", [])
            if not hits:
                return None
            return self._extract_gene(hits[0])

    def search(self, query: str, species: str = "human", limit: int = 10) -> list[dict]:
        """Search genes by keyword."""
        data = self._cached("search", "/query", {
            "q": query,
            "species": species,
            "fields": "symbol,name,entrezgene,ensembl.gene",
            "size": str(min(limit, 50)),
        })
        hits = data.get("hits", [])
        return [
            {
                "entrez_id": str(h.get("entrezgene", h.get("_id", ""))),
                "symbol": h.get("symbol"),
                "name": h.get("name"),
                "ensembl_gene": self._first_ensembl(h.get("ensembl")),
            }
            for h in hits
        ]

    def _extract_gene(self, data: dict) -> dict:
        """Extract structured gene info from raw MyGene response."""
        result = {
            "entrez_id": str(data.get("entrezgene", data.get("_id", ""))),
            "symbol": data.get("symbol"),
            "name": data.get("name"),
            "type_of_gene": data.get("type_of_gene"),
            "aliases": data.get("alias", []) if isinstance(data.get("alias"), list) else [data["alias"]] if data.get("alias") else [],
        }

        # Ensembl
        ensembl = data.get("ensembl")
        result["ensembl_gene"] = self._first_ensembl(ensembl)

        # UniProt
        uniprot = data.get("uniprot")
        if isinstance(uniprot, dict):
            swiss = uniprot.get("Swiss-Prot")
            result["uniprot"] = swiss if isinstance(swiss, str) else swiss[0] if isinstance(swiss, list) and swiss else None
        else:
            result["uniprot"] = None

        # KEGG pathways
        kegg = data.get("pathway", {}).get("kegg") if isinstance(data.get("pathway"), dict) else None
        if kegg:
            if isinstance(kegg, list):
                result["kegg_pathways"] = [{"id": p.get("id"), "name": p.get("name")} for p in kegg[:15]]
            elif isinstance(kegg, dict):
                result["kegg_pathways"] = [{"id": kegg.get("id"), "name": kegg.get("name")}]
        else:
            result["kegg_pathways"] = []

        # GO terms (grouped by category)
        go = data.get("go")
        if isinstance(go, dict):
            go_summary = {}
            for category in ("BP", "CC", "MF"):
                terms = go.get(category, [])
                if isinstance(terms, dict):
                    terms = [terms]
                if isinstance(terms, list):
                    go_summary[category] = [
                        {"id": t.get("id"), "term": t.get("term")}
                        for t in terms[:10]
                    ]
            result["go_terms"] = go_summary
        else:
            result["go_terms"] = {}

        # Genomic position
        gpos = data.get("genomic_pos")
        if isinstance(gpos, dict):
            result["genomic_position"] = {
                "chr": gpos.get("chr"),
                "start": gpos.get("start"),
                "end": gpos.get("end"),
                "strand": gpos.get("strand"),
            }
        elif isinstance(gpos, list) and gpos:
            g = gpos[0]
            result["genomic_position"] = {
                "chr": g.get("chr"),
                "start": g.get("start"),
                "end": g.get("end"),
                "strand": g.get("strand"),
            }

        return result

    @staticmethod
    def _first_ensembl(ensembl) -> str | None:
        """Extract first Ensembl gene ID (can be str, dict, or list)."""
        if isinstance(ensembl, str):
            return ensembl
        if isinstance(ensembl, dict):
            return ensembl.get("gene")
        if isinstance(ensembl, list) and ensembl:
            first = ensembl[0]
            return first.get("gene") if isinstance(first, dict) else str(first)
        return None
