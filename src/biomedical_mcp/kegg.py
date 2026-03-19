"""KEGG REST API — pathway, compound, and cross-reference lookups."""

import hashlib
import logging
import re

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

from biomedical_mcp.cache import Cache

log = logging.getLogger(__name__)

KEGG_URL = "https://rest.kegg.jp"


def _transient(exc: BaseException) -> bool:
    return isinstance(exc, (httpx.TimeoutException, httpx.HTTPStatusError))


class KEGG:
    def __init__(self, cache: Cache):
        self.cache = cache
        self.client = httpx.Client(timeout=30.0)

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=1, max=10),
           retry=retry_if_exception(_transient))
    def _get_text(self, path: str) -> str:
        resp = self.client.get(f"{KEGG_URL}/{path.lstrip('/')}")
        resp.raise_for_status()
        return resp.text

    def _parse_flat(self, text: str) -> list[dict]:
        """Parse KEGG tab-delimited list output."""
        results = []
        for line in text.strip().split("\n"):
            if "\t" in line:
                parts = line.split("\t", 1)
                results.append({"id": parts[0].strip(), "description": parts[1].strip() if len(parts) > 1 else ""})
        return results

    def gene_pathways(self, gene_symbol: str, organism: str = "hsa") -> dict:
        """Get KEGG pathways for a human gene."""
        key = hashlib.md5(f"kegg:gene_pathways:{gene_symbol}:{organism}".encode()).hexdigest()
        cached = self.cache.get(key, max_age_days=14)
        if cached is not None:
            return cached

        # First find the KEGG gene ID — match exact symbol in the description field
        text = self._get_text(f"find/genes/{gene_symbol}")
        gene_id = None
        for line in text.strip().split("\n"):
            if not line.startswith(f"{organism}:"):
                continue
            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue
            # Description field has format: "SYMBOL, ALIAS1, ALIAS2; description"
            desc = parts[1]
            symbols = desc.split(";")[0].split(",")
            symbols = [s.strip().upper() for s in symbols]
            if gene_symbol.upper() in symbols:
                gene_id = parts[0].strip()
                break

        if not gene_id:
            return {"error": f"Gene not found in KEGG: {gene_symbol}", "organism": organism}

        # Get pathways linked to this gene
        pathway_text = self._get_text(f"link/pathway/{gene_id}")
        pathways = []
        pathway_ids = []
        for line in pathway_text.strip().split("\n"):
            if "\t" in line:
                pathway_id = line.split("\t")[1].strip()
                pathway_ids.append(pathway_id)

        # Get pathway names via get endpoint (handles path: prefix)
        for pid in pathway_ids[:25]:
            clean_id = pid.replace("path:", "")
            try:
                entry_text = self._get_text(f"get/{clean_id}")
                name = ""
                for line in entry_text.split("\n"):
                    if line.startswith("NAME"):
                        name = line.replace("NAME", "").strip()
                        name = re.sub(r"\s*-\s*Homo sapiens.*$", "", name)
                        break
                pathways.append({"pathway_id": clean_id, "name": name})
            except Exception:
                pathways.append({"pathway_id": clean_id, "name": ""})

        result = {"gene": gene_symbol, "kegg_gene_id": gene_id, "organism": organism,
                  "count": len(pathways), "pathways": pathways}
        self.cache.set(key, result)
        return result

    def pathway_info(self, pathway_id: str) -> dict:
        """Get detailed info for a KEGG pathway."""
        key = hashlib.md5(f"kegg:pathway:{pathway_id}".encode()).hexdigest()
        cached = self.cache.get(key, max_age_days=14)
        if cached is not None:
            return cached

        text = self._get_text(f"get/{pathway_id}")
        info: dict = {"pathway_id": pathway_id}
        current_section = ""
        genes = []
        compounds = []
        diseases = []

        for line in text.split("\n"):
            if line and not line[0].isspace():
                current_section = line.split()[0] if line.split() else ""
                if current_section == "NAME":
                    info["name"] = line.replace("NAME", "").strip()
                elif current_section == "DESCRIPTION":
                    info["description"] = line.replace("DESCRIPTION", "").strip()
                elif current_section == "CLASS":
                    info["class"] = line.replace("CLASS", "").strip()
                elif current_section == "ORGANISM":
                    info["organism"] = line.replace("ORGANISM", "").strip()
            elif current_section == "GENE" and line.strip():
                parts = line.strip().split(None, 1)
                if len(parts) >= 2:
                    gene_entry = parts[1].split(";")
                    genes.append({"id": parts[0], "symbol": gene_entry[0].strip(),
                                  "description": gene_entry[1].strip() if len(gene_entry) > 1 else ""})
            elif current_section == "COMPOUND" and line.strip():
                parts = line.strip().split(None, 1)
                if len(parts) >= 2:
                    compounds.append({"id": parts[0], "name": parts[1].strip()})
            elif current_section == "DISEASE" and line.strip():
                parts = line.strip().split(None, 1)
                if len(parts) >= 2:
                    diseases.append({"id": parts[0], "name": parts[1].strip()})

        if genes:
            info["genes"] = genes[:50]
            info["gene_count"] = len(genes)
        if compounds:
            info["compounds"] = compounds
        if diseases:
            info["diseases"] = diseases

        self.cache.set(key, info)
        return info

    def find(self, query: str, database: str = "pathway") -> dict:
        """Search KEGG database. database: pathway, genes, compound, disease, drug."""
        key = hashlib.md5(f"kegg:find:{database}:{query}".encode()).hexdigest()
        cached = self.cache.get(key, max_age_days=7)
        if cached is not None:
            return cached

        text = self._get_text(f"find/{database}/{query}")
        results = self._parse_flat(text)
        result = {"query": query, "database": database, "count": len(results), "results": results[:50]}
        self.cache.set(key, result)
        return result
