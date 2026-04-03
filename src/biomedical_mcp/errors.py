"""Structured error types with actionable suggestions.

Pattern from BioMCP: errors include the next step, not just what went wrong.
"""

from __future__ import annotations


class BiomedicalError(Exception):
    """Base for all biomedical-mcp errors."""

    def to_dict(self) -> dict:
        return {"error": str(self), "type": type(self).__name__}


class NotFoundError(BiomedicalError):
    """Entity not found — includes suggested next action."""

    def __init__(self, entity: str, query: str, suggestion: str):
        self.entity = entity
        self.query = query
        self.suggestion = suggestion
        super().__init__(f"{entity} '{query}' not found. {suggestion}")

    def to_dict(self) -> dict:
        return {
            "error": f"{self.entity} '{self.query}' not found",
            "suggestion": self.suggestion,
            "type": "NotFoundError",
        }


class ApiKeyRequiredError(BiomedicalError):
    """API key needed — tells user exactly which env var to set."""

    def __init__(self, api: str, env_var: str, docs_url: str = ""):
        self.api = api
        self.env_var = env_var
        self.docs_url = docs_url
        msg = f"{api} requires {env_var} environment variable."
        if docs_url:
            msg += f" See: {docs_url}"
        super().__init__(msg)

    def to_dict(self) -> dict:
        return {
            "error": f"{self.api} requires API key",
            "env_var": self.env_var,
            "docs_url": self.docs_url,
            "type": "ApiKeyRequiredError",
        }


class SourceUnavailableError(BiomedicalError):
    """Source API is down or unreachable — suggests alternatives."""

    def __init__(self, source: str, reason: str, suggestion: str = ""):
        self.source = source
        self.reason = reason
        self.suggestion = suggestion
        msg = f"{source} unavailable: {reason}"
        if suggestion:
            msg += f". Try: {suggestion}"
        super().__init__(msg)

    def to_dict(self) -> dict:
        return {
            "error": f"{self.source} unavailable",
            "reason": self.reason,
            "suggestion": self.suggestion,
            "type": "SourceUnavailableError",
        }
