from __future__ import annotations

from urllib.error import HTTPError, URLError

import pytest

from genomics import gnomad


SUCCESS_RESPONSE = {
    "data": {
        "region": {
            "variants": [
                {
                    "chrom": "1",
                    "pos": 100,
                    "ref": "A",
                    "alt": "G",
                }
            ]
        }
    }
}


def disable_retry_wait(monkeypatch: pytest.MonkeyPatch) -> list[float]:
    sleeps: list[float] = []
    monkeypatch.setattr(gnomad.random, "uniform", lambda _start, _end: 0.0)
    monkeypatch.setattr(gnomad.time, "sleep", sleeps.append)
    return sleeps


def test_transient_timeout_is_retried(monkeypatch: pytest.MonkeyPatch) -> None:
    responses = iter([TimeoutError("timed out"), SUCCESS_RESPONSE])

    def execute(_query: str, _variables: dict) -> dict:
        response = next(responses)
        if isinstance(response, Exception):
            raise response
        return response

    monkeypatch.setattr(gnomad, "execute_graphql", execute)
    sleeps = disable_retry_wait(monkeypatch)

    variants = gnomad.fetch_region_variants_recursive("1", 100, 100)

    assert variants == SUCCESS_RESPONSE["data"]["region"]["variants"]
    assert sleeps == [5.0]


def test_transient_http_error_is_retried(monkeypatch: pytest.MonkeyPatch) -> None:
    responses = iter(
        [
            HTTPError("https://example.test", 503, "unavailable", None, None),
            SUCCESS_RESPONSE,
        ]
    )

    def execute(_query: str, _variables: dict) -> dict:
        response = next(responses)
        if isinstance(response, Exception):
            raise response
        return response

    monkeypatch.setattr(gnomad, "execute_graphql", execute)
    sleeps = disable_retry_wait(monkeypatch)

    variants = gnomad.fetch_region_variants_recursive("1", 100, 100)

    assert variants == SUCCESS_RESPONSE["data"]["region"]["variants"]
    assert sleeps == [5.0]


def test_non_retryable_http_error_fails_immediately(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = 0

    def execute(_query: str, _variables: dict) -> dict:
        nonlocal calls
        calls += 1
        raise HTTPError("https://example.test", 400, "bad request", None, None)

    monkeypatch.setattr(gnomad, "execute_graphql", execute)
    sleeps = disable_retry_wait(monkeypatch)

    with pytest.raises(HTTPError, match="400"):
        gnomad.fetch_region_variants_recursive("1", 100, 100)

    assert calls == 1
    assert sleeps == []


def test_transient_error_uses_bounded_exponential_backoff(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = 0

    def execute(_query: str, _variables: dict) -> dict:
        nonlocal calls
        calls += 1
        raise URLError("connection reset")

    monkeypatch.setattr(gnomad, "execute_graphql", execute)
    sleeps = disable_retry_wait(monkeypatch)

    with pytest.raises(URLError, match="connection reset"):
        gnomad.fetch_region_variants_recursive("1", 100, 100)

    assert calls == 10
    assert sleeps == [5.0, 10.0, 20.0, 40.0, 60.0, 60.0, 60.0, 60.0, 60.0]
