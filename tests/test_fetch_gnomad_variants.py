from __future__ import annotations

import io
import json
from http.client import IncompleteRead, RemoteDisconnected
from urllib.error import HTTPError, URLError

import pytest

from genomics import gnomad


@pytest.fixture(autouse=True)
def no_network_pacing(monkeypatch):
    monkeypatch.setattr(gnomad, "wait_for_request_slot", lambda: None)


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


@pytest.mark.parametrize("response", [{}, {"data": {"region": None}}, {"data": {"region": {"variants": None}}}])
def test_incomplete_success_response_is_not_absence(monkeypatch, response):
    monkeypatch.setattr(gnomad, "execute_graphql", lambda *args: response)
    with pytest.raises(RuntimeError, match="Incomplete gnomAD response"):
        gnomad.fetch_region_variants_recursive("1", 1, 25000)


def test_invalid_json_response_is_pending_transport_evidence(monkeypatch):
    monkeypatch.setattr(gnomad.urllib.request, "urlopen", lambda *args, **kwargs: io.BytesIO(b"<html>unavailable</html>"))
    with pytest.raises(RuntimeError, match="Invalid JSON"):
        gnomad.execute_graphql("query", {})


def disable_retry_wait(monkeypatch: pytest.MonkeyPatch) -> list[float]:
    sleeps: list[float] = []
    monkeypatch.setattr(gnomad.random, "uniform", lambda _start, _end: 0.0)
    monkeypatch.setattr(gnomad.time, "sleep", sleeps.append)
    return sleeps


@pytest.mark.parametrize("error", [
    TimeoutError("timed out"),
    IncompleteRead(b"partial", 100),
    RemoteDisconnected("Remote end closed connection without response"),
])
def test_transient_transport_error_is_retried(
    monkeypatch: pytest.MonkeyPatch, error: Exception,
) -> None:
    responses = iter([error, SUCCESS_RESPONSE])

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


@pytest.mark.parametrize("status", [200, 400, 500])
def test_region_limit_response_splits_without_retrying_large_request(
    monkeypatch: pytest.MonkeyPatch, status: int,
) -> None:
    calls = []

    def urlopen(request, timeout):
        variables = json.loads(request.data)["variables"]
        start, stop = variables["start"], variables["stop"]
        calls.append((start, stop))
        if (start, stop) == (1, 1000):
            body = json.dumps({"errors": [{"message":
                "This region has too many variants to display. "
                "Select a smaller region to view variants."
            }], "data": None}).encode()
            if status != 200:
                raise HTTPError(request.full_url, status, "region limit", None, io.BytesIO(body))
            return io.BytesIO(body)
        return io.BytesIO(json.dumps({"data": {"region": {"variants": [
            {"chrom": "1", "pos": start, "ref": "A", "alt": "G"}
        ]}}}).encode())

    monkeypatch.setattr(gnomad.urllib.request, "urlopen", urlopen)
    sleeps = disable_retry_wait(monkeypatch)

    variants = gnomad.fetch_region_variants_recursive("1", 1, 1000)

    assert [v["pos"] for v in variants] == [1, 501]
    assert calls == [(1, 1000), (1, 500), (501, 1000)]
    assert sleeps == []


def test_http_failure_retains_bounded_provider_detail(monkeypatch: pytest.MonkeyPatch) -> None:
    def urlopen(request, timeout):
        body = b'{"errors":[{"message":"Unknown dataset"}]}' + b"x" * 10_000
        raise HTTPError(request.full_url, 400, "Bad Request", None, io.BytesIO(body))

    monkeypatch.setattr(gnomad.urllib.request, "urlopen", urlopen)
    sleeps = disable_retry_wait(monkeypatch)

    with pytest.raises(HTTPError, match="Unknown dataset") as failure:
        gnomad.fetch_region_variants_recursive("1", 1, 1000)

    assert len(str(failure.value)) < 4200
    assert failure.value.code == 400
    assert sleeps == []


def test_region_limit_at_minimum_window_remains_a_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    def execute(_query, _variables):
        return {"errors": [{"message": "Select a smaller region"}], "data": None}

    monkeypatch.setattr(gnomad, "execute_graphql", execute)
    sleeps = disable_retry_wait(monkeypatch)

    with pytest.raises(RuntimeError, match="smaller region"):
        gnomad.fetch_region_variants_recursive("1", 1, gnomad.GNOMAD_REGION_MIN_WINDOW_BP)

    assert sleeps == []


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
